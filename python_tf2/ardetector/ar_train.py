# ar_train.py
# Andreas Brink-Kjaer
# Spring 2018
#
# Based on scripts by Caspar Aleksander Bang Jespersen
#
''' 
ar_train iterates over all training data using ar_reader.py and ar_config.py for data loading and 
trains the neural network graph, which is created in ar_network.py, models.py, 
and ar_resnet.py. The performance is validated on 10 PSG's every 5000 batches and on 30 every 50000 iterations.
'''

from datetime import datetime
import time
import numpy as np
import tensorflow as tf
from tensorflow.python.platform import gfile
from tensorflow.core.framework import summary_pb2
import ar_config
import ar_reader
import ar_network
import ar_perf

flags = tf.compat.v1.app.flags
FLAGS = flags.FLAGS
flags.DEFINE_string('model', 'resnet', 'Which model to build')
flags.DEFINE_integer('proceed', 200000, 'Whether or not to continue previous')
flags.DEFINE_float('lr',0.001,'Learning rate')
flags.DEFINE_integer('batch_size',5*60,'Batch size')
flags.DEFINE_float('kp',1.0,'Dropout keep probability')
flags.DEFINE_integer('NH',128,'Number of hidden units')
flags.DEFINE_integer('RL',20,'Number of ResNet layers')

LOGDIR = "/tmp/train_stats/W_NORM_SW_test_c9_"

def rand_hparam(randomgen=False):
    '''Generates random hyper-parameters or uses hyper-parameter values specified in FLAGS.

    Args: 
        randomgen: boolean variable to use random or predefined hyper-parameters.
    '''
    if randomgen:
        lr = 10**np.random.uniform(-1,-5)
        kp = np.random.choice(np.array([0.5, 0.6, 0.7, 0.8, 0.9, 1.0]))
        bs = 60*np.random.choice(np.array([5, 10, 20, 60]))
        NH = np.random.choice(np.array([64, 128, 192, 256]))
        RL = np.random.choice(np.array([20, 32, 44]))
    else:
        lr = FLAGS.lr
        kp = FLAGS.kp
        bs = FLAGS.batch_size
        NH = FLAGS.NH
        RL = FLAGS.RL
    return lr, kp, bs, NH, RL

def make_hparam_string(lr,kp,bs,NH,RL):
    '''Prints a string that describes used hyper-parameters.

    Args: 
        lr: learning rate.
        kp: keep probability for dropout.
        bs: batch size.
        NH: number of hidden units.
        RL: Resnet layers.
    '''
    return "lr_%.7f,kp_%.2f,bs_%.0f,nh_%.0f,rl_%.0f" % (lr, kp, bs, NH, RL)

def main(argv=None):
    '''The main function sets configurations and initiates the training'''

    # Set hyper-parameters
    lr, kp, bs, NH, RL = rand_hparam(randomgen = False)

    # Configurations for test and validation data
    tf.compat.v1.reset_default_graph()
    config_train = ar_config.ARConfig(num_hidden = NH, lr = lr, kp = kp, batch_size = bs, resnet_size = RL, model_name='resnet', is_training=True)
    config_val = ar_config.ARConfig(num_hidden = NH, lr = lr, kp = kp, batch_size = bs, resnet_size = RL, model_name='resnet', is_training=False)

    # Creates a folder
    if gfile.Exists(config_train.train_model):
        gfile.DeleteRecursively(config_train.train_model)
        gfile.MakeDirs(config_train.train_model)

    # Start model training
    train(config_train,config_val,make_hparam_string(lr,kp,bs,NH,RL))

def train(config,config_val,hparam):
    '''Trains model using data specified in config, validates on data specified in config_val.

    Args: 
        config: configuartions for network and training data.
        config_val: configurations for network ad validation data.
        hparam: string of hyper-parameters generated by make_hparam_string.
    '''

    # Read data on CPU to avoid memory issues on GPU
    with tf.device('/cpu:0'):
        train_data = ar_reader.ArousalData(config.train_dir, config, num_steps=config.max_steps)
        val_data = ar_reader.ArousalData(config_val.val_dir, config_val, num_steps=config.max_steps)

    # Create network graph
    sess_config = tf.compat.v1.ConfigProto(log_device_placement=False)
    sess_config.gpu_options.allow_growth = True
    with tf.Graph().as_default(), tf.compat.v1.Session(config=sess_config) as session:
        with tf.compat.v1.variable_scope('model', reuse=False) as scope:
            m = ar_network.ARModel(config)
            session.run(tf.compat.v1.local_variables_initializer())
            glob_vars = tf.compat.v1.global_variables()
            model_vars = [k for k in glob_vars if k.name.startswith(scope.name)]
            s = tf.compat.v1.train.Saver(model_vars, max_to_keep=None)

        # Load trained model parameters or random initialization
        if FLAGS.proceed > 0:
            s.restore(session, config.checkpoint_file(FLAGS.proceed))
            train_data.iter_steps += FLAGS.proceed
        else:
            tf.compat.v1.global_variables_initializer().run()

        # Tensorboard writer
        writer = tf.compat.v1.summary.FileWriter(LOGDIR + hparam)
        writer.add_graph(session.graph)

        # Variables for hold metric averages
        loss_avg = 0
        tp_l_ar = [0]*500
        fp_l_ar = [0]*500
        fn_l_ar = [0]*500
        tp_l_w = [0]*500
        fp_l_w = [0]*500
        fn_l_w = [0]*500
        
        # Iterate training data
        for batch_input, batch_target, batch_mask in train_data:

            # Record batch iteration time
            step = train_data.iter_steps + 1
            start_time = time.time()

            #print(train_data.batch_size)
            #print(batch_input.shape)
            #print(batch_target.shape)
            
            # Session output
            operations = [m.train_op, m.loss, m.softmax, m.summ, m.TP_ar, m.FP_ar, m.FN_ar, m.TP_w, m.FP_w, m.FN_w]
            # Session input
            params = {
                m.features: batch_input,
                m.targets: batch_target,
                m.mask: batch_mask,
                m.batch_size: train_data.batch_size,
                m.step: train_data.iter_steps
            }
            # Run session
            _, loss, softmax, summa, TP_ar, FP_ar, FN_ar, TP_w, FP_w, FN_w = session.run(operations, feed_dict=params)

            # update performance metrics
            # loss
            loss_avg = 0.9*loss_avg + (1-0.9)*loss
            # ar
            tp_l_ar.insert(0,TP_ar)
            tp_l_ar.pop()
            fp_l_ar.insert(0,FP_ar)
            fp_l_ar.pop()
            fn_l_ar.insert(0,FN_ar)
            fn_l_ar.pop()
            # w
            tp_l_w.insert(0,TP_w)
            tp_l_w.pop()
            fp_l_w.insert(0,FP_w)
            fp_l_w.pop()
            fn_l_w.insert(0,FN_w)
            fn_l_w.pop()
            # batch iteration time
            batch_time = time.time() - start_time

            # Print performance and model validation
            if not step % 100:
                now = datetime.now().time()
                writer.add_summary(summa, step)
                if (sum(tp_l_ar) + sum(fp_l_ar)) > 0:
                    precision_ar = sum(tp_l_ar)/(sum(tp_l_ar) + sum(fp_l_ar))
                else:
                    precision_ar = 0
                if (sum(tp_l_ar) + sum(fn_l_ar)) > 0:
                    recall_ar = sum(tp_l_ar)/(sum(tp_l_ar) + sum(fn_l_ar))
                else:
                    recall_ar = 0
                if (precision_ar + recall_ar) > 0:
                    F1_ar = 2*precision_ar*recall_ar/(precision_ar + recall_ar)
                else:
                    F1_ar = 0
                if (sum(tp_l_w) + sum(fp_l_w)) > 0:
                    precision_w = sum(tp_l_w)/(sum(tp_l_w) + sum(fp_l_w))
                else:
                    precision_w = 0
                if (sum(tp_l_w) + sum(fn_l_w)) > 0:
                    recall_w = sum(tp_l_w)/(sum(tp_l_w) + sum(fn_l_w))
                else:
                    recall_w = 0
                if (precision_w + recall_w) > 0:
                    F1_w = 2*precision_w*recall_w/(precision_w + recall_w)
                else:
                    F1_w = 0
                print('%s | Step %4d | Cross Ent %.5f | Loss_avg %.5f | Ar %.0f | sum(softmax) %.5f | W %.0f | sum(softmax) %.5f | F1_ar %.2f | F1_w %.2f ' % (now, step, loss, loss_avg, np.sum(batch_target[:,1],0),np.sum(softmax[:,1],0), np.sum(batch_target[:,-1],0),np.sum(softmax[:,-1],0),F1_ar,F1_w))
                # Add summary stats
                # ckappa
                #summa_test = summary_pb2.Summary(value=[value])
                #writer.add_summary(summa_test, step)
                # Cost
                value = summary_pb2.Summary.Value(tag="train_loss", simple_value = loss_avg)
                summa_test = summary_pb2.Summary(value=[value])
                writer.add_summary(summa_test, step)
                # recall_ar
                value = summary_pb2.Summary.Value(tag="train_Recall_ar", simple_value = recall_ar)
                summa_test = summary_pb2.Summary(value=[value])
                writer.add_summary(summa_test, step)
                # precision_ar
                value = summary_pb2.Summary.Value(tag="train_precision_ar", simple_value = precision_ar)
                summa_test = summary_pb2.Summary(value=[value])
                writer.add_summary(summa_test, step)
                # F1_ar
                value = summary_pb2.Summary.Value(tag="train_F1_ar", simple_value = F1_ar)
                summa_test = summary_pb2.Summary(value=[value])
                writer.add_summary(summa_test, step)
                # recall_w
                value = summary_pb2.Summary.Value(tag="train_Recall_w", simple_value = recall_w)
                summa_test = summary_pb2.Summary(value=[value])
                writer.add_summary(summa_test, step)
                # precision_w
                value = summary_pb2.Summary.Value(tag="train_precision_w", simple_value = precision_w)
                summa_test = summary_pb2.Summary(value=[value])
                writer.add_summary(summa_test, step)
                # F1_w
                value = summary_pb2.Summary.Value(tag="train_F1_w", simple_value = F1_w)
                summa_test = summary_pb2.Summary(value=[value])
                writer.add_summary(summa_test, step)
                # Iteration time
                value = summary_pb2.Summary.Value(tag="Iter_time", simple_value = batch_time)
                summa_test = summary_pb2.Summary(value=[value])
                writer.add_summary(summa_test, step)

                # Validation
                if not step % 50000:
                    n_val_files = 30
                    val_data.iter_rewind = -1 + 10
                else:
                    n_val_files = 10
                    val_data.iter_rewind = -1
                if not step % 5000:
                    # Save checkpoint
                    s.save(session, config.checkpoint_file(), global_step=step)
                    # Validation
                    cost = 0
                    n_batches = 0
                    TP_ar = 0
                    FP_ar = 0
                    FN_ar = 0
                    TPe_ar = np.full(shape=(39), fill_value=0, dtype=np.float32)
                    FPe_ar = np.full(shape=(39), fill_value=0, dtype=np.float32)
                    FNe_ar = np.full(shape=(39), fill_value=0, dtype=np.float32)
                    TP_w = 0
                    FP_w = 0
                    FN_w = 0
                    for batch_input, batch_target, batch_mask in val_data:
                        if val_data.iter_rewind > n_val_files:
                            break
                        # Session output
                        operations = [m.loss, m.softmax, m.TP_ar, m.FP_ar, m.FN_ar, m.TP_w, m.FP_w, m.FN_w]
                        # Session input
                        params = {
                            m.features: batch_input,
                            m.targets: batch_target,
                            m.mask: batch_mask,
                            m.batch_size: val_data.batch_size
                        }
                        # Run session
                        cross_ent, softmax, TP_ar_, FP_ar_, FN_ar_, TP_w_, FP_w_, FN_w_ = session.run(operations, feed_dict=params)
                        # Calculate sample and event TP, FP, and FN
                        TPe_ar_, FPe_ar_, FNe_ar_ = ar_perf.CalcHits(softmax[:,1],batch_target[:,1])
                        cost += cross_ent
                        TP_ar += TP_ar_
                        FP_ar += FP_ar_
                        FN_ar += FN_ar_
                        TPe_ar += TPe_ar_
                        FPe_ar += FPe_ar_
                        FNe_ar += FNe_ar_
                        TP_w += TP_w_
                        FP_w += FP_w_
                        FN_w += FN_w_
                        n_batches += 1
                    #cost = cost/(val_data.test_split*n_val_files)
                    cost = cost / n_batches
                    precision = TP_ar/(TP_ar + FP_ar)
                    recall = TP_ar/(TP_ar + FN_ar)
                    F1 = 2*precision*recall/(precision + recall)
                    precision_e, recall_e, F1_e, threshold = ar_perf.CalcPerf(TPe_ar,FPe_ar,FNe_ar)
                    precision_w = TP_w/(TP_w + FP_w)
                    recall_w = TP_w/(TP_w + FN_w)
                    F1_w = 2*precision_w*recall_w/(precision_w + recall_w)
                    # Add test summary statistics
                    # Cost
                    value = summary_pb2.Summary.Value(tag="Test_cost", simple_value = cost)
                    summa_test = summary_pb2.Summary(value=[value])
                    writer.add_summary(summa_test, step)
                    # Precision
                    value = summary_pb2.Summary.Value(tag="Test_sample_precision", simple_value = precision)
                    summa_test = summary_pb2.Summary(value=[value])
                    writer.add_summary(summa_test, step)
                    # Recall
                    value = summary_pb2.Summary.Value(tag="Test_sample_recall", simple_value = recall)
                    summa_test = summary_pb2.Summary(value=[value])
                    writer.add_summary(summa_test, step)
                    # F1
                    value = summary_pb2.Summary.Value(tag="Test_sample_F1_score", simple_value = F1)
                    summa_test = summary_pb2.Summary(value=[value])
                    writer.add_summary(summa_test, step)
                    # Precision Wake
                    value = summary_pb2.Summary.Value(tag="Test_sample_precision_wake", simple_value = precision_w)
                    summa_test = summary_pb2.Summary(value=[value])
                    writer.add_summary(summa_test, step)
                    # Recall Wake
                    value = summary_pb2.Summary.Value(tag="Test_sample_recall_wake", simple_value = recall_w)
                    summa_test = summary_pb2.Summary(value=[value])
                    writer.add_summary(summa_test, step)
                    # F1 Wake
                    value = summary_pb2.Summary.Value(tag="Test_sample_F1_score_wake", simple_value = F1_w)
                    summa_test = summary_pb2.Summary(value=[value])
                    writer.add_summary(summa_test, step)
                    # Precision event
                    value = summary_pb2.Summary.Value(tag="Test_precision", simple_value = precision_e)
                    summa_test = summary_pb2.Summary(value=[value])
                    writer.add_summary(summa_test, step)
                    # Recall event
                    value = summary_pb2.Summary.Value(tag="Test_recall", simple_value = recall_e)
                    summa_test = summary_pb2.Summary(value=[value])
                    writer.add_summary(summa_test, step)
                    # F1 event
                    value = summary_pb2.Summary.Value(tag="Test_F1_score", simple_value = F1_e)
                    summa_test = summary_pb2.Summary(value=[value])
                    writer.add_summary(summa_test, step)
                    # Threshold
                    value = summary_pb2.Summary.Value(tag="Test_threshold", simple_value = threshold)
                    summa_test = summary_pb2.Summary(value=[value])
                    writer.add_summary(summa_test, step)
                    
if __name__ == '__main__':
    tf.compat.v1.app.run()
