# ar_predict.py
# Andreas Brink-Kjaer
# Spring 2018
#
# Based on scripts by Caspar Aleksander Bang Jespersen
#
''' 
ar_predict iterates over all data in directory specified by flag 'pathname' using ar_reader.py and ar_config.py for data loading and 
uses a trained model specified by flags 'model' and 'ckpt' to predict arousal and wake. The model predictions are saved in 'output_dir'. 
The 'overwrite' flag can be set to 1 to overwrite predictions in output_dir.
'''

import os
import numpy as np
import tensorflow as tf
import tempfile

import ar_network
import ar_config
import ar_reader

flags = tf.app.flags
FLAGS = flags.FLAGS
flags.DEFINE_string('pathname', 'L:\\LovbeskyttetMapper\\Narkolepsi - Matlab Julie AEC\\Pseudoanonymized data - new\\ProcessedDataAndreasBK2019\\Training', 'Files to execute.')
flags.DEFINE_string('model', 'resnet', 'Which model to use')
flags.DEFINE_integer('ckpt', '350000', 'Which checkpoint to use (blank for none)')#70000
flags.DEFINE_string('output_dir','L:\\LovbeskyttetMapper\\Narkolepsi - Matlab Julie AEC\\Pseudoanonymized data - new\\ProcessedDataAndreasBK2019\\Predictions','Directory for predictions')
flags.DEFINE_integer('overwrite',0,'Overwrite previous predictions (1 = Overwrite)')

def main(argv=None):
    '''The main function sets configurations and runs predict that computes and writes model predictions'''

    # Settings
    config = ar_config.ARConfig(num_hidden = 128, lr = 0.001, kp = 1.0, batch_size = 5*60, resnet_size = 20, model_name='resnet', is_training=False)
    
    _, _ = predict(config,FLAGS.ckpt,FLAGS.pathname)


def predict(config,ckpt,file):
    '''predict uses the config and ckpt to create and load a network graph. The model reads all files in FLAGS.pathname and makes predictions.

    Args: 
        config: configurations for network and data set.
        ckpt: checkpoint number for model weights to load.
        file: Directory for data.
    '''

    # Load data and adjust batch size in CPU
    with tf.device('/cpu:0'):
        data = ar_reader.ArousalData(FLAGS.pathname, config, num_steps=config.max_steps, overwrite = FLAGS.overwrite, output_dir = FLAGS.output_dir)

    # Creates network graph
    sess_config = tf.ConfigProto(log_device_placement=False)
    sess_config.gpu_options.allow_growth = True
    with tf.Graph().as_default(), tf.Session(config=sess_config) as session:
        with tf.variable_scope('model', reuse=False) as scope:
            m = ar_network.ARModel(config)
            glob_vars = tf.global_variables()
            model_vars = [k for k in glob_vars if k.name.startswith(scope.name)]
            s = tf.train.Saver(model_vars, max_to_keep=None)

        # Loads model weights
        s.restore(session, config.checkpoint_file(ckpt))

        current_batch = 0
        for batch_input, batch_target, batch_mask in data:
            #test_step += 1
            current_batch += 1
            #if test_step > data.test_split * 10:
            #    break
            operations = [m.softmax]
            params = {
                m.features: batch_input,
                m.targets: batch_target,
                m.mask: batch_mask,
                m.batch_size: data.batch_size,
                m.step: data.iter_steps
            }
            softmax = session.run(operations, feed_dict=params)
            # Get ar probability
            softmax = softmax[0]
            ar = np.array([x[1] for m, x in enumerate(softmax)])
            # Get W probability
            wake = np.array([x[3] for m, x in enumerate(softmax)])
            # Reshape and save
            if current_batch == 1:
                output_ar = ar
                output_w = wake
            else:
                output_ar = np.append(output_ar,ar,0)
                output_w = np.append(output_w,wake,0)
            if current_batch == data.num_batches:
                # Save output
                current_batch = 0
                filename = data.filename
                file = filename.split('\\')[-1]
                output_file = str.join('\\',[FLAGS.output_dir, file])
                if data.batch_shift == 0:
                    pred_file = open(output_file,'w')
                else:
                    pred_file = open(output_file,'a')
                np.savetxt(pred_file,(output_ar,output_w),delimiter=',',fmt='%.2f')
                pred_file.close()

        return output_ar, output_w


if __name__ == '__main__':
    tf.app.run()
