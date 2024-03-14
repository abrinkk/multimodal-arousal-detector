# ar_network.py
# Andreas Brink-Kjaer
# Spring 2018
#
# Based on scripts by Caspar Aleksander Bang Jespersen
#
''' 
ar_network generates the network graph using models.py and ar_resnet.py.
'''

import ar_config
#import ar_resnet
import tensorflow as tf
import numpy as np
import models

class ARModel(object):
    def __init__(self, config):
        '''Creates a new network using the set configurations.

        Args:
            config: configurations from ar_config.py.
        '''
        

        self.is_training = config.is_training

        # Placeholders
        self._features = tf.compat.v1.placeholder(tf.float32, [None, config.num_features, 1], name='ModelInput')
        self._targets = tf.compat.v1.placeholder(tf.float32, [None, config.num_classes], name='ModelOutput')
        self._mask = tf.compat.v1.placeholder(tf.float32, [None, int(config.num_classes/2)], name='ModelMask')
        self._batch_size = tf.compat.v1.placeholder(tf.int32, name='BatchSize')
        self._step = tf.compat.v1.placeholder(tf.float32, [])
        

        batch_size_int = tf.reshape(self._batch_size, []) 

        #self._initial_state = tf.placeholder_with_default(tf.zeros([batch_size_int,  ,config.num_hidden*2],dtype=tf.float32), [None, None, config.num_hidden*2],name='InitialState')

        # Reshape input data
        with tf.compat.v1.variable_scope('input_hidden') as scope:
            inputs = self._features
            eeg = tf.expand_dims(inputs[:,:128,:],3)
            eog = tf.reshape(inputs[:,128:384,:],shape=[batch_size_int,128,2,1])
            emg = tf.expand_dims(inputs[:,384:,:],3)

        # Create resnet cnn structure for each modality
        with tf.compat.v1.variable_scope('eeg_resnet') as scope:
            print('eeg resnet: ')
            hidden_eeg = models.resnet(eeg,config.num_resnet,1)
        with tf.compat.v1.variable_scope('eog_resnet') as scope:
            print('eog resnet: ')
            hidden_eog = models.resnet(eog,config.num_resnet,2)
        with tf.compat.v1.variable_scope('emg_resnet') as scope:
            print('emg resnet: ')
            hidden_emg = models.resnet(emg,config.num_resnet,1)

        # Concatenate output features for each resnet structure
        hidden_concat = tf.concat([hidden_eeg,hidden_eog,hidden_emg],1)
        hidden_concat = tf.expand_dims(hidden_concat,0)

        nHid = hidden_concat.get_shape()
        print('ResNet concat shape: ',nHid)
       
        # Regularization
        if config.is_training and config.keep_prob < 1.0:
            iKeepProb = config.keep_prob
            oKeepProb = config.keep_prob
        else:
            iKeepProb = 1
            oKeepProb = 1

        # Bi-directional LSTM network
        with tf.compat.v1.variable_scope('blstm') as scope:
            # Forward cells
            cell_fw = tf.compat.v1.nn.rnn_cell.BasicLSTMCell(config.num_hidden, forget_bias=1.0,state_is_tuple=True)
            cell_fw = tf.compat.v1.nn.rnn_cell.DropoutWrapper(cell_fw, input_keep_prob=iKeepProb, output_keep_prob=oKeepProb)
            # Backward cells
            cell_bw = tf.compat.v1.nn.rnn_cell.BasicLSTMCell(config.num_hidden, forget_bias=1.0,state_is_tuple=True)
            cell_bw = tf.compat.v1.nn.rnn_cell.DropoutWrapper(cell_bw, input_keep_prob=iKeepProb, output_keep_prob=oKeepProb)
            #initial_state = tf.nn.rnn_cell.LSTMStateTuple(self._initial_state[:,:config.num_hidden],self._initial_state[:,config.num_hidden:])
            outputs,final_state = tf.compat.v1.nn.bidirectional_dynamic_rnn(cell_fw, cell_bw, hidden_concat, dtype=tf.float32)#, initial_state = initial_state)
            print('BLSTM OUTPUT; ',outputs)
            outputs = tf.concat(outputs,2)
            print('BLSTM CONCAT: ',outputs)
            outputs = tf.reshape(outputs, [-1,2*config.num_hidden])
            print('BLSTM RESHAPE: ',outputs.shape)

        # Fully connected layer 1 and 2
        with tf.compat.v1.variable_scope('hidden_reduce') as scope:
            hidden_r1 = models.fc_layer(outputs,[2*config.num_hidden, config.num_hidden])
            hidden_r2 = models.fc_layer(hidden_r1,[config.num_hidden, int(config.num_hidden/2)])

        # Fully connected output layer
        with tf.compat.v1.variable_scope('hidden_output') as scope:
            logits = models.output(hidden_r2,[int(config.num_hidden/2),config.num_classes])

        # Evaluate
        cross_ent, TP_ar, FP_ar, FN_ar, TP_w, FP_w, FN_w = self.intelligent_cost(logits)
        self._loss = cross_ent
        self._learning_rate = config.learning_rate
        self._TP_ar = TP_ar
        self._FP_ar = FP_ar
        self._FN_ar = FN_ar
        self._TP_w = TP_w
        self._FP_w = FP_w
        self._FN_w = FN_w
        self._logits = logits
        self._cross_ent = cross_ent
        self._softmax_ar = tf.nn.softmax(logits[:,:2])
        self._softmax_w = tf.nn.softmax(logits[:,2:])
        self._softmax = tf.concat((self._softmax_ar,self._softmax_w),1)
        self._correct_w = tf.equal(tf.argmax(logits[:,2:], 1), tf.argmax(self._targets[:,2:], 1))
        self._accuracy_w = tf.reduce_mean(tf.cast(self._correct_w, tf.float32))

        tf.compat.v1.summary.scalar('Cross_ent',self._cross_ent)
        tf.compat.v1.summary.scalar('Learning Rate',self._learning_rate)
        tf.compat.v1.summary.scalar('Acuracy_w',self._accuracy_w)

        # If not training, do not optimize weights
        if not config.is_training:
            return

        # Optimize
        update_ops = tf.compat.v1.get_collection(tf.compat.v1.GraphKeys.UPDATE_OPS)
        with tf.control_dependencies(update_ops):
            optimizer = tf.compat.v1.train.AdamOptimizer(learning_rate = self._learning_rate)

        self._train_op = optimizer.minimize(self._cross_ent)

        # Summary tensorboard
        # Variable histograms
        vars = tf.compat.v1.trainable_variables()
        for v in vars:
            tf.compat.v1.summary.histogram(v.name, v)
        self._summ = tf.compat.v1.summary.merge_all()

    def intelligent_cost(self, logits):
        '''This function calculates the loss function and summary statistics.

        Args: 
            logits: Network output logits.
        '''
        # Clip logits to avoid extreme values
        logits = tf.clip_by_value(logits,-1e10,1e+10)
        # Evalute loss for arousals and wake
        cross_ent_ar = tf.nn.softmax_cross_entropy_with_logits(logits=logits[:,0:2], labels=self._targets[:,0:2])
        cross_ent_w  = tf.nn.softmax_cross_entropy_with_logits(logits=logits[:,2:], labels=self._targets[:,2:])
        cross_ent = tf.reduce_sum(tf.multiply(cross_ent_ar,self._mask[:,0])) / tf.reduce_sum(self._mask[:,0]) + tf.reduce_sum(tf.multiply(cross_ent_w,self._mask[:,1])) / tf.reduce_sum(self._mask[:,1])
        # Calculate TP, FP, FN
        pred_ar = tf.argmax(logits[:,0:2],1)
        target_ar = tf.argmax(self._targets[:,0:2],1)
        TP_ar = tf.reduce_sum(tf.cast(tf.equal(tf.boolean_mask(pred_ar,tf.equal(target_ar,1)),1),tf.float32))
        FP_ar = tf.reduce_sum(tf.cast(tf.equal(tf.boolean_mask(pred_ar,tf.equal(target_ar,0)),1),tf.float32))
        FN_ar = tf.reduce_sum(tf.cast(tf.equal(tf.boolean_mask(pred_ar,tf.equal(target_ar,1)),0),tf.float32))
        # Calculate TP, FP, FN
        pred_w = tf.argmax(logits[:,2:],1)
        target_w = tf.argmax(self._targets[:,2:],1)
        TP_w = tf.reduce_sum(tf.cast(tf.equal(tf.boolean_mask(pred_w,tf.equal(target_w,1)),1),tf.float32))
        FP_w = tf.reduce_sum(tf.cast(tf.equal(tf.boolean_mask(pred_w,tf.equal(target_w,0)),1),tf.float32))
        FN_w = tf.reduce_sum(tf.cast(tf.equal(tf.boolean_mask(pred_w,tf.equal(target_w,1)),0),tf.float32))

        return cross_ent, TP_ar, FP_ar, FN_ar, TP_w, FP_w, FN_w
        
    def gather_loss(self):
        '''Obsolete function'''
        loss_averages = tf.train.ExponentialMovingAverage(0.9,name='avg_loss')
        losses = tf.compat.v1.get_collection('losses')
        total_loss = tf.add_n(losses, name='total_loss')
        loss_averages_op = loss_averages.apply(losses + [total_loss])

        return total_loss

    @property
    def features(self):
        return self._features

    @property   
    def final_state(self):
        return self._final_state

    @property
    def initial_state(self):
        return self._initial_state
    @property
    def targets(self):
        return self._targets
        
    @property
    def mask(self):
        return self._mask
        
    @property
    def batch_size(self):
        return self._batch_size

    @property
    def learning_rate(self):
        return self._learning_rate

    @property
    def loss(self):
        return self._loss
        
    @property
    def cross_ent(self):
        return self._cross_ent

    @property
    def correct_w(self):
        return self._correct_w 

    @property
    def accuracy(self):
        return self._accuracy_w
     
    @property
    def baseline(self):
        return self._baseline

    @property
    def train_op(self):
        return self._train_op

    @property
    def predict(self):
        return self._predict

    @property
    def logits(self):
        return self._logits
        
    @property
    def confidence(self):
        return self._confidence    

    @property
    def ar_prob(self):
        return self._ar_prob

    @property
    def softmax(self):
        return self._softmax

    @property
    def softmax_ar(self):
        return self._softmax_ar

    @property
    def softmax_w(self):
        return self._softmax_w

    @property
    def ckappa(self):
        return self._ckappa

    @property
    def TP_ar(self):
        return self._TP_ar

    @property
    def FP_ar(self):
        return self._FP_ar

    @property
    def FN_ar(self):
        return self._FN_ar

    @property
    def TP_w(self):
        return self._TP_w

    @property
    def FP_w(self):
        return self._FP_w

    @property
    def FN_w(self):
        return self._FN_w

    @property
    def precision(self):
        return self._precision

    @property
    def recall(self):
        return self._recall

    @property
    def F1(self):
        return self._F1

    @property
    def summ(self):
        return self._summ

    @property
    def step(self):
        return self._step