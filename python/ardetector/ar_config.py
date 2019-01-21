# ar_config.py
# Andreas Brink-Kjaer
# Spring 2018
#
# Based on scripts by Caspar Aleksander Bang Jespersen
#
'''
This file defines the config class used to distinguish between the different network types defined.
'''

import os
import numpy as np


class Config(object):

    @staticmethod
    def get(scope, model_name):
        if scope == 'ar':
            return ARConfig(model_name)
        else:
            raise Exception

    def __init__(self, scope, num_features, num_hidden, kp, lr, num_classes, batch_size, resnet_size, wake_def, model_name='resnet', is_train=False):
        '''This function initializes the config class by specifying hyper-parameters and data directories.

        Args: 
            scope: set to 'ar', obsolete.
            num_featres: number of samples in each bin (fs = 128 Hz, bin = 1 s -> num_features = 128).
            num_hidden: number of hidden unit in LSTM and fully connected layers.
            kp: keep probability for dropout.
            lr: learning rate.
            num_classes: number of classes ([non-arousal,arousal],[non-wake,wake]) -> num_classes = 4
            batch_size: batch size for training
            resnet_size: parameter for specifying number of layers in resnet structure.
            wake_def: definition of wake (wake_def = 0 -> wake = W, wake_def = 1 -> wake = {W,N1})
            model_name: model name.
            is_train: boolean option for training / testing.
        '''

        # Model folder
        root_python = os.path.dirname(os.path.realpath(__file__))
        root_base = str.join('\\', root_python.split('\\')[:-2])
        self.model_dir = os.path.join(root_python, 'model', scope, model_name)
        if not os.path.isdir(self.model_dir):
            os.mkdir(self.model_dir)

        # Data
        #data_dir = os.path.join(root_base, 'matlab', 'data')
        data_dir = 'D:\\ardetector\\data'
        self.train_model = os.path.join(self.model_dir,'train')
        self.train_dir = os.path.join(data_dir, 'train')
        self.val_dir = os.path.join(data_dir, 'val')
        self.test_dir = os.path.join(data_dir,'test')
        if not os.path.isdir(self.train_model):
            os.mkdir(self.train_model)

        # Configuration
        self.model_name = model_name
        self.scope = scope
        self.is_training = is_train
        self.num_features = num_features
        self.num_classes = num_classes
        self.batch_size = batch_size
        self.wake_def = wake_def

        if model_name == 'resnet':
            self.type = self.TYPE_RESNET
            self.bidirectional = True
            self.learning_rate = lr
            self.num_hidden = num_hidden
            self.keep_prob = kp
            self.num_resnet = resnet_size
        else:
            raise Exception

        # Training
        self.max_steps = 500000

    def checkpoint_file(self, ckpt=0):
        '''Generates directory name for checkpoints.

        Args: 
            ckpt: checkpoint number.
        '''
        if ckpt == 0:
            return os.path.join(self.model_dir, 'model.ckpt')
        else:
            return os.path.join(self.model_dir, 'model.ckpt-%.0f' % ckpt)

    @property
    def TYPE_RESNET(self):
        return 'RESNET'


class ARConfig(Config):

    def __init__(self, num_hidden = 128, lr = 10**(-4), kp = 0.75, batch_size = 4*60, resnet_size = 32, model_name='resnet', is_training=False):
        '''This function calls the config class to set hyper-parameters and data directories.

        Args: 
            num_hidden: number of hidden unit in LSTM and fully connected layers.
            lr: learning rate.
            kp: keep probability for dropout.
            batch_size: batch size for training
            resnet_size: parameter for specifying number of layers in resnet structure.
            model_name: model name.
            is_training: boolean option for training / testing.
        '''
        scope = 'ar'
        num_features = 4 * 128
        num_classes = 4
        wake_def = 0

        super(ARConfig, self).__init__(scope, num_features, num_hidden, kp, lr, num_classes, batch_size, resnet_size, wake_def, model_name, is_training)