# models.py
# Andreas Brink-Kjaer
# Spring 2018
# Modified from 'https://github.com/xuyuwei/resnet-tf'
#
'''
This script creates full resnet structure using ar_resnet.
'''

import tensorflow as tf
from ar_resnet import fc_layer_relu, conv_layer, residual_block, weight_variable

n_dict = {20:1, 32:2, 44:3, 56:4}
# ResNet architectures used for CIFAR-10
def resnet(inpt, n, nsig):
    '''This function creates the resnet strucutre similar to that used for the CIFAR-10 dataset in the official ResNet 
    publication 'https://arxiv.org/abs/1512.03385'. The function uses the script ar_resnet to create substrucutres.

    Args: 
        inpt: Input data
        n, number specifying the number of layers [20, 32, 44, 56] -> [7, 13, 19, 25]
        nsig, number of signals (1 for EEG and EMG, 2 for EOG)
    '''

    n = int(n)
    if n < 20 or (n - 20) % 12 != 0:
        print("ResNet depth invalid.")
        return

    num_conv = int((n - 20) / 12 + 1)
    layers = []

    with tf.compat.v1.variable_scope('conv1'):
        conv1 = conv_layer(inpt, [9, nsig, 1, 16], 1, nsig)
        layers.append(conv1)

    for i in range (num_conv):
        with tf.compat.v1.variable_scope('conv2_%d' % (i+1)):
            conv2_x = residual_block(layers[-1], 16, False)
            conv2 = residual_block(conv2_x, 16, False)
            layers.append(conv2_x)

            layers.append(conv2)
        #assert conv2.get_shape().as_list()[1:] == [128, 1, 16]
        print(conv2.get_shape().as_list()[1:])

    for i in range (num_conv):
        down_sample = True if i == 0 else False
        with tf.compat.v1.variable_scope('conv3_%d' % (i+1)):
            conv3_x = residual_block(layers[-1], 32, down_sample)
            conv3 = residual_block(conv3_x, 32, False)
            layers.append(conv3_x)
            layers.append(conv3)

        #assert conv3.get_shape().as_list()[1:] == [64, 1, 32]
        print(conv3.get_shape().as_list()[1:])
    
    for i in range (num_conv):
        down_sample = True if i == 0 else False
        with tf.compat.v1.variable_scope('conv4_%d' % (i+1)):
            conv4_x = residual_block(layers[-1], 64, down_sample)
            conv4 = residual_block(conv4_x, 64, False)
            layers.append(conv4_x)
            layers.append(conv4)

        #assert conv4.get_shape().as_list()[1:] == [32, 1, 64]
        print(conv4.get_shape().as_list()[1:])

    with tf.compat.v1.variable_scope('fc'):
        global_pool = tf.reduce_mean(layers[-1], [1, 2])
        #assert global_pool.get_shape().as_list()[1:] == [64]
        print(global_pool.get_shape().as_list()[1:])
        
        layers.append(global_pool)
        #out = softmax_layer(global_pool, [64, 10])
        #layers.append(out)

    return layers[-1]

def output(inpt,shape):
    '''This function creates a fully connected output layer.

    Args:
        inpt: input data.
        shape: weight shape.
    '''

    fc_w = weight_variable(shape)
    fc_b = tf.Variable(tf.zeros([shape[1]]))

    return tf.matmul(inpt, fc_w) + fc_b

def fc_layer(inpt,shape):
    '''Fully connected layer with ReLU activation

    Args:
        inpt: input data.
        shape: weight shape.
    '''

    return fc_layer_relu(inpt,shape)