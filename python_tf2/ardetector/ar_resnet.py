# ar_resnet.py
# Andreas Brink-Kjaer
# Spring 2018
# Modified from 'https://github.com/xuyuwei/resnet-tf'
#
'''
This script creates the resnet blocks..
'''

import numpy as np
import tensorflow as tf


def weight_variable(shape, name=None):
    '''This functiton creates and specifies initialization of weights.

    Args:
        shape: weight shape.
        name: variable name.
    '''

    initial = tf.random.truncated_normal(shape, stddev=np.sqrt(2/shape[0]))
    #initial = tf.truncated_normal(shape, stddev=0.1)
    return tf.Variable(initial, name=name)

def fc_layer_relu(inpt, shape):
    '''This function creates a fully connected layer with batch normalization and a ReLU activation.

    Args:
        inpt: input data.
        shape: weight shape.
    '''

    fc_w = weight_variable(shape)
    #fc_b = tf.Variable(tf.zeros([shape[1]]))

    fc_z = tf.matmul(inpt, fc_w)# + fc_b
    # Batch norm
    mean, var = tf.nn.moments(fc_z, axes=[0])
    beta = tf.Variable(tf.zeros([shape[1]]), name="beta")
    gamma = tf.Variable(tf.ones([shape[1]]), name="gamma")
    #gamma = 1 + weight_variable([shape[1]], name="gamma")
    #gamma = tf.Variable(tf.truncated_normal([shape[1]],mean=1.0,stddev=0.1),name="gamma")

    batch_norm = tf.nn.batch_normalization(
        fc_z, mean, var, beta, gamma, 0.0001)


    fc_h = tf.nn.relu(batch_norm)
    #fc_h = tf.nn.relu(fc_z)

    return fc_h

def conv_layer(inpt, filter_shape, stride1=1, stride2=1):
    '''This function creates a convolutional layer with batch normalization and a ReLU activation.

    Args:
        inpt: input data.
        filter_shape: shape of convolutional filter.
        stride1: width stride for convolution.
        stride2: channel height stride for convolution.
    '''
    out_channels = filter_shape[3]

    filter_ = weight_variable(filter_shape)
    conv = tf.nn.conv2d(inpt, filters=filter_, strides=[1, stride1, stride2, 1], padding="SAME")
    mean, var = tf.nn.moments(conv, axes=[0,1,2])
    beta = tf.Variable(tf.zeros([out_channels]), name="beta")
    #gamma = tf.Variable(tf.ones([out_channels]), name="gamma")
    gamma = weight_variable([out_channels], name="gamma")
    #gamma = tf.Variable(tf.truncated_normal([out_channels],mean=1.0,stddev=0.1),name="gamma")

    batch_norm = tf.nn.batch_normalization(
        conv, mean, var, beta, gamma, 0.0001)
    out = tf.nn.relu(batch_norm)

    return out

def basic_conv_layer(inpt, filter_shape, stride1=1, stride2=1):
    '''This function creates a simple convolutional layer.

    Args:
        inpt: input data.
        filter_shape: shape of convolutional filter.
        stride1: width stride for convolution.
        stride2: channel height stride for convolution.
    '''
    filter_ = weight_variable(filter_shape)
    conv = tf.nn.conv2d(inpt, filters=filter_, strides=[1, stride1, stride2, 1], padding="SAME")
    return conv

def resnet_batch_norm(conv, out_channels):
    '''This function performs batch normalization for a convolutional layer.

    Args:
        conv: input data.
        out_channels: number of output channels.
    '''
    mean, var = tf.nn.moments(conv, axes=[0,1,2])
    beta = tf.Variable(tf.zeros([out_channels]), name="beta")
    #gamma = tf.Variable(tf.ones([out_channels]), name="gamma")
    gamma = weight_variable([out_channels], name="gamma")
    #gamma = tf.Variable(tf.truncated_normal([out_channels],mean=1.0,stddev=0.1),name="gamma")

    batch_norm = tf.nn.batch_normalization(
        conv, mean, var, beta, gamma, 0.0001)
    return batch_norm

def residual_block(inpt, output_depth, down_sample, projection=False):
    '''This function creates a residual block.

    Args:
        inpt: input data.
        output_depth: output depth of block.
        down_sample: option to downsample width.
        projection: option to use projection shortcut or zero-padding.
    '''
    
    input_width = inpt.get_shape().as_list()[2]
    input_depth = inpt.get_shape().as_list()[3]
    if down_sample:
        stride = 2
    else:
        stride = 1

    # Option full pre-activation
    #BN1 = resnet_batch_norm(inpt, input_depth)
    #ReLU1 = tf.nn.relu(BN1)
    #conv1 = basic_conv_layer(ReLU1, [3, 1, input_depth, output_depth], stride)
    #BN2 = resnet_batch_norm(conv1, output_depth)
    #ReLU2 = tf.nn.relu(BN2)
    #out = basic_conv_layer(ReLU2, [3, 1, output_depth, output_depth], 1)
    # Original
    conv1 = basic_conv_layer(inpt, [9, 1, input_depth, output_depth], stride)
    BN1 = resnet_batch_norm(conv1, output_depth)
    ReLU1 = tf.nn.relu(BN1)
    conv2 = basic_conv_layer(ReLU1, [9, 1, output_depth, output_depth])
    out = resnet_batch_norm(conv2, output_depth)

    if input_depth != output_depth:
        if projection:
            # Option B: Projection shortcut
            input_layer = conv_layer(inpt, [1, 1, input_depth, output_depth], 2)
        else:
            # Option A: Zero-padding
            filter_ = [1,2,input_width,1]
            input_layer = tf.nn.max_pool2d(input=inpt, ksize=filter_, strides=filter_, padding='SAME')
            input_layer = tf.pad(input_layer, [[0,0], [0,0], [0,0], [0, output_depth - input_depth]])
    else:
        input_layer = inpt

    res = tf.nn.relu(out + input_layer)
    return res