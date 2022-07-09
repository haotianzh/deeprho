import logging
import tensorflow as tf
from tensorflow.keras import layers, models

def _block1(x, features, kernel_size, pool_size, regularizer=False):
    x = layers.Conv2D(features,
                      kernel_size=kernel_size,
                      padding='same',
                      kernel_initializer='he_normal')(x)
    x = layers.BatchNormalization()(x)
    x = layers.Activation('tanh')(x)
    x = layers.MaxPooling2D(pool_size=pool_size)(x)
    x = layers.Dropout(.1)(x)
    return x


def RecombNet1(num_layers, num_features, input_shape, output_shape):
    inputs = layers.Input(shape=input_shape)
    x = layers.Conv2D(num_features[0],
                      kernel_size=(3,3),
                      padding='same',
                      kernel_initializer='he_normal')(inputs)
    # stack CNN blocks
    for i in range(num_layers-1):
        x = _block1(x, num_features[i+1])
    # feed-forward, this part is ignored when perform transfer learning
    x = layers.Flatten()(x)
    x = layers.Dense(128)(x)
    x = layers.Activation('tanh')(x)
    x = layers.Dense(64)(x)
    x = layers.Activation('tanh')(x)
    x = layers.Dense(output_shape)(x)
    return x


def RecombNet2(num_layers, num_features, input_shape, output_shape):
    inputs = layers.Input(shape=input_shape)
    x = layers.Conv2D(num_features[0],
                      kernel_size=(3,3),
                      padding='same',
                      kernel_initializer='he_normal')(inputs)
    # stack CNN blocks
    for i in range(num_layers-1):
        x = _block1(x, num_features[i+1])
    # feed-forward, this part is ignored when perform transfer learning
    x = layers.Flatten()(x)
    x = layers.Dense(128)(x)
    x = layers.Activation('tanh')(x)
    x = layers.Dense(64)(x)
    x = layers.Activation('tanh')(x)
    x = layers.Dense(output_shape)(x)
    return x