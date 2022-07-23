"""

This script is used for transfer learning when user seeks to retrain the model.
Input: data in pickle.
Output: model
Parameters:
Steps:
    1. load pre-trained model
    2. load data
    3. train

"""
import argparse
from deeprho.popgen.utils import check_file_existence
from deeprho.net import recomb_net_1, recomb_net_2
import tensorflow as tf
import logging
import os
import pickle

logger = logging.getLogger(__name__)


def load_data(filename):
    check_file_existence(filename)
    with open(filename, 'rb') as file:
        data = pickle.load(file)
    return data

#
# def load_weights(model, weights):


def _train_large():
    pass


def _train_fine():
    pass



def train(args):
    assert args.out is not None, f'output directory should be specified.'
    assert args.data is not None, f'there is no dataset.'
    assert args.m1 is not None and args.m2 is not None, f'pretrain models should be specified.'
    # check availability of GPU and set to allow memory growth.
    gpus = tf.config.list_physical_devices("GPU")
    if gpus:
        try:
            tf.config.set_visible_devices(gpus[args.gpu], 'GPU')
            tf.config.experimental.set_memory_growth(gpus[args.gpu], True)
        except RuntimeError as e:
            logger.error(e)
            exit(1)
    else:
        logger.warning('no gpu found, use cpu instead.')

    x_train, x_test, y_train, y_test = load_data(args.data)
    print(f'dataset view: train({x_train.shape}), test({x_test.shape})')
    y_train = y_train / args.scale_factor
    y_test = y_test / args.scale_factor
    input_shape = x_train.shape[1:]
    inputs, outputs = recomb_net_1(4, [64,256,256,512], input_shape, 1)
    model = tf.keras.models.Model(inputs=inputs, outputs=outputs)
    optimizer = tf.keras.optimizers.Adam(0.001)
    model.compile(loss='mse', optimizer=optimizer, metrics=['mae'])
    training_epochs = 200
    if not os.path.exists(args.out):
        os.mkdir(args.out)
    for e in range(training_epochs):
        print(f'epoch: {e}')
        hist = model.fit(x_train, y_train, verbose=True, batch_size=200, validation_data=(x_test, y_test))
        model.save(f'{args.out}/model_epoch_{e}.h5')


def gt_args(parser):
    parser.add_argument('--m1', type=str, help='pretrain fine model', default=None)
    parser.add_argument('--m2', type=str, help='pretrain large model', default=None)
    parser.add_argument('--data', type=str, help='dataset', default=None)
    parser.add_argument('--scale-factor', type=float, help='scale factor for large rho', default=10)
    parser.add_argument('--gpu', type=int, help='gpu id', default=0)
    parser.add_argument('--batch-size', type=int, help='batch size for training', default=200)
    parser.add_argument('--out', type=str, help='output directory')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='train')
    gt_args(parser)
    args = parser.parse_args(['--m1', '2',
                              '--m2', '50',
                              '--data', '../garbo/dataset_ld_rf_tri.5',
                              '--scale-factor', '10',
                              '--out', '../garbo/models'
                              ])
    train(args)


