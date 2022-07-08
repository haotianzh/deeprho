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
import tensorflow as tf
import logging
import os
import pickle

def load_data(filename):
    check_file_existence(filename)

tf.keras.applications.ResNet50()




def train(args):
    assert args.out is not None, f'output directory should be specified.'
    assert args.data is not None, f'there is no dataset.'
    assert args.m1 is not None and args.m2 is not None, f'pretrain models should be specified.'

    x_train, x_test, y_train, y_test = load_data(args.data)


def gt_args(parser):
    parser.add_argument('--m1', type=str, help='pretrain fine model', default=None)
    parser.add_argument('--m2', type=str, help='pretrain large model', default=None)
    parser.add_argument('--data', type=str, help='dataset', default=None)
    parser.add_argument('--scale-y', type=float, help='scale factor for large rho', default=10)
    parser.add_argument('--out', type=str, help='output directory')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='train')
    gt_args(parser)
    args = parser.parse_args()

