"""
Author: Haotian Z
Entrance of deeprho.
All ops of deeprho can be executed by typing: deeprho [subcommand] [arguments]

There are 4 subcommands so far: estimate, simulate, train, test.

"""
import logging

import argparse
from deeprho import __version__
from deeprho.estimate import gt_args as est_args
from deeprho.estimate import run as est_run
from deeprho.data_provider_parallel import gt_args as dpp_args
from deeprho.data_provider_parallel import run as dpp_run


subcommands_entries = { 'estimate': est_run, 'simulate': dpp_run}
subcommands_args = {'estimate': est_args, 'simulate': dpp_args}
subcommands_helps = {'estimate': 'recombination rate estimation',
                     'simulate': 'simulation for re-training'}


def main():
    parser = argparse.ArgumentParser(
        description='deeprho is used for estimating recombination rate given population genetic data')
    parser.add_argument('-v', '--version', action='store_true', help='show version')
    subparser = parser.add_subparsers(dest='command')
    # adding arguments for subcommands
    for subcommand in subcommands_args:
        sbc = subparser.add_parser(subcommand, help=subcommands_helps[subcommand])
        subcommands_args[subcommand](sbc)

    args = parser.parse_args()
    if args.version:
        print(f'deeprho v{__version__}')
        exit()
    if args.command is not None:
        subcommands_entries[args.command](args)
    else:
        parser.print_help()
        exit()


if __name__ == '__main__':
    main()