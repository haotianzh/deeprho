import logging
import argparse

from deeprho.estimate import gt_args as est_args
from deeprho.estimate import run as est_run
from deeprho.data_provider_parallel import gt_args as dpp_args
from deeprho.data_provider_parallel import run as  dpp_run


subcommands_entries = { 'estimate': est_run, 'simulate': dpp_run}
subcommands_args = {'estimate': est_args, 'simulate': dpp_args}
subcommands_helps = {'estimate': 'recombination rate estimation',
                     'simulate': 'simulation for re-training'}

def main():
    parser = argparse.ArgumentParser(
        description='deeprho is used for estimating recombination rate given population genetic data')
    subparser = parser.add_subparsers(dest='command')
    # adding arguments for subcommands
    for subcommand in subcommands_args:
        sbc = subparser.add_parser(subcommand, help=subcommands_helps[subcommand])
        subcommands_args[subcommand](sbc)
    args = parser.parse_args()
    entry = subcommands_entries[args.command]
    entry(args)

if __name__ == '__main__':
    main()