import argparse
import sys

import solve_sat

parser = argparse.ArgumentParser()
subparsers = parser.add_subparsers()

# 3SAT
sat_argparse = subparsers.add_parser('sat')
sat_argparse.add_argument('--filename', help='name of SAT file')
sat_argparse.set_defaults(func=solve_sat.run)

args = parser.parse_args(sys.argv[1:])
args.func(args)
