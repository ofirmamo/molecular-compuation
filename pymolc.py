import argparse
import sys

import sole_even_b
import solve_sat

parser = argparse.ArgumentParser()
subparsers = parser.add_subparsers()

# 3SAT
sat_argparse = subparsers.add_parser('sat')
sat_argparse.add_argument('--filename', help='name of SAT file')
sat_argparse.set_defaults(func=solve_sat.run)

sat_argparse = subparsers.add_parser('even-b')
sat_argparse.add_argument('--filename', help='name of input file')
sat_argparse.set_defaults(func=sole_even_b.run)

args = parser.parse_args(sys.argv[1:])
args.func(args)
