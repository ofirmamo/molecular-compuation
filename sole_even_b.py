from typing import List

from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord

import lab.common, lab.gel_electrophoresis, lab.ligation
from lab.cleave import FOKI_BASE_SEQ

START_MSG = 'DNA sequence is:'

BASE = 'DNAComputationSim/simulator/b-even'
BP_LENGTH = 20
TERMINAL = Dseqrecord(Dseq(watson='TGTCGC', crick='GCGACA', ovhg=0, linear=True))
letter_mapping = {
    'a': Dseqrecord(Dseq(watson='CTGGCT', crick='AGCCAG', ovhg=0, linear=True)),
    'b': Dseqrecord(Dseq(watson='CGCAGC', crick='GCTGCG', ovhg=0, linear=True))
}
SOFTWARE = [
    Dseqrecord(Dseq(watson='GGATGTAC', crick='AGCCGTACATCC', ovhg=0, linear=True)),  # q0 => q0
    Dseqrecord(Dseq(watson='GGATGACGAC', crick='GCTGGTCGTCATCC', ovhg=0, linear=True)),  # q0 => q1
    Dseqrecord(Dseq(watson='GGATGG', crick='TGCGCCATCC', ovhg=0, linear=True)),  # q1 => q0
    Dseqrecord(Dseq(watson='GGATGACG', crick='CCAGCGTCATCC', ovhg=0, linear=True)),  # q1 => q1
]


def merge(s1, s2):
    return Dseqrecord(Dseq(
        watson=s1.seq.watson + s2.seq.watson,
        crick=s2.seq.crick + s1.seq.crick,
        linear=True,
        ovhg=0))


def create_start_seq(buf):
    return Dseqrecord(Dseq(
        watson=FOKI_BASE_SEQ.seq.watson + lab.common.COMPLEMENTS[lab.common.A] * buf,
        crick=(lab.common.A * buf) + FOKI_BASE_SEQ.seq.crick,
        ovhg=0,
        linear=True))


def create_message(*args, upper=True, lower=True):
    max_len = max([len(arg) for arg in args])

    if upper:
        print(max_len * "*")

    for arg in args:
        print(arg)

    if lower:
        print(max_len * "*")


def crick(seq):
    return seq.seq.crick[::-1]


def run(args):
    with open(f'{BASE}/{args.filename}') as f:
        letters = f.read().strip()
    create_message(f'Input is: {letters}')

    seq = create_start_seq(buf=1 if letters[0] == 'a' else 3)
    for letter in letters:
        seq = merge(seq, letter_mapping[letter])
    seq = merge(seq, TERMINAL)
    create_message(f'DNA sequence is:', seq.seq.watson, crick(seq))

    tube = [Dseqrecord(seq)]
    while tube_length := len(lab.gel_electrophoresis.gel_electrophoresis(tube, variant='min')) > 6:
        tube = lab.cleave.FokI(tube)
        tube_length = len(lab.gel_electrophoresis.gel_electrophoresis(tube, variant='min'))
        if tube_length <= 6:
            break

        tube = tube + (SOFTWARE * 5)
        lab.ligation.ligate(tube, limit=4, rounds=100)
        tube = [lab.gel_electrophoresis.gel_electrophoresis(tube, variant='max')]

    create_message('Final state:',
                   'q0 (There is even occurrences of "b")'
                   if tube_length == 4
                   else 'q1 (There is odd occurrences of "b")')
