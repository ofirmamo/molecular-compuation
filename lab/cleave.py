from typing import List

from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord

FOKI_BASE_SEQ = Dseqrecord(Dseq(watson='GGATG', crick='CATCC', ovhg=0, linear=True))
WATSON_CLEAVE_BASES = 9
CRICK_CLEAVE_BASES = 13


def FokI(tube: List[Dseqrecord]):
    result = []
    for molecule in tube:
        if FOKI_BASE_SEQ == Dseqrecord(Dseq(
                watson=molecule.seq.watson[:5],
                crick=molecule.seq.crick[-5:],
                ovhg=0,
                linear=True)):

            e = Dseqrecord(Dseq(
                watson=molecule.seq.watson[len(FOKI_BASE_SEQ) + WATSON_CLEAVE_BASES:],
                crick=molecule.seq.crick[:-len(FOKI_BASE_SEQ) - CRICK_CLEAVE_BASES],
                ovhg=-4,
                linear=True))
            result.append(e)
        else:
            result.append(molecule)
    return result
