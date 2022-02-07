from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord

import lab.cleave


def site():
    def procedure(molecule):
        return Dseqrecord(Dseq(
            watson=molecule.seq.watson[:5],
            crick=molecule.seq.crick[-5:],
            ovhg=0, linear=True))

    base = lab.cleave.FOKI_BASE_SEQ
    assert base == procedure(Dseqrecord(Dseq(watson='GGATGAATT', crick='AATTCATCC', ovhg=0, linear=True)))
    assert base == procedure(base)
    print('** End site **')


def FokI():
    lab.cleave.FAIL_PROB = 0
    m = Dseqrecord(Dseq(watson='GGATG' + (14 * 'T'), crick=(14 * 'T' + 'CATCC'), ovhg=0, linear=True))

    result = lab.cleave.FokI([m])
    assert len(result) == 2
    print('** End FokI **')


site()
FokI()
