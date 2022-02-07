from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord
from lab.magnetic_beads_filter import magnetic_beads_filter


def sanity():
    e1 = Dseqrecord(Dseq(watson='AAGG', crick='TT', ovhg=0, linear=True))
    e2 = Dseqrecord(Dseq(watson='AAGG', crick='CT', ovhg=-1, linear=True))
    e3 = Dseqrecord(Dseq(watson='AAGG', crick='T', ovhg=-1, linear=True))
    e4 = Dseqrecord(Dseq(watson='AAGG', crick='TTT', ovhg=1, linear=True))

    tube = [e1, e2, e3, e4]
    assert len(magnetic_beads_filter(tube, 'AAG')) == 4
    assert len(magnetic_beads_filter(tube, 'T')) == 4
    assert len(magnetic_beads_filter(tube, 'CT')) == 1

    print('** End sanity **')


sanity()
