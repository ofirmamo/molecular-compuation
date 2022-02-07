from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord

from lab import pcr


def amplify_watson():
    watson = 'AAGGAG'
    crick = 'CTCCTT'

    success, result = pcr.anneal_and_amplify_watson(watson, 'T')
    assert success
    assert result == Dseqrecord(Dseq(watson=watson, crick=crick, ovhg=0, linear=True))

    success, result = pcr.anneal_and_amplify_watson(watson, 'TT')
    assert success
    assert result == Dseqrecord(Dseq(watson=watson, crick=crick, ovhg=0, linear=True))

    success, result = pcr.anneal_and_amplify_watson(watson, 'C')
    assert not success

    success, result = pcr.anneal_and_amplify_watson(watson, 'CC')
    assert not success

    success, result = pcr.anneal_and_amplify_watson(watson, 'TTC')
    assert success
    assert result == Dseqrecord(Dseq(watson=watson, crick=crick, ovhg=0, linear=True))

    print('** End amplify_watson **')


def amplify_crick():
    watson = 'AAGGAG'
    crick = 'CTCCTT'

    success, result = pcr.anneal_and_amplify_crick(crick, 'A')
    assert success
    assert result == Dseqrecord(Dseq(watson=watson, crick=crick, ovhg=0, linear=True))

    success, result = pcr.anneal_and_amplify_crick(crick, 'AA')
    assert success
    assert result == Dseqrecord(Dseq(watson=watson, crick=crick, ovhg=0, linear=True))

    success, result = pcr.anneal_and_amplify_crick(crick, 'G')
    assert not success

    success, result = pcr.anneal_and_amplify_crick(crick, 'GG')
    assert not success

    success, result = pcr.anneal_and_amplify_crick(crick, 'AAG')
    assert success
    assert result == Dseqrecord(Dseq(watson=watson, crick=crick, ovhg=0, linear=True))

    print('** End amplify_crick **')


def heat_tube():
    pcr.fail_prob = 0
    e1 = Dseqrecord(Dseq(watson='AAGG', crick='CCTT', ovhg=0, linear=0))
    heated = pcr.heat_tube([e1])

    assert len(heated) == 2
    assert Dseqrecord(Dseq(watson='AAGG', crick='', ovhg=0, linear=0)) in heated
    assert Dseqrecord(Dseq(watson='', crick='CCTT', ovhg=0, linear=0)) in heated

    print('** End heat_tube **')


def pcr_round():
    pcr.fail_prob = 0
    e1 = Dseqrecord(Dseq(watson='AAGG', crick='', ovhg=0, linear=0))
    e2 = Dseqrecord(Dseq(watson='', crick='CCTT', ovhg=0, linear=0))
    pcred = pcr.pcr_round(tube=[e1, e2], primers=['TT', 'AA'])

    assert len(pcred) == 2
    assert all(x == Dseqrecord(Dseq(watson='AAGG', crick='CCTT', ovhg=0, linear=0)) for x in pcred)

    pcred = pcr.pcr_round(tube=[e1, e2], primers=['TT', 'CC'])
    assert len(pcred) == 2
    assert Dseqrecord(Dseq(watson='AAGG', crick='CCTT', ovhg=0, linear=0)) in pcred
    assert Dseqrecord(Dseq(watson='', crick='CCTT', ovhg=0, linear=0)) in pcred

    pcred = pcr.pcr_round(tube=[e1, e2], primers=['GG', 'CC'])
    assert len(pcred) == 2
    assert Dseqrecord(Dseq(watson='AAGG', crick='', ovhg=0, linear=0)) in pcred
    assert Dseqrecord(Dseq(watson='', crick='CCTT', ovhg=0, linear=0)) in pcred

    print('** End pcr_round **')


amplify_watson()
amplify_crick()
heat_tube()
pcr_round()
