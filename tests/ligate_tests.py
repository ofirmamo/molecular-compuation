import lab.ligation

from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord

A = "A"
T = "T"
C = "C"
G = "G"
base_list = [A, T, C, G]
COMPLEMENTS = {
    A: T,
    G: C,
    T: A,
    C: G
}


def should_ligate():
    assert lab.ligation.should_ligate(Dseqrecord(Dseq(watson='AAAAA', crick='TTTTT'))) is False
    assert lab.ligation.should_ligate(Dseqrecord(Dseq(watson='GAAAA', crick='TTTT'))) is True
    assert lab.ligation.should_ligate(Dseqrecord(Dseq(watson='AAAAG', crick='TTTT'))) is True

    print('** End should_ligate **')


def end_watson_match():
    a1 = Dseqrecord(Dseq(watson='ACCC', crick='GTT', ovhg=1))
    assert lab.ligation.end_watson_match(a1, 'C')
    assert lab.ligation.end_watson_match(a1, 'CC')
    assert not lab.ligation.end_watson_match(a1, 'CCC')

    a2 = Dseqrecord(Dseq(watson='ACCC', crick='GGGT', ovhg=0))
    assert not lab.ligation.end_watson_match(a2, 'C')
    assert not lab.ligation.end_watson_match(a2, 'CC')
    assert not lab.ligation.end_watson_match(a2, 'CCC')

    a3 = Dseqrecord(Dseq(watson='ACCC', crick='GTT', ovhg=0))
    assert lab.ligation.end_watson_match(a3, 'C')
    assert not lab.ligation.end_watson_match(a3, 'CC')
    assert not lab.ligation.end_watson_match(a3, 'CCC')

    a4 = Dseqrecord(Dseq(watson='ACCC', crick='GGGT', ovhg=-1))
    assert not lab.ligation.end_watson_match(a4, 'C')
    assert not lab.ligation.end_watson_match(a4, 'CC')
    assert not lab.ligation.end_watson_match(a4, 'CCC')

    a5 = Dseqrecord(Dseq(watson='AGG', crick='C', ovhg=-2))
    assert not lab.ligation.end_watson_match(a5, 'AG')

    print('** End end_watson_match **')


def start_watson_match():
    a1 = Dseqrecord(Dseq(watson='ACCC', crick='TGG', ovhg=-2))
    assert lab.ligation.start_watson_match(a1, 'A')
    assert lab.ligation.start_watson_match(a1, 'AC')
    assert not lab.ligation.start_watson_match(a1, 'C')
    assert not lab.ligation.start_watson_match(a1, 'CCC')

    print('** End start_watson_match **')


def end_crick_match():
    # lookup = TTC
    a1 = Dseqrecord(Dseq(watson='AAG', crick='CTTCC', ovhg=-1))
    assert lab.ligation.end_crick_match(a1, 'TTC')
    assert lab.ligation.end_crick_match(a1, 'TC')
    assert lab.ligation.end_crick_match(a1, 'C')
    assert not lab.ligation.end_crick_match(a1, 'TT')
    assert not lab.ligation.end_crick_match(a1, 'TTCT')

    print('** End end_crick_match **')


def start_crick_match():
    a1 = Dseqrecord(Dseq(watson='GG', crick='CCTTCT', ovhg=4))
    assert lab.ligation.start_crick_match(a1, 'TCTT')
    assert lab.ligation.start_crick_match(a1, 'TCT')
    assert lab.ligation.start_crick_match(a1, 'TC')
    assert lab.ligation.start_crick_match(a1, 'T')
    assert not lab.ligation.start_crick_match(a1, 'TT')
    assert not lab.ligation.start_crick_match(a1, 'CTT')

    print('** End start_crick_match **')


def ligation_single_watson_end():
    fragment = Dseqrecord(Dseq(watson='AAGG', crick='CCTTCT', ovhg=2))
    candidates = [Dseqrecord(Dseq(watson='GAGG', crick='C', ovhg=3)), Dseqrecord(Dseq(watson='GAG', crick='C', ovhg=0))]
    assert lab.ligation.ligation_single_watson_end(fragment, candidates, 1)
    assert candidates[0] == Dseqrecord(Dseq(watson='GAGG', crick='C', ovhg=3))
    assert candidates[1] == Dseqrecord(Dseq(watson='GAGAAGG', crick='CCTTCTC', ovhg=0))

    candidates = [Dseqrecord(Dseq(watson='GAGG', crick='C', ovhg=3)), Dseqrecord(Dseq(watson='GA', crick='C', ovhg=0))]
    assert lab.ligation.ligation_single_watson_end(fragment, candidates, 1)
    assert candidates[0] == Dseqrecord(Dseq(watson='GAGG', crick='C', ovhg=3))
    assert candidates[1] == Dseqrecord(Dseq(watson='GA AAGG', crick='CCTTCTC', ovhg=0))

    candidates = [Dseqrecord(Dseq(watson='GAGG', crick='C', ovhg=3))]
    assert not lab.ligation.ligation_single_watson_end(fragment, candidates, 1)
    assert candidates[0] == Dseqrecord(Dseq(watson='GAGG', crick='C', ovhg=3))


def ligation_single_crick_end():
    fragment = Dseqrecord(Dseq(watson='AAGGAGG', crick='CCTC', ovhg=-3))
    candidates = [Dseqrecord(Dseq(watson='AG', crick='', ovhg=0)),
                  Dseqrecord(Dseq(watson='AG', crick='CTTCTT', ovhg=1))]
    assert lab.ligation.ligation_single_crick_end(fragment, candidates, limit=1)
    assert candidates[0] == Dseqrecord(Dseq(watson='AG', crick='', ovhg=0))
    assert candidates[1] == Dseqrecord(Dseq(watson='AGAAGGAGG', crick='CCTCCTTCTT', ovhg=1))

    candidates = [Dseqrecord(Dseq(watson='AG', crick='', ovhg=0)),
                  Dseqrecord(Dseq(watson='AG', crick='TTCT', ovhg=0))]
    assert lab.ligation.ligation_single_crick_end(fragment, candidates, limit=1)
    assert candidates[0] == Dseqrecord(Dseq(watson='AG', crick='', ovhg=0))
    assert candidates[1] == Dseqrecord(Dseq(watson='AGAAGGAGG', crick='CCTC TTCT', ovhg=0))

    candidates = [Dseqrecord(Dseq(watson='AG', crick='', ovhg=0)),
                  Dseqrecord(Dseq(watson='AG', crick='TC', ovhg=-1))]
    assert lab.ligation.ligation_single_crick_end(fragment, candidates, limit=1)
    assert candidates[0] == Dseqrecord(Dseq(watson='AG', crick='', ovhg=0))
    assert candidates[1] == Dseqrecord(Dseq(watson='AGAAGGAGG', crick='CCTC  TC', ovhg=-1))

    print('** End ligation_single_crick_end **')


def ligation_single_crick_start():
    fragment = Dseqrecord(Dseq(watson='AAGGAG', crick='CCTT', ovhg=0))
    candidates = [Dseqrecord(Dseq(watson='AA', crick='TTCC', ovhg=2)),
                  Dseqrecord(Dseq(watson='AGG', crick='CCTCT', ovhg=2))]
    assert lab.ligation.ligation_single_crick_start(fragment, candidates, limit=1)
    assert candidates[0] == Dseqrecord(Dseq(watson='AA', crick='TTCC', ovhg=2))
    assert candidates[1] == Dseqrecord(Dseq(watson='AAGGAGAGG', crick='CCTCTCCTT', ovhg=0))

    candidates = [Dseqrecord(Dseq(watson='AGG', crick='CCTC', ovhg=1))]
    assert lab.ligation.ligation_single_crick_start(fragment, candidates, limit=1)
    assert candidates[0] == Dseqrecord(Dseq(watson='AAGGAGAGG', crick='CCTC CCTT', ovhg=0))

    print('** End ligation_single_crick_start **')


def ligation_single_watson_start():
    fragment = Dseqrecord(Dseq(watson='AAGG', crick='CTCTCCTTT', ovhg=1))
    candidates = [Dseqrecord(Dseq(watson='AGAA', crick='TTCT', ovhg=0)),
                  Dseqrecord(Dseq(watson='AGAGAA', crick='TT', ovhg=-4))]
    assert lab.ligation.ligation_single_watson_start(fragment, candidates, limit=1)
    assert candidates[0] == Dseqrecord(Dseq(watson='AGAA', crick='TTCT', ovhg=0))
    assert candidates[1] == Dseqrecord(Dseq(watson='AAGGAGAGAA', crick='TTCTCTCCTTT', ovhg=1))

    candidates = [Dseqrecord(Dseq(watson='AGAA', crick='TTCT', ovhg=0)),
                  Dseqrecord(Dseq(watson='GAGAA', crick='TT', ovhg=-3))]
    assert lab.ligation.ligation_single_watson_start(fragment, candidates, limit=1)
    assert candidates[0] == Dseqrecord(Dseq(watson='AGAA', crick='TTCT', ovhg=0))
    assert candidates[1] == Dseqrecord(Dseq(watson='AAGG GAGAA', crick='TTCTCTCCTTT', ovhg=1))

    print('** End ligation_single_watson_start **')


def sanity_sat_edges():
    a0_x1 = Dseqrecord(Dseq(watson='AAAAGG', crick='', ovhg=0))
    a0_comp = Dseqrecord(Dseq(watson='', crick='TTTT', ovhg=0))
    x1_comp = Dseqrecord(Dseq(watson='', crick='CCCC', ovhg=0))

    candidates = [x1_comp, a0_comp]
    assert lab.ligation.ligation_single_fragment(a0_x1, candidates, limit=2)
    assert candidates[0] == a0_comp
    assert candidates[1] == Dseqrecord(Dseq(watson='AAAAGG', crick='CCCC', ovhg=-4, linear=True))

    candidates_copy = candidates.copy()
    fragment = candidates_copy[0]
    candidates_new = candidates_copy[1:]
    assert lab.ligation.ligation_single_fragment(fragment, candidates_new, limit=2)
    assert candidates_new[0] == Dseqrecord(Dseq(watson='AAAAGG', crick='CCCCTTTT', ovhg=0))

    candidates_copy = candidates.copy()
    fragment = candidates_copy[1]
    candidates_new = candidates_copy[:1]
    assert lab.ligation.ligation_single_fragment(fragment, candidates_new, limit=2)
    assert candidates_new[0] == Dseqrecord(Dseq(watson='AAAAGG', crick='CCCCTTTT', ovhg=0))


def complex_sat_edges():
    A0 = "GGGG"
    A1 = "GAGA"
    X1 = "AAAA"
    getCrickComp = lambda x: "".join(reversed([COMPLEMENTS[l] for l in x]))
    A0COMP = getCrickComp(A0)
    A1COMP = getCrickComp(A1)
    X1COMP = getCrickComp(X1)
    start = lambda x: x[:2]
    end = lambda x: x[2:]
    edge = lambda u, v: end(u) + start(v)
    edge_start = lambda u, v: u + start(v)

    a0_x1 = Dseqrecord(Dseq(watson=A0 + start(X1), crick='', ovhg=0))
    x1_a1 = Dseqrecord(Dseq(watson=start(X1) + end(A1), crick='', ovhg=0))
    a0_comp = Dseqrecord(Dseq(watson='', crick=A0COMP, ovhg=0))
    a1_comp = Dseqrecord(Dseq(watson='', crick=A1COMP, ovhg=0))
    x1_comp = Dseqrecord(Dseq(watson='', crick=X1COMP, ovhg=0))

    candidates = [x1_a1, a0_x1, a0_comp, a1_comp]
    assert lab.ligation.ligation_single_fragment(x1_comp, candidates, limit=2)
    assert candidates[-1] == Dseqrecord(Dseq(watson=edge(X1, A1), crick=X1COMP, ovhg=2))

    fragment = candidates[0]
    candidates = candidates[1:]
    assert lab.ligation.ligation_single_fragment(fragment, candidates, limit=2)
    assert candidates[-1] == Dseqrecord(Dseq(watson=edge_start(A0, X1) + edge(X1, A1), crick=X1COMP, ovhg=-4))

    fragment = candidates[0]
    candidates = candidates[1:]
    assert lab.ligation.ligation_single_fragment(fragment, candidates, limit=2)
    assert candidates[-1] == Dseqrecord(
        Dseq(watson=edge_start(A0, X1) + edge(X1, A1), crick=(A0COMP + X1COMP)[::-1], ovhg=0))

    fragment = candidates[0]
    candidates = candidates[1:]
    assert lab.ligation.ligation_single_fragment(fragment, candidates, limit=2)
    assert candidates[-1] == Dseqrecord(
        Dseq(watson=edge_start(A0, X1) + edge(X1, A1), crick="".join(reversed([A0COMP, X1COMP, A1COMP])), ovhg=0))

    print('** complex_sat_edges **')


def full_batch_test():
    A0 = "GGGG"
    A1 = "GAGA"
    X1 = "AAAA"
    getCrickComp = lambda x: "".join(reversed([COMPLEMENTS[l] for l in x]))
    A0COMP = getCrickComp(A0)
    A1COMP = getCrickComp(A1)
    X1COMP = getCrickComp(X1)
    start = lambda x: x[:2]
    end = lambda x: x[2:]
    edge = lambda u, v: end(u) + start(v)
    edge_start = lambda u, v: u + start(v)

    a0_x1 = Dseqrecord(Dseq(watson=A0 + start(X1), crick='', ovhg=0))
    x1_a1 = Dseqrecord(Dseq(watson=start(X1) + end(A1), crick='', ovhg=0))
    a0_comp = Dseqrecord(Dseq(watson='', crick=A0COMP, ovhg=0))
    a1_comp = Dseqrecord(Dseq(watson='', crick=A1COMP, ovhg=0))
    x1_comp = Dseqrecord(Dseq(watson='', crick=X1COMP, ovhg=0))

    tube = [a0_x1, x1_a1, a0_comp, a1_comp, x1_comp]
    lab.ligation.ligate(tube, limit=2)
    assert tube[-1] == Dseqrecord(
        Dseq(watson=edge_start(A0, X1) + edge(X1, A1), crick="".join(reversed([A0COMP, X1COMP, A1COMP])), ovhg=0))
    print('** End full_batch_test **')


should_ligate()
end_watson_match()
start_watson_match()
end_crick_match()
start_crick_match()
ligation_single_watson_end()
ligation_single_crick_end()
ligation_single_crick_start()
ligation_single_watson_start()
sanity_sat_edges()
complex_sat_edges()
full_batch_test()
