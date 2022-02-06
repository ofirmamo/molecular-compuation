import lab.synthesizer
from lab import synthesizer

from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord


def random_sanity():
    elements = ['a0', 'x1']
    result = synthesizer.synthesize_watson_random(elements, bp_length=6)
    assert len(result) == 2
    assert all([element.seq.crick == '' for element in result.values()])

    print('** End random_sanity **')


def join_sanity():
    elements = ['a0', 'x1']
    synthesized = synthesizer.synthesize_watson_random(elements, bp_length=10)
    assert len(synthesized) == 2

    key = 'a0_x1'
    data = (key, synthesized[elements[0]], synthesized[elements[1]])
    result = synthesizer.join_watson(
        [data],
        take_first=5,
        take_second=5)

    assert len(result) == 1
    assert result[key] is not None
    assert result[key].seq.watson[:5] == data[1].seq.watson[-5:]
    assert result[key].seq.watson[5:] == data[2].seq.watson[:5]
    assert result[key].seq.crick == ''

    print('** End join_sanity **')


def complement_tests():
    a1 = 'ATCG'
    a2 = 'AAGG'
    comps = lab.synthesizer.synthesize_crick_complements(
        [('a1', Dseqrecord(Dseq(watson=a1, linear=True))), ('a2', Dseqrecord(Dseq(watson=a2, linear=True)))])

    assert comps['a1'].seq.crick == 'CGAT'
    assert comps['a1'].seq.watson == ''
    assert comps['a2'].seq.crick == 'CCTT'
    assert comps['a2'].seq.watson == ''

    print('** End complement_tests **')


def synthesize_sat_graph():
    nodes = ['a0', 'x1', '~x1', 'a2', 'x2', '~x2', 'a3']
    synthesized_nodes = lab.synthesizer.synthesize_watson_random(nodes, bp_length=20)
    for node in nodes:
        assert node in synthesized_nodes.keys()
        assert len(synthesized_nodes[node]) == 20

    edges = list(map(lambda edge: ('_'.join(edge), synthesized_nodes[edge[0]], synthesized_nodes[edge[1]]),
                     [('x1', 'a2'), ('a2', 'x2')]))
    synthesized_edges = lab.synthesizer.join_watson(edges, take_first=10, take_second=10)
    edges = list(map(lambda edge: ('_'.join(edge), synthesized_nodes[edge[0]], synthesized_nodes[edge[1]]),
                     [('a0', 'x1')]))
    synthesized_edges.update(lab.synthesizer.join_watson(edges, take_first=20, take_second=10))
    edges = list(map(lambda edge: ('_'.join(edge), synthesized_nodes[edge[0]], synthesized_nodes[edge[1]]),
                     [('~x2', 'a3')]))
    synthesized_edges.update(lab.synthesizer.join_watson(edges, take_first=10, take_second=20))

    assert len(synthesized_edges['x1_a2']) == 20
    assert synthesized_edges['x1_a2'].seq.watson[:10] == synthesized_nodes['x1'].seq.watson[10:] and \
           synthesized_edges['x1_a2'].seq.watson[10:] == synthesized_nodes['a2'].seq.watson[:10]
    assert synthesized_edges['x1_a2'].seq.crick == ''

    assert len(synthesized_edges['a2_x2']) == 20
    assert synthesized_edges['a2_x2'].seq.watson[:10] == synthesized_nodes['a2'].seq.watson[10:] and \
           synthesized_edges['a2_x2'].seq.watson[10:] == synthesized_nodes['x2'].seq.watson[:10]
    assert synthesized_edges['a2_x2'].seq.crick == ''

    assert len(synthesized_edges['a0_x1']) == 30
    assert synthesized_edges['a0_x1'].seq.watson[:20] == synthesized_nodes['a0'].seq.watson and \
           synthesized_edges['a0_x1'].seq.watson[20:] == synthesized_nodes['x1'].seq.watson[:10]
    assert synthesized_edges['a0_x1'].seq.crick == ''

    assert len(synthesized_edges['~x2_a3']) == 30
    assert synthesized_edges['~x2_a3'].seq.watson[:10] == synthesized_nodes['~x2'].seq.watson[10:] and \
           synthesized_edges['~x2_a3'].seq.watson[10:] == synthesized_nodes['a3'].seq.watson
    assert synthesized_edges['~x2_a3'].seq.crick == ''

    print('** End synthesize_sat_graph **')


random_sanity()
join_sanity()
synthesize_sat_graph()
complement_tests()