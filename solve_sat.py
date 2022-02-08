import networkx
import networkx as nx
from matplotlib import pyplot as plt

from lab.gel_electrophoresis import gel_electrophoresis
from lab import synthesizer, ligation
from lab.magnetic_beads_filter import magnetic_beads_filter
from lab.pcr import pcr, create_primers
from sat.instance import SATInstance

SAT_BASE = 'DNAComputationSim/sat'
SAT_BP_LENGTH = 20


def synthesize_nodes(graph):
    nodes = synthesizer.synthesize_watson_random(graph.nodes(), bp_length=SAT_BP_LENGTH)
    print(f'Nodes synthesized')
    return nodes


def synthesize_edges(graph, synthesized_nodes):
    edges = list(map(lambda edge: ('_'.join(edge), synthesized_nodes[edge[0]], synthesized_nodes[edge[1]]),
                     list(graph.edges)[2:-2]))
    synthesized_edges = synthesizer.join_watson(edges, take_first=SAT_BP_LENGTH // 2, take_second=SAT_BP_LENGTH // 2)
    edges = list(map(lambda edge: ('_'.join(edge), synthesized_nodes[edge[0]], synthesized_nodes[edge[1]]),
                     list(graph.edges)[:2]))
    synthesized_edges.update(synthesizer.join_watson(edges, take_first=SAT_BP_LENGTH, take_second=SAT_BP_LENGTH // 2))
    edges = list(map(lambda edge: ('_'.join(edge), synthesized_nodes[edge[0]], synthesized_nodes[edge[1]]),
                     list(graph.edges)[-2:]))
    synthesized_edges.update(synthesizer.join_watson(edges, take_first=SAT_BP_LENGTH // 2, take_second=SAT_BP_LENGTH))
    print(f'Edges synthesized')
    return synthesized_edges


def create_graph(args, out_file):
    graph = networkx.read_edgelist(out_file, create_using=networkx.DiGraph(), nodetype=str, data=[('to', str)])
    nx.draw_networkx(graph, with_labels=True)
    plt.savefig(f'{SAT_BASE}/{args.filename}.png')
    plt.clf()
    print('3SAT graph image saved..')
    return graph


def init(args, out_file):
    instance = SATInstance.from_file(f'{SAT_BASE}/{args.filename}')
    instance.to_assignment_graph_file(instance.variables, out_file)
    print('SAT input file parsed..')
    return instance


def run(args):
    out_file = f'{SAT_BASE}/{args.filename}.out'

    instance = init(args, out_file)
    graph = create_graph(args=args, out_file=out_file)
    nodes = synthesize_nodes(graph=graph)
    edges = synthesize_edges(graph=graph, synthesized_nodes=nodes)
    complement_nodes = synthesizer.synthesize_crick_complements(list(nodes.items()))
    reversed_nodes_dict = {value.seq.watson: key for (key, value) in nodes.items()}

    a_node_names = sorted([node_name for node_name in nodes.keys() if node_name.startswith('a')])
    a0 = nodes[a_node_names[0]]
    an_complement = complement_nodes[a_node_names[-1]]

    # Create sufficient DNA to generate candidates.
    tube = (list(edges.values()) + list(complement_nodes.values()))
    ligation.ligate(tube, limit=10, rounds=150000)
    for clause in instance.clauses:
        Ti = []
        for literal, state in clause.items():
            Ti.extend(magnetic_beads_filter(tube, nodes[literal if state else f"~{literal}"].seq.watson))
        for mol in Ti:
            if len(mol.seq.watson) == 140 and len(mol.seq.crick):
                print(f'First iteration path: {sequence_result(mol, reversed_nodes_dict)}')
        tube = pcr(Ti, create_primers(a0, an_complement), rounds=1)
    res = gel_electrophoresis(tube)

    if any(res):
        print('SAT is satisfiable')
        path = gel_electrophoresis(res, variant='max')
        dna_desequenced = sequence_result(path, reversed_nodes_dict)
        print(f'assignement: {dna_desequenced}')
    else:
        print('SAT is not satisfiable')


def sequence_result(path, reversed_nodes_dict):
    dna_desequenced = []
    for i in range(0, len(path), 20):
        node = reversed_nodes_dict.get(path[i:i + 20].seq.watson)
        if not node.startswith('a'):
            dna_desequenced.append(node)
    return dna_desequenced