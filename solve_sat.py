import networkx
import networkx as nx
from matplotlib import pyplot as plt

from lab import synthesizer, ligation
from sat.instance import SATInstance

SAT_BASE = 'DNAComputationSim/simulator/sat'
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

    # Create sufficient DNA to generate candidates.
    tube = (list(edges.values()) + list(complement_nodes.values()))
    ligation.ligate(tube, limit=10)
    pass
