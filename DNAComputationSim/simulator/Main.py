# This is the main method to run the DNA encoding algorithm
#
#
# author: Jack Burns
# create date: 11/13/2017
# version 1.0
import sys

import matplotlib.pyplot as plt
import networkx as nx
from pydna.amplify import Anneal
from pydna.assembly import Assembly
from pydna.dseqrecord import Dseqrecord as DnaSeqRecord

from Graph import Graph
from MoleculeEncoder import Encoder

enc = Encoder()


def main():
    # arguments: graph_name, start_v, end_v
    graphs = Graph()
    if len(sys.argv) < 4:
        print("Input is not valid")
        return exit(1)
    else:
        try:
            graph = getattr(graphs, sys.argv[1])
        except FileNotFoundError:
            print("Wrong file or file path")
            return exit(1)

        else:
            start_vertex = str(sys.argv[2]).upper()
            end_vertex = str(sys.argv[3]).upper()
            encoded_nodes = enc.encodeNodes(graph)

            if start_vertex not in encoded_nodes:
                print("could not identify start vertex")
                return exit(1)
            elif end_vertex not in encoded_nodes:
                print("could not identify end vertex")
                return exit(1)
            else:
                ham_graph_path = 'graphs/ham_path'
                encoded_edges = enc.encodeEdges(graph, encoded_nodes)
                node_names = [node for node in encoded_nodes if node != start_vertex and node != end_vertex]
                nx.draw_circular(graph, with_labels=True)
                plt.savefig('graphs/graph.png')
                plt.clf()

                x = createTubeSolution(encoded_nodes, encoded_edges, start_vertex, end_vertex)
                filtered_paths = filter(x, encoded_nodes, node_names)
                results = getResults(encoded_nodes, encoded_edges, filtered_paths, ham_graph_path)

                if results[0] == True:
                    print("\nTHERE IS A HAMILTONIAN PATH " + str(results[1]))
                    return exit(0)

                else:
                    print("\nTHERE IS NOT HAMILTONIAN PATH")
                    return exit(0)


def getResults(nodes, edges, filtered_results, file_name):
    ham_path = False
    ham_path_tup = (0, 0)
    if len(filtered_results) != 0:
        for result in filtered_results:
            print("\nHAMILTONIAN PATH: " + str(result.seq))
            ham_path_tup = extractEdges(str(result.seq), nodes, edges)
            nx.draw_circular(ham_path_tup[0], with_labels=True)
            plt.savefig(file_name)
            plt.clf()
        ham_path = True
    return (ham_path, ham_path_tup[1])


def createTubeSolution(encoded_nodes, encoded_edges, start, end):
    dna_sequences = enc.toDnaSequence(encoded_nodes, encoded_edges)
    p1 = DnaSeqRecord(encoded_nodes[start])
    p2 = DnaSeqRecord(enc.getSeqComplement(encoded_nodes[end]))
    tube_solution = Assembly(dna_sequences, limit=10)
    print("\n" + str(tube_solution) + "\n")
    candidates = []

    for i in range(len(tube_solution.linear_products)):
        product = tube_solution.linear_products[i]
        template = DnaSeqRecord(product)
        pcr = Anneal([p1, p2], template, limit=10)
        gel = len(encoded_nodes) * enc.SEQ_LEN

        if len(pcr.products) != 0:
            print(product.detailed_figure())
            print(product.figure())
            for p in pcr.products:
                if len(p.seq) == gel:
                    p.seq = p.seq[10:]
                    p.seq = p.seq[:-10]
                    candidates.append(p)

    # print("\n" +str(nodes))
    # print(str(edges) +"\n")
    return candidates


def filter(paths, nodes, node_name_list):
    n = len(node_name_list)
    unique_paths = list(set(paths))
    gp = []
    pp = {}

    for path in unique_paths:
        pp[path] = 0

    while len(node_name_list) > 0:
        x = node_name_list.pop()
        for path in unique_paths:
            if nodes[x] in path:
                pp[path] += 1

    for path in pp:
        if pp[path] == n:
            gp.append(path)

    return gp


def extractEdges(path, node_list, edge_list):
    g = nx.DiGraph()
    order_edges = [0] * int(len(path) / enc.SEQ_LEN)
    for node in node_list:
        g.add_node(node)

    for edge in edge_list:
        if edge_list[edge] in path:
            index = int(path.index(edge_list[edge]) / enc.SEQ_LEN)
            order_edges[index] = edge
            g.add_edge(edge[0], edge[1])

    return (g, order_edges)

if __name__ == '__main__':
    main()
