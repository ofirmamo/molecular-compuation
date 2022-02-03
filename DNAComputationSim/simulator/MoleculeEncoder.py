import itertools
import random as rand

from pydna.dseq import Dseq as DnaSeq
from pydna.dseqrecord import Dseqrecord as DnaSeqRecord


class Encoder:
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
    SEQ_LEN = 20

    def synthesizeNodes(self, graph):
        return {node: "".join(rand.choices(self.base_list, k=self.SEQ_LEN))
                for node in graph.nodes()}

    def synthesizeEdges(self, graph, encoded_nodes):
        return {(u, v): encoded_nodes[u][-10:] + encoded_nodes[v][:10]
                for (u, v) in graph.edges()}

    def getSeqComplement(self, seq):
        return "".join(self.COMPLEMENTS[base] for base in seq)[::-1]

    def generateComplements(self, encoded_nodes):
        return {node: self.getSeqComplement(encoded_node)
                for node, encoded_node in encoded_nodes.items()}

    def toDnaSequence(self, encoded_nodes, encoded_edges):
        return [DnaSeqRecord(
            DnaSeq(encoded_edge, self.generateComplements(encoded_nodes)[edge[1]], ovhg=(offset := -10)),
            name=f"{edge[0]}_{edge[1]}")
                for edge, encoded_edge in encoded_edges.items()]
