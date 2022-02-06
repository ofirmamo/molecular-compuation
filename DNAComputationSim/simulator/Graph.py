# This class loads txt files into usable graph data structures
#
#
# author: Jack Burns
# create date: 11/13/2017
# version: 1.0

import networkx


class Graph:

    def __init__(self):
        self.file_path_cities = 'Networks/cities.txt'
        self.file_path_nums = 'Networks/nums.txt'
        self.file_path_peterson = 'Networks/peterson.txt'
        self._citiesGraph = networkx.read_edgelist(self.file_path_cities,
                                                  create_using=networkx.DiGraph(),
                                                  nodetype=str,
                                                  data=[('to', str)])
        self._numsGraph = networkx.read_edgelist(self.file_path_nums,
                                                create_using=networkx.DiGraph(),
                                                nodetype=str,
                                                data=[('to', str)])
        self._petersonGraph = networkx.read_edgelist(self.file_path_peterson,
                                                    create_using=networkx.DiGraph(),
                                                    nodetype=str,
                                                    data=[('to', str)])
    @property
    def citiesGraph(self):
        return self._citiesGraph

    @property
    def numsGraph(self):
        return self._numsGraph

    @property
    def petersonGraph(self):
        return self._petersonGraph