import json
import pdb
import pickle
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pkg_resources
from networkx.drawing.nx_agraph import graphviz_layout, write_dot

"""
Things to keep in mind:
1. good naming conventions
2. global variables = uppercase
3. file path handling
"""

JSON_PATH = pkg_resources.resource_filename("subway", "data/auto_nodes.json")
# JSON_PATH = pkg_resources.resource_filename("subway", "gephi/auto_nodes.json")
# JSON_PATH = pkg_resources.resource_filename("subway", "data/auto_nodes_stepscore.json")
PKL_PATH = pkg_resources.resource_filename("subway", "data/adj_list.pkl")
# PKL_PATH = pkg_resources.resource_filename("subway", "gephi/adj_list_gephi.pkl")
GRAPH_PKL_PATH = pkg_resources.resource_filename("subway", "data/graph.pkl")
# GEPHI_PATH = pkg_resources.resource_filename(
#     "subway", "data/full_graph_nodes_step_score.gexf"
# )


def graph_to_pkl(graph, pkl_path):
    file = open(pkl_path, "wb")
    pickle.dump(graph, file)
    file.close()


class Graph:
    """
    Class initialized to contain json and adjacency list (nodes and edges) information, 
    and the functions for creating the graph.
    """

    def __init__(self, json_data, adj_data):
        self.graph = nx.DiGraph()
        self.color_map = []
        # create all nodes in json
        for node in json_data:
            node_id = node["id"]
            self.graph.add_nodes_from([(node_id, node)])
            if node["type"] == "molecule":
                self.color_map.append("purple")
            elif node["type"] == "reaction":
                self.color_map.append("green")
        # connect all nodes with edges from adjacency list
        # a lot less edges in graph than adj_data because duplicated edges are not added
        # duplicates come from same reactions being used but later down the route it becomes different
        self.graph.add_edges_from(adj_data)

    @staticmethod
    def create_graph_from_file(json_path, adjacency_pkl):
        # load node data
        json_data = json.load(open(json_path))

        # load adjacency list
        file = open(adjacency_pkl, "rb")
        adj_data = pickle.load(file)
        file.close()

        return Graph(json_data, adj_data)

    def create_gexf_from_file(self, gephi_path):
        nx.write_gexf(self.graph, gephi_path)


graph = Graph.create_graph_from_file(JSON_PATH, PKL_PATH)
graph_to_pkl(graph.graph, GRAPH_PKL_PATH)
# graph.create_gexf_from_file(GEPHI_PATH)

# nx.write_gexf(graph.graph, "test.gexf")
