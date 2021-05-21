import networkx as nx
import json
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from typing import List, Dict, Tuple
import matplotlib.pyplot as plt
import numpy as np
from networkx.drawing.nx_agraph import write_dot, graphviz_layout
import pickle

path_to_molecule = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_graphs\subway\data\molecule_node_auto.json"
with open(path_to_molecule) as f:
    mol_json = json.load(f)

path_to_reaction = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_graphs\subway\data\reaction_node_auto.json"
with open(path_to_reaction) as f2:
    rxn_json = json.load(f2)

path_to_adj_lists = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_graphs\subway\data\auto_adj_list.pkl"
file = open(path_to_adj_lists, "rb")
adj_data = pickle.load(file)

# get multiple index positions; returns a list
def get_index_positions(list_of_elems, element):
    """ Returns the indexes of all occurrences of give element in
    the list- listOfElements """
    index_pos_list = []
    index_pos = 0
    while True:
        try:
            # Search for item in list from indexPos to the end of list
            index_pos = list_of_elems.index(element, index_pos)
            # Add the index position in list
            index_pos_list.append(index_pos)
            index_pos += 1
        except ValueError as e:
            break
    return index_pos_list


# create graph
class Graph:
    def __init__(self):
        self.G = nx.DiGraph()
        self.color_map = []

    # fill part of graph that pertains to product
    # 1 : create product node
    def create_product(self, SMILES: str, input_from: List[int], input_to: List[int]):
        data = list(filter(lambda x: x["SMILES"] == SMILES, mol_json))
        index = data[0]["id"]
        print(index)
        self.G.add_node(
            index,
            Name=data[0]["Name"],
            SMILES=data[0]["SMILES"],
            gmol=data[0]["g/mol"],
            CAD=data[0]["CAD"],
            Quantity=data[0]["Quantity"],
        )
        self.color_map.append("purple")
        self.create_reaction(index, input_from, input_to)

    # 2 : create reaction node(s)
    def create_reaction(self, mol_id: int, copy_from: List[int], copy_to: List[int]):
        index_list = get_index_positions(copy_from, mol_id)
        if len(index_list) == 0:
            return "last molecule"
        else:
            rxn_list = []
            for i in index_list:
                rxn_index = copy_to[i]
                r_data = rxn_json[rxn_index]
                rxn_id = r_data["id"]
                # check if product is made from reaction
                # get molecule smiles
                data = list(filter(lambda x: x["id"] == mol_id, mol_json))
                mol_SMILES = data[0]["SMILES"]
                can_mol_SMILES = Chem.CanonSmiles(mol_SMILES)
                # get product from rxn smiles
                rxn_SMILES = r_data["rxn_SMILES"]
                rxn = rdChemReactions.ReactionFromSmarts(rxn_SMILES)
                product = rxn.GetProductTemplate(0)
                pdt_SMILES = Chem.MolToSmiles(product)
                can_pdt_SMILES = Chem.CanonSmiles(pdt_SMILES)
                if rxn_id not in self.G and can_mol_SMILES == can_pdt_SMILES:
                    rxn_list.append(rxn_index)
                    self.G.add_node(
                        rxn_id,
                        rxn_SMILES=r_data["rxn_SMILES"],
                        tH=r_data["tH"],
                        tM=r_data["tM"],
                        yld=r_data["yld"]
                        # max_sites  = r_data["max_sites"]
                    )
                    self.color_map.append("green")
                    self.G.add_edge(mol_id, rxn_id)
            if len(index_list) != 0:
                copy_from = [j for i, j in enumerate(copy_from) if i not in index_list]
                copy_to = [j for i, j in enumerate(copy_to) if i not in index_list]
            for i in rxn_list:
                self.create_molecule(i, copy_to, copy_from)

    # 3 : create molecule node(s)
    def create_molecule(
        self, rxn_index: int, copy_to: List[int], copy_from: List[int],
    ):
        r_data = rxn_json[rxn_index]
        rxn_id = r_data["id"]
        rxn_index_list = get_index_positions(copy_to, rxn_index)
        if len(rxn_index_list) == 0:
            return "last molecule"
        else:
            for i in rxn_index_list:
                mol_id = copy_from[i]
                m_data = mol_json[mol_id]
                if mol_id not in self.G:
                    self.G.add_node(
                        mol_id,
                        Name=m_data["Name"],
                        SMILES=m_data["SMILES"],
                        gmol=m_data["g/mol"],
                        CAD=m_data["CAD"],
                        Quantity=m_data["Quantity"],
                    )
                    self.color_map.append("#2fb8f7")
                self.G.add_edge(rxn_id, mol_id)
            if len(rxn_index_list) != 0:
                copy_from = [
                    j for i, j in enumerate(copy_from) if i not in rxn_index_list
                ]
                copy_to = [j for i, j in enumerate(copy_to) if i not in rxn_index_list]
            self.create_reaction(mol_id, copy_from, copy_to)

    def run(self, pdt_SMILES: str, original_from: List[int], original_to: List[int]):
        copy_from = original_from.copy()
        copy_to = original_to.copy()
        self.create_product(pdt_SMILES, copy_from[::-1], copy_to[::-1])
        self.G = self.G.reverse(True)


graph = Graph()
graph.run(
    "c1cnc(-c2ccc3c(ccn3-c3cnccn3)c2)c(-c2cc(-n3c4ccccc4c4ccccc43)c(-c3cccnc3-c3ccc4c(ccn4-c4cnccn4)c3)cc2-n2c3ccccc3c3ccccc32)c1",
    adj_data[0],
    adj_data[1],
)

print(graph.G.number_of_nodes())
print(graph.G.number_of_edges())

with open(
    r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_graphs\subway\data\graph_auto.pkl",
    "wb",
) as output:
    pickle.dump(graph, output)


write_dot(graph.G, "test.dot")
plt.title("subway_graph")
pos = graphviz_layout(graph.G, prog="dot")
nx.draw(graph.G, pos, node_color=graph.color_map, with_labels=True, arrows=True)
plt.savefig("directed_graph_auto.png")

"""
plt.subplot(121)
nx.draw_circular(graph.G, node_color = graph.color_map, with_labels = True)
plt.show()
"""

