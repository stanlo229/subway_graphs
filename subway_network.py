from re import I
import networkx as nx
import json
from networkx.classes.function import all_neighbors
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from index_adj_list import *
from typing import List, Dict, Tuple
import matplotlib.pyplot as plt
import numpy as np

path_to_molecule = r'C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_graphs\subway\data\molecule_node_basic.json'
with open(path_to_molecule) as f:
    mol_json = json.load(f)

path_to_reaction = r'C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_graphs\subway\data\reaction_node_basic.json'
with open(path_to_reaction) as f2:
    rxn_json = json.load(f2)

#get multiple index positions; returns a list
def get_index_positions(list_of_elems, element):
    ''' Returns the indexes of all occurrences of give element in
    the list- listOfElements '''
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

#create graph
class Graph:

    def __init__(self):
        self.G = nx.DiGraph()
        self.color_map = []
    
    #fill part of graph that pertains to product
    #1 : create product node
    def create_product(self, SMILES: str, input_from: List[int], input_to: List[int]):
        data = list(filter(lambda x:x["SMILES"] == SMILES, mol_json))
        index = data[0]["id"]
        self.G.add_node(index, name       = data[0]["name"],
                               SMILES     = data[0]["SMILES"],
                               gmol       = data[0]["g/mol"],
                               purchasable = data[0]["purchasable"],
                               cost       = data[0]["cost"],
                               quantity   = data[0]["quantity"])
        self.color_map.append("purple")
        self.create_reaction(index, input_from, input_to)
    
    #2 : create reaction node(s)
    def create_reaction(self, mol_id: int, copy_from: List[int], copy_to: List[int]):
        index_list = get_index_positions(copy_from, mol_id)
        if len(index_list) == 0:
            return "last molecule"
        else:
            for i in index_list:
                rxn_index = copy_to[i]
                r_data = rxn_json[rxn_index]
                rxn_id = r_data["id"]
                self.G.add_node(rxn_id, name       = r_data["name"],
                                        rxn_SMILES = r_data["rxn_SMILES"],
                                        tH         = r_data["tH"],
                                        tM         = r_data["tM"],
                                        yld        = r_data["yield"],
                                        max_sites  = r_data["max_sites"])
                self.color_map.append("green")
                self.G.add_edge(mol_id, rxn_id)
            if len(index_list) != 0:
                del copy_from[index_list[0]:(index_list[-1]+1)]
                del copy_to[index_list[0]:(index_list[-1]+1)]
            self.create_molecule(rxn_index, copy_to, copy_from)
            
    #3 : create molecule node(s)
    def create_molecule(self, rxn_index: int, copy_to: List[int], copy_from: List[int],):
        r_data = rxn_json[rxn_index]
        rxn_id = r_data["id"]
        rxn_index_list = get_index_positions(copy_to, rxn_index)
        if len(rxn_index_list) == 0:
            return "last molecule"
        else:
            for i in rxn_index_list:
                mol_index = copy_from[i]
                m_data = mol_json[mol_index]
                self.G.add_node(mol_index, name       = m_data["name"],
                                           SMILES     = m_data["SMILES"],
                                           gmol       = m_data["g/mol"],
                                           purchasable = m_data["purchasable"],
                                           cost       = m_data["cost"],
                                           quantity   = m_data["quantity"])
                self.color_map.append("#2fb8f7")
                self.G.add_edge(rxn_id, mol_index)
            if len(rxn_index_list) != 0:
                del copy_from[rxn_index_list[0]:(rxn_index_list[-1]+1)]
                del copy_to[rxn_index_list[0]:(rxn_index_list[-1]+1)]
            self.create_reaction(mol_index, copy_from, copy_to)
        
    def run(self, pdt_SMILES: str, original_from: List[int], original_to: List[int]):
        copy_from = original_from.copy()
        copy_to = original_to.copy()
        self.create_product(pdt_SMILES, copy_from[::-1], copy_to[::-1])
        self.G = self.G.reverse(True)
    
    def route_search(self, node_index: int) -> Dict:
        #organizes all routes and returns routescore of each route
        #NOTE: start backtracking at first reaction
        # 5 ways: 1) search through edges, 2) only 1 reaction node connected for all SM, 
        # 3) implement own DFS, 4) implement it implicity with adjacency list, 5) directed graph
        # option 5
        #look for last reaction node, record all the nodes passed
        route_dict = {}
        main_child = list(self.G.predecessors(node_index))
        if len(main_child) == 0:
            return "reagent, pick an intermediate/product"
        else:
            counter = 1
            for main_rxn in main_child:
                for route in self.route_traversal(main_rxn, [], 0, 0.0001, []):
                    route_dict[counter] = route
                    counter += 1
        print(route_dict)
    
    def route_traversal(self, rxn_node: str, visited: List[int or str], score: float, n_target: float, main_list: List[Tuple]) -> List[Tuple]:
        #key assumption: 1 route per reaction node
        parent = list(self.G.successors(rxn_node))
        child = list(self.G.predecessors(rxn_node))
        stepscore = self.step_score(child, rxn_node, n_target)
        score += stepscore[0]
        n_target = stepscore[1]*n_target
        visited.append(rxn_node)
        visited.extend(child)
        next_mol = []
        for mol in child:
            child_child = list(self.G.predecessors(mol))
            if len(child_child) != 0:
                next_mol.append(mol)
        if len(next_mol) != 0:
            for n in next_mol:
                child_child = list(self.G.predecessors(n))
                for c in child_child:
                    return(self.route_traversal(c, visited, score, n_target, main_list))
        else:
            main_list.append((score/n_target, visited))
            return main_list

        
    def step_score(self, nodes: List[int], rxn_id: str, scale: float) -> float:
        cH = 52.97
        cM = 17.13
        a = cH/cM
        a_ = cM/cH
        sum_cost = 0
        sum_mw = 0
        ttc = 0
        sum_timecost = 0
        for i in nodes:
            gmol = self.G.nodes[i]["gmol"]
            purchasable = self.G.nodes[i]["purchasable"]
            cost = self.G.nodes[i]["cost"]
            quantity = self.G.nodes[i]["quantity"]
            sum_cost += quantity*cost*scale
            sum_mw += quantity*gmol*scale
        rxn_node = self.G.nodes[rxn_id]
        tH = rxn_node["tH"]
        tM = rxn_node["tM"]
        yld = 1
        ttc += np.sqrt((a * tH)**2 + (a_ * tM)**2) #not purely manual synthesis
        sum_timecost += (tH*cH) + (tM*cM)
        #print("ttc: ", ttc)
        #print("sum_mw: ", sum_mw) #small discrepancies to get_man_synth of routescore.py
        #print("sum_mancost: ", sum_timecost+sum_cost) #small discrepancies to get_man_synth of routescore.py

        return (ttc*(sum_timecost+sum_cost)*sum_mw, yld)


graph = Graph()
graph.run("Brc1ccc2oc(Br)cc2c1", mol_to_rxn_from, mol_to_rxn_to)
graph.route_search(7)

'''
plt.subplot(121)
nx.draw_circular(graph.G, node_color = graph.color_map, with_labels = True)
plt.show()
'''
