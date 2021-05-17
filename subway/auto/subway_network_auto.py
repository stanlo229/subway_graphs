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
from networkx.drawing.nx_agraph import write_dot, graphviz_layout
import pickle

path_to_molecule = r'C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_graphs\subway\data\molecule_node_auto.json'
with open(path_to_molecule) as f:
    mol_json = json.load(f)

path_to_reaction = r'C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_graphs\subway\data\reaction_node_auto.json'
with open(path_to_reaction) as f2:
    rxn_json = json.load(f2)

path_to_adj_lists = r'C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_graphs\auto_adj_list.pkl'
file = open(path_to_adj_lists, "rb")
adj_data = pickle.load(file)

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
        data = list(filter(lambda x:x['SMILES'] == SMILES, mol_json))
        index = data[0]["id"]
        self.G.add_node(index, Name        = data[0]["Name"],
                               SMILES      = data[0]["SMILES"],
                               gmol        = data[0]["g/mol"],
                               CAD         = data[0]["CAD"],
                               Quantity    = data[0]["Quantity"])
        self.color_map.append("purple")
        self.create_reaction(index, input_from, input_to)
    
    #2 : create reaction node(s)
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
                if rxn_id not in self.G:
                    rxn_list.append(rxn_index)
                    self.G.add_node(rxn_id, rxn_SMILES = r_data["rxn_SMILES"],
                                            tH         = r_data["tH"],
                                            tM         = r_data["tM"],
                                            yld        = r_data["yld"]
                                            #max_sites  = r_data["max_sites"]
                                            )
                    self.color_map.append("green")
                    self.G.add_edge(mol_id, rxn_id)
            if len(index_list) != 0:
                copy_from = [j for i, j in enumerate(copy_from) if i not in index_list]
                copy_to = [j for i, j in enumerate(copy_to) if i not in index_list]
            for i in rxn_list:
                self.create_molecule(i, copy_to, copy_from)
            
    #3 : create molecule node(s)
    def create_molecule(self, rxn_index: int, copy_to: List[int], copy_from: List[int],):
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
                    self.G.add_node(mol_id, Name        = m_data["Name"],
                                            SMILES      = m_data["SMILES"],
                                            gmol        = m_data["g/mol"],
                                            CAD         = m_data["CAD"],
                                            Quantity    = m_data["Quantity"])
                    self.color_map.append("#2fb8f7")
                self.G.add_edge(rxn_id, mol_id)
            if len(rxn_index_list) != 0:
                copy_from = [j for i, j in enumerate(copy_from) if i not in rxn_index_list]
                copy_to = [j for i, j in enumerate(copy_to) if i not in rxn_index_list]
            self.create_reaction(mol_id, copy_from, copy_to)
        
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
                route_lists = self.route_traversal(main_rxn, [], 0, 0.0001, [])
                for route in route_lists:
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
        new_visited = []
        new_visited.extend(visited)
        new_visited.append(rxn_node)
        next_mol = []
        for mol in child:
            child_child = list(self.G.predecessors(mol))
            if len(child_child) != 0:
                next_mol.append(mol)
            else:
                new_visited.append(mol)
        if len(next_mol) != 0:
            for n in next_mol:
                n_child_child = list(self.G.predecessors(n))
                new_visited.append(n)
                for c in n_child_child:
                    self.route_traversal(c, new_visited, score, n_target, main_list)
        else:
            main_list.append((score/n_target, new_visited))
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
            gmol = float(self.G.nodes[i]["gmol"])
            cost = float(self.G.nodes[i]["CAD"])
            quantity = float(self.G.nodes[i]["Quantity"])
            sum_cost += quantity*cost*scale
            sum_mw += quantity*gmol*scale
        rxn_node = self.G.nodes[rxn_id]
        tH = float(rxn_node["tH"])
        tM = float(rxn_node["tM"])
        yld = 1
        ttc += np.sqrt((a * tH)**2 + (a_ * tM)**2) #not purely manual synthesis
        sum_timecost += (tH*cH) + (tM*cM)
        #print("ttc: ", ttc)
        #print("sum_mw: ", sum_mw) #small discrepancies to get_man_synth of routescore.py
        #print("sum_mancost: ", sum_timecost+sum_cost) #small discrepancies to get_man_synth of routescore.py

        return (ttc*(sum_timecost+sum_cost)*sum_mw, yld)



graph = Graph()
graph.run("CC(C)(C)OC(=O)n1ccc2cc(-c3cccc(F)c3[B-]34OC(=O)C[N+]3(C)CC(=O)O4)ccc21", adj_data[0], adj_data[1])
graph.route_search(51)

write_dot(graph.G,'test.dot')
plt.title("subway_graph")
pos = graphviz_layout(graph.G, prog = 'dot')
nx.draw(graph.G, pos, node_color = graph.color_map, with_labels = True, arrows = True)
plt.savefig('directed_graph_auto.png')

'''
plt.subplot(121)
nx.draw_circular(graph.G, node_color = graph.color_map, with_labels = True)
plt.show()
'''

