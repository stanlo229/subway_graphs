from re import I
import pickle
from subway.data.index_adj_list import *
from typing import List, Dict, Tuple
import matplotlib.pyplot as plt
import numpy as np
from networkx.drawing.nx_agraph import write_dot, graphviz_layout
import pdb
from subway_network_build import Graph

# pdb.set_trace() for stopping while debugging

# graph pickle
graph_pkl = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_graphs\subway\data\graph.pkl"


def route_search(graph, node_index: int) -> Dict:
    # organizes all routes and returns routescore of each route
    # NOTE:look for last reaction node, record all the nodes passed
    route_dict = {}
    main_child = list(graph.G.predecessors(node_index))
    if len(main_child) == 0:
        return "reagent, pick an intermediate/product"
    else:
        counter = 1
        for main_rxn in main_child:
            route_lists = route_traversal(graph, main_rxn, [], 0, 0.0001, [])
            for route in route_lists:
                route_dict[counter] = route
                counter += 1
    print(route_dict)


def route_traversal(
    graph,
    rxn_node: str,
    visited: List[int or str],
    score: float,
    n_target: float,
    main_list: List[Tuple],
) -> List[Tuple]:
    # key assumption: 1 product per reaction node
    parent = list(graph.G.successors(rxn_node))
    child = list(graph.G.predecessors(rxn_node))
    stepscore = step_score(graph, child, rxn_node, n_target)
    score += stepscore[0]
    n_target = stepscore[1] * n_target
    new_visited = []
    new_visited.extend(visited)
    new_visited.append(rxn_node)
    next_mol = []
    for mol in child:
        child_child = list(graph.G.predecessors(mol))
        if len(child_child) != 0:
            next_mol.append(mol)
        else:
            new_visited.append(mol)
    if len(next_mol) != 0:
        for n in next_mol:
            n_child_child = list(graph.G.predecessors(n))
            new_visited.append(n)
            for c in n_child_child:
                route_traversal(graph, c, new_visited, score, n_target, main_list)
    else:
        main_list.append((score / n_target, new_visited))
    return main_list


def step_score(graph, nodes: List[int], rxn_id: str, scale: float) -> float:
    cH = 52.97
    cM = 17.13
    a = cH / cM
    a_ = cM / cH
    sum_cost = 0
    sum_mw = 0
    ttc = 0
    sum_timecost = 0
    for i in nodes:
        gmol = graph.G.nodes[i]["gmol"]
        purchasable = graph.G.nodes[i]["purchasable"]
        cost = graph.G.nodes[i]["CAD"]
        quantity = graph.G.nodes[i]["Quantity"]
        sum_cost += quantity * cost * scale
        sum_mw += quantity * gmol * scale
    rxn_node = graph.G.nodes[rxn_id]
    tH = rxn_node["tH"]
    tM = rxn_node["tM"]
    yld = 1
    ttc += np.sqrt((a * tH) ** 2 + (a_ * tM) ** 2)  # not purely manual synthesis
    sum_timecost += (tH * cH) + (tM * cM)
    # print("ttc: ", ttc)
    # print("sum_mw: ", sum_mw) #small discrepancies to get_man_synth of routescore.py
    # print("sum_mancost: ", sum_timecost+sum_cost) #small discrepancies to get_man_synth of routescore.py

    return (ttc * (sum_timecost + sum_cost) * sum_mw, yld)


file = open(graph_pkl, "rb")
graph = pickle.load(file)
file.close()
route_search(graph, 7)
