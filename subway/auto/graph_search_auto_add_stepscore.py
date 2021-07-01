import copy
import csv
import json
import pickle
from typing import List, Dict

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from IPython.display import display
from networkx.drawing.nx_agraph import graphviz_layout, write_dot
from pandas.core.algorithms import unique
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, rdChemReactions
from rdkit.Chem.Draw import IPythonConsole

IPythonConsole.ipython_useSVG = True
import pkg_resources

JSON_PATH = pkg_resources.resource_filename("subway", "data/auto_nodes.json")
GRAPH_PKL_PATH = pkg_resources.resource_filename("subway", "data/graph.pkl")
CSV_PATH = pkg_resources.resource_filename("subway", "data/rxn_stepscore.csv")
FULL_PROPS_PATH = pkg_resources.resource_filename(
    "subway", "data/subway_maps/full_props.csv"
)
ADJ_PATH = pkg_resources.resource_filename("subway", "data/adj_list.pkl")
PNG_PATH = pkg_resources.resource_filename("subway", "auto/visualized_graph_tree.png")


class Search:
    """
    Class for all functions needed to search through the graph to 
    produce routescores for a desired product.
    """

    def __init__(self, graph_path, json_file, csv_path, adj_path):
        """
        Parameters
        ----------
        - graph: generated graph from json file

        """
        file = open(graph_path, "rb")
        self.graph = pickle.load(file)
        file.close()
        self.json_path = json_file
        file = open(json_file, "rb")
        self.json_data = json.load(file)
        file.close()
        self.csv_path = csv_path
        file = open(adj_path, "rb")
        self.adj_list = pickle.load(file)
        file.close()

    def update_json_with_stepscore(self, rxn_node, step_score):
        for node in self.json_data:
            if node["id"] == rxn_node:
                node["step_score"] = step_score
                node["inverse_step_score"] = 1 / step_score

        with open(self.json_path, "r+") as file:
            json.dump(self.json_data, file, indent=4)

    def add_stepscore_to_csv(self, rxn_SMILES, step_score):
        with open(self.csv_path, "a", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow([rxn_SMILES, step_score])

    def route_search(self, product_smiles: str, final_scale: float):
        """ Determines routescore of all routes used to produce product_smiles

        Parameters
        ----------
        - product_smiles: SMILES of desired product

        Returns
        -------
        - columns of csv file:
            - product_name: name of the product of each route
            - route_name: Index of the route that get to the desired product
            - route_score: routescore of each route
            - visited_nodes: list of all nodes that were visited
        """
        # find product_node_index
        node_list = list(
            self.graph.nodes.data("SMILES")
        )  # compile all nodes with SMILES data
        for node in node_list:
            if node[1] == product_smiles:
                product_node_index = node[0]

        route_data = []
        child_of_product = list(self.graph.predecessors(product_node_index))
        route_name = 1
        for rxn in child_of_product:
            route_lists, visited, scores, tracker, scale = self.recursive_traversal(
                rxn, [], 0, final_scale, [], True, [], 0
            )
            for route in route_lists:
                route_dict = {
                    "product_SMILES": product_smiles,
                    "route_name": route_name,
                    "routescore": route[0] / final_scale,
                    "visited_nodes": route[1],
                }
                route_name += 1
                route_data.append(route_dict)

        # return route_data

        # column names
        # headers = [
        #     "index",
        #     "product_SMILES",
        #     "route_name",
        #     "routescore",
        #     "visited_nodes",
        # ]
        # # read csv_data (check if empty)
        # try:
        #     # file.csv is an empty csv file
        #     df = pd.read_csv(self.csv_path)
        # except pd.errors.EmptyDataError:
        #     # print("The CSV file is empty")
        #     index_count = 0
        #     with open(self.csv_path, "w", newline="") as csvfile:
        #         for route in route_data:
        #             route["index"] = index_count
        #             index_count += 1
        #         writer = csv.DictWriter(csvfile, fieldnames=headers)
        #         writer.writeheader()
        #         writer.writerows(route_data)
        # else:
        #     index_count = len(df.index)
        #     # write route_dict data to csv_data
        #     with open(self.csv_path, "a", newline="") as csvfile:
        #         for route in route_data:
        #             route["index"] = index_count
        #             index_count += 1
        #         writer = csv.DictWriter(csvfile, fieldnames=headers)
        #         writer.writerows(route_data)
        # print("Dataframe loaded successfully!!")

    def recursive_traversal(
        self,
        rxn_node: int,
        visited: List[int],
        score: float,
        final_scale: float,
        route_list: List[List],
        separate_route: bool,
        route_tracker: List,
        man_scale: float,
    ) -> List[List]:
        """ Traverses through graph to find all routes that produce product_node

        Parameters
        ----------
        - rxn_node: current reaction node
        - visited: list of nodes that have already been visited in that specific route
        - score: sum of stepscores in the recursion search
        - final_scale: final scale of reaction in mols
        - route_list: main list of routescores and visited nodes for all routes
        - separate_route: boolean for differentiating between breadth search and depth search, breadth search are separate routes, depth search is one route.
        - route_tracker: list of routes that are tracked for depth search

        Returns
        -------
        - route_list: listing containing all routescores and visited_nodes for each route
        format: [(routescore, [visited_nodes]), ...]
        """
        # NOTE: recursively traverse but calculate stepscore bottom up instead of top down
        # calculate stepscore
        child_of_rxn = list(self.graph.predecessors(rxn_node))
        step_score, scale, man_scale = Calculate().StepScore(
            self.graph, rxn_node, child_of_rxn, final_scale, man_scale
        )
        score += step_score

        # add current node to visited list
        new_visited = []
        new_visited.extend(visited)
        new_visited.append(rxn_node)

        # list of molecules that will continue the traversal
        not_last_mols = []
        for mol in child_of_rxn:
            child_of_mol = list(self.graph.predecessors(mol))
            if len(child_of_mol) != 0:
                not_last_mols.append(mol)
            new_visited.append(mol)

        # re-order not_last_mols so manual molecules are searched through first
        # NOTE: important because of scale in StepScore
        for mol in not_last_mols:
            if self.graph.nodes[mol]["Manual?"] == "Yes":
                not_last_mols.insert(0, not_last_mols.pop(not_last_mols.index(mol)))

        if len(not_last_mols) != 0:
            for mol in not_last_mols:
                child_of_not_last_mol = list(
                    self.graph.predecessors(mol)
                )  # visited one neighbour so the next node in the same level will have one less neighbour
                if len(child_of_not_last_mol) > 1:
                    for rxn_node in child_of_not_last_mol:
                        (
                            intermediate_routes,
                            inter_visited,
                            inter_score,
                            inter_tracker,
                            inter_scale,
                        ) = self.recursive_traversal(
                            rxn_node,
                            new_visited,
                            score,
                            scale,
                            route_list,
                            True,
                            route_tracker,
                            man_scale,
                        )
                elif len(child_of_not_last_mol) == 1:
                    for rxn_node in child_of_not_last_mol:
                        (
                            route_list,
                            new_visited,
                            score,
                            route_tracker,
                            scale,
                        ) = self.recursive_traversal(
                            rxn_node,
                            new_visited,
                            score,
                            scale,
                            route_list,
                            False,
                            route_tracker,
                            man_scale,
                        )
        elif not separate_route:
            if route_list == []:
                route_list.append([score, new_visited])
                route_tracker = [0]
            else:
                for route in route_tracker:
                    route_list[route][0] += step_score
                    list_difference = new_visited[len(route_list[route][1]) :]
                    route_list[route][1].extend(list_difference)
        else:
            route_list.append([score, new_visited])
            route_tracker.append(len(route_list) - 1)
        return route_list, new_visited, score, route_tracker, scale

    def run(self, full_props_path):
        """ Opens full_props.csv and runs route_search for all target molecules

        Parameters
        ----------
        - full_props_path: path to full_props.csv

        Returns
        -------
        - None
        """
        df = pd.read_csv(full_props_path)
        # find all unique smiles so routes are not repeated
        unique_smiles = df.smiles.unique()
        count = 0
        while count < len(unique_smiles):
            self.route_search(
                unique_smiles[count], (0.0001 / 6),
            )
            count += 1
        print(count)

    # NOTE: for debugging
    def route_visualizer(self, visited_nodes):
        """ Produces .png of graph of visited_nodes in a route

        Parameters
        ----------
        - visited_nodes: list of nodes that were visited

        Returns
        -------
        - route_graph: .png of graph with visited nodes
        """
        # find edges from list of visited_nodes
        visualized_graph = nx.DiGraph()
        color_map = []
        reaction_list = []
        molecule_list = []
        for visit in visited_nodes:
            visualized_graph.add_node(visit)
            for node in self.json_data:
                if node["id"] == visit:
                    if node["type"] == "reaction":
                        reaction_list.append(visit)
                    elif node["type"] == "molecule":
                        molecule_list.append(visit)

        for rxn in reaction_list:
            for mol in molecule_list:
                if [rxn, mol] in self.adj_list:
                    visualized_graph.add_edge(rxn, mol)
                if [mol, rxn] in self.adj_list:
                    visualized_graph.add_edge(mol, rxn)

        # add colour
        for node in visualized_graph.nodes:
            for data in self.json_data:
                if data["id"] == node:
                    if data["type"] == "reaction":
                        color_map.append("green")
                    elif data["type"] == "molecule":
                        color_map.append("blue")
        """
        # tree graph
        write_dot(visualized_graph, "test.dot")
        plt.title("route_visualized_tree")
        pos = graphviz_layout(visualized_graph, prog="dot")
        nx.draw(visualized_graph, node_color=color_map, with_labels=True, arrows=True)
        plt.savefig(PNG_PATH)
        """
        # draw_networkx
        plt.title("route_visualized")
        nx.draw_networkx(
            visualized_graph, node_color=color_map, with_labels=True, arrows=True
        )
        plt.savefig("visualized_graph.png")


class Calculate:
    """
    Class containing all functions required to perform any calculations for StepScore
    """

    # constants
    cH = 52.97
    cM = 17.13
    a = cH / cM
    a_ = cM / cH

    def __init__(self):
        self.calculate = None

    def CustomError(self, func, message: str):
        """Raise a custom error and print message using assert."""
        print(message)
        return func

    def draw_mols(self, smiles_list: List[str]):
        """Draw molecular structures inline.

        Parameters
        ----------
        smiles_list: List of SMILES to draw

        Returns
        -------
        Nothing
        """
        mol_list = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]
        img = Draw.MolsToGridImage(mol_list, molsPerRow=3)
        display(img)

    def stoichiometry(self, smiles: str, rxn: str) -> int:
        """Determine the number of reaction sites.

        Parameters
        ----------
        smiles:    SMILES of the molecule with reaction site substructures
        rxn:    Name of reaction to match to dictionary with substructures to count

        Returns
        -------
        mult: Number of reaction sites on the molecule
        """
        patts: dict = {
            "Suzuki": [Chem.MolFromSmarts("cBr"), Chem.MolFromSmarts("cI")],
            "deBoc": [Chem.MolFromSmiles("O=C(OC(C)(C)C)n1c2c(cc1)cccc2")],
            "BHA-H": [Chem.MolFromSmiles("[H]n1c2c(cc1)cccc2")],
            "BHA-Py": [Chem.MolFromSmiles("n1(c2nccnc2)ccc3ccccc13")],
            "SNAr-F": [Chem.MolFromSmarts("cF")],
            "SNAr-Cz": [Chem.MolFromSmarts("cn(c1c2cccc1)c3c2cccc3")],
        }
        mol = Chem.MolFromSmiles(smiles)
        matches_list: List[int] = [
            len(mol.GetSubstructMatches(substruct)) for substruct in patts[rxn]
        ]
        mult: int = sum(matches_list)
        assert mult != 0, Calculate().CustomError(
            Calculate().draw_mols([smiles]), "mult = 0!"
        )
        return mult

    def TTC(self, time_H: float, time_M: float) -> float:
        """Calculate total time cost (TTC) for the reaction.

        Parameters
        ----------
        time_H: Human labor time (in h)
        time_M: Machine labor time (in h)

        Returns
        -------
        cost_t: Result of TTC calculation
        """
        cost_t: float = np.sqrt((self.a * time_H) ** 2 + (self.a_ * time_M) ** 2)
        return cost_t

    def molecule_classifier(self, mol_node: Dict, rxn_type: str):
        """Determines molecule type by substructure and gives appropriate equivalents
        
        Parameters
        ----------
        - mol_node: node of molecule with info stored as Dict
        - rxn_type: type of reaction
        
        Returns
        -------
        - mol_info: list that contains equivalent info and name of molecule type
        """
        # figure out what substructures correspond to which smiles
        # NOTE: might have to look at block type
        if rxn_type == "wingSuzuki":
            patts: list = [
                {
                    "name": "smiles_halo-BMIDA",
                    "substructure": [
                        Chem.MolFromSmarts("Br"),
                        Chem.MolFromSmarts("I"),
                    ],
                    "eq": 1,
                },
                {
                    "name": "smiles_boronic_acid",
                    "substructure": [Chem.MolFromSmarts("cB(O)O")],
                    "eq": 3,
                },
            ]
        elif rxn_type == "pentamerSuzuki":
            patts: list = [
                {
                    "name": "smiles_BMIDA",
                    "substructure": [Chem.MolFromSmiles("B1OC(=O)CN(C)CC(=O)O1")],
                    "eq": 3,
                },
                {
                    "name": "smiles_dihalide",
                    "substructure": [
                        Chem.MolFromSmarts("cBr"),
                        Chem.MolFromSmarts("cI"),
                    ],
                    "eq": 1,
                },
            ]
        elif rxn_type == "deBoc":
            patts: list = [
                {
                    "name": "smiles_NBoc",
                    "substructure": [Chem.MolFromSmiles("CC(C)(C)OC(=O)")],
                    "eq": 1,
                }
            ]
        elif rxn_type == "BHA":
            patts: list = [
                {
                    "name": "smiles_NH",
                    "substructure": [Chem.MolFromSmarts("[nH]")],
                    "eq": 1,
                },
                {
                    "name": "smiles_halide",
                    "substructure": [Chem.MolFromSmiles("Br")],
                    "eq": 3,  # not sure if all halides
                },
            ]
        elif rxn_type == "SNAr":
            patts: list = [
                {
                    "name": "smiles_ArF",
                    "substructure": [Chem.MolFromSmarts("cF")],
                    "eq": 1,
                },
                {
                    "name": "smiles_Nu",
                    "substructure": [Chem.MolFromSmiles("c1ccc2c(c1)[nH]c1ccccc12")],
                    "eq": 2,
                },
            ]
        else:
            patts: list = []

        mol_info = []
        mol_smiles = mol_node["SMILES"]

        # make sure mol is converted from smiles or smarts
        if Chem.MolFromSmiles(mol_smiles) is None:
            mol = Chem.MolFromSmarts(mol_smiles)
        else:
            mol = Chem.MolFromSmiles(mol_smiles)

        for structure in patts:
            for substruct in structure["substructure"]:
                if mol.HasSubstructMatch(substruct):
                    mol_info.append(structure["name"])
                    mol_info.append(structure["eq"])
        if len(mol_info) == 0:
            mol_info.append("manual_or_reagent")
            mol_info.append(mol_node["eq_per_site"])

        return mol_info

    def StepScore(
        self, graph, rxn_node: int, mol_nodes: List[int], scale: float, man_scale: float
    ):
        """Perform calculations for the StepScore.

        Parameters
        ----------
        rxn_node: node of the reaction taking place
        mol_nodes: list of molecule nodes required to execute the reaction (reactants and reagents)

        Returns
        -------
        stepscore_results: tuple containing relevant information for stepscore (stepscore, yield)
        """

        # setup cost variables
        cost_time = 0
        cost_money = 0
        cost_materials = 0

        # reaction information
        rxn_json_data = graph.nodes[rxn_node]
        n_parr = float(rxn_json_data["n_parr"])
        tH = float(rxn_json_data["tH"]) / n_parr
        tM = float(rxn_json_data["tM"]) / n_parr
        yld = float(rxn_json_data["yield"])
        reaction_smiles = str(rxn_json_data["rxn_SMILES"])
        reaction_type = str(rxn_json_data["reaction_type"])
        # print(reaction_smiles, reaction_type)

        # molecule information
        # copy molecule nodes
        copy_mol_nodes = []
        for mol in mol_nodes:
            copy_mol_nodes.append(copy.copy(graph.nodes[mol]))

        # prepare molecules by adding smiles_type and equivalents to molecule_node info
        # NOTE: check if this actually modifies the graph's node: YES
        for mol_node in copy_mol_nodes:
            mol_node["smiles_type"] = Calculate().molecule_classifier(
                mol_node, reaction_type
            )[0]
            mol_node["eq_per_site"] = Calculate().molecule_classifier(
                mol_node, reaction_type
            )[1]

        # find reaction sites
        if reaction_type == "BHA":
            reaction = rdChemReactions.ReactionFromSmarts(reaction_smiles)
            product = reaction.GetProductTemplate(0)
            product_smiles = Chem.MolToSmiles(product)
            reaction_sites = self.stoichiometry(product_smiles, "BHA-Py")
        elif reaction_type == "SNAr":
            reaction = rdChemReactions.ReactionFromSmarts(reaction_smiles)
            product = reaction.GetProductTemplate(0)
            product_smiles = Chem.MolToSmiles(product)
            reaction_sites = self.stoichiometry(product_smiles, "SNAr-Cz")
        else:
            for mol_node in copy_mol_nodes:
                # print(mol_node["smiles_type"])
                if (
                    mol_node["smiles_type"] == "smiles_halo-BMIDA"
                    and reaction_type == "wingSuzuki"
                ):
                    reaction_sites = self.stoichiometry(mol_node["SMILES"], "Suzuki")
                elif (
                    mol_node["smiles_type"] == "smiles_dihalide"
                    and reaction_type == "pentamerSuzuki"
                ):
                    reaction_sites = self.stoichiometry(mol_node["SMILES"], "Suzuki")
                elif (
                    mol_node["smiles_type"] == "smiles_NBoc"
                    and reaction_type == "deBoc"
                ):
                    reaction_sites = self.stoichiometry(mol_node["SMILES"], "deBoc")

        # process molecules by multiplying reaction site with equivalents
        # reagents and necessary reactants are multiplied
        # NOTE: molecule node is changed, so changes stay
        # print(reaction_sites)
        for mol_node in copy_mol_nodes:
            if reaction_type == "wingSuzuki":
                if mol_node["smiles_type"] != "smiles_halo-BMIDA":
                    mol_node["eq_per_site"] = reaction_sites * float(
                        mol_node["eq_per_site"]
                    )
            elif reaction_type == "pentamerSuzuki":
                if mol_node["smiles_type"] != "smiles_dihalide":
                    mol_node["eq_per_site"] = reaction_sites * float(
                        mol_node["eq_per_site"]
                    )
            elif reaction_type == "deBoc":
                if mol_node["smiles_type"] != "smiles_NBoc":
                    mol_node["eq_per_site"] = reaction_sites * float(
                        mol_node["eq_per_site"]
                    )
            elif reaction_type == "BHA":
                if mol_node["smiles_type"] != "smiles_NH":
                    mol_node["eq_per_site"] = reaction_sites * float(
                        mol_node["eq_per_site"]
                    )
            elif reaction_type == "SNAr":
                if mol_node["smiles_type"] != "smiles_ArF":
                    mol_node["eq_per_site"] = reaction_sites * float(
                        mol_node["eq_per_site"]
                    )

        # calculate money and materials
        for mol_node in copy_mol_nodes:
            gmol = float(mol_node["g/mol"])
            cost_per_mol = float(mol_node["$/mol"])
            eq = float(mol_node["eq_per_site"])
            cost_money += cost_per_mol * eq
            cost_materials += gmol * eq

        cost_money = cost_money * scale
        cost_materials = cost_materials * scale
        cost_time = self.TTC(tH, tM)
        cost_money += tH * self.cH + tM * self.cM
        cost = cost_time * cost_money * cost_materials
        # print("----------")
        # print(reaction_type)
        # print("pre: ", scale)
        # print("cost:", cost_money, "time:", cost_time, "materials:", cost_materials)
        # print("stepscore: ", cost)
        # process scale
        # (NOTE: inverse design because recursive search follows backward reaction steps)
        # (NOTE: manual comes in between so order matters a lot because scale is modified in order)
        # check if there are manually made molecules
        manual_bool = False
        for mol_node in copy_mol_nodes:
            if mol_node["Manual?"] == "Yes":
                manual_bool = True

        if not manual_bool:
            for mol_node in copy_mol_nodes:
                if reaction_type == "pentamerSuzuki":
                    if mol_node["smiles_type"] == "smiles_BMIDA":
                        scale = scale * mol_node["eq_per_site"]
                elif reaction_type == "deBoc":
                    if mol_node["smiles_type"] == "smiles_NBoc":
                        scale = scale * mol_node["eq_per_site"]
                elif reaction_type == "BHA":
                    if mol_node["smiles_type"] == "smiles_NH":
                        scale = scale * mol_node["eq_per_site"]
                elif reaction_type == "SNAr":
                    if mol_node["smiles_type"] == "smiles_ArF":
                        scale = scale * mol_node["eq_per_site"]
                elif reaction_type == "manual":
                    scale = man_scale
        else:
            for mol_node in copy_mol_nodes:
                if reaction_type == "pentamerSuzuki":
                    if mol_node["smiles_type"] == "smiles_BMIDA":
                        man_scale = scale * mol_node["eq_per_site"]
                elif reaction_type == "deBoc":
                    if mol_node["smiles_type"] == "smiles_NBoc":
                        man_scale = scale * mol_node["eq_per_site"]
                elif reaction_type == "BHA":
                    if mol_node["smiles_type"] == "smiles_NH":
                        man_scale = scale * mol_node["eq_per_site"]
                elif reaction_type == "SNAr":
                    if mol_node["smiles_type"] == "smiles_ArF":
                        man_scale = scale * mol_node["eq_per_site"]
        # print("post: ", scale)
        # print("man: ", man_scale)
        Search(GRAPH_PKL_PATH, JSON_PATH, CSV_PATH, ADJ_PATH).add_stepscore_to_csv(
            rxn_json_data["rxn_SMILES"], cost
        )
        return cost, scale, man_scale


tester = Search(GRAPH_PKL_PATH, JSON_PATH, CSV_PATH, ADJ_PATH)


# route_dict = tester.route_search(
#     "Fc1cc(-c2ccsc2-c2ccc(-c3sccc3-c3ccnc(F)c3)o2)ccn1", 0.00001666666666666666666666,
# )

# print(route_dict)

# tester.route_visualizer(route_dict["visited_nodes"])


tester.run(FULL_PROPS_PATH)

# qp = Chem.AdjustQueryParameters()
# qp.makeDummiesQueries = True
# qp.adjustDegree = True
# qp.adjustDegreeFlags = Chem.ADJUST_IGNOREDUMMIES

# qm = Chem.AdjustQueryProperties(substruct, qp)

# check = mol.HasSubstructMatch(substruct)

# mol_del = Chem.DeleteSubstructs(mol, substruct)
# # mol_del_smiles = Chem.MolToSmiles(mol_del)
# img = Draw.MolToImageFile(mol_del, "mol_del.png")

