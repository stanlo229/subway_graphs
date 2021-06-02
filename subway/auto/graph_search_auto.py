import csv
import json
import pickle
import pandas as pd
from typing import Tuple, List
import copy

from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors, Draw
from IPython.display import display
import numpy as np

IPythonConsole.ipython_useSVG = True
import pkg_resources

JSON_PATH = pkg_resources.resource_filename("subway", "data/auto_nodes.json")
GRAPH_PKL_PATH = pkg_resources.resource_filename("subway", "data/graph.pkl")
CSV_PATH = pkg_resources.resource_filename("subway", "data/results.csv")
FULL_PROPS_PATH = pkg_resources.resource_filename(
    "subway", "data/subway_maps/full_props.csv"
)


class Search:
    """
    Class for all functions needed to search through the graph to 
    produce routescores for a desired product.
    """

    def __init__(self, graph_path, json_file, csv_path):
        """
        Parameters
        ----------
        - graph: generated graph from json file

        """
        file = open(graph_path, "rb")
        self.graph = pickle.load(file)
        file.close()
        file = open(json_file, "rb")
        self.json_data = json.load(file)
        file.close()
        self.csv_path = csv_path

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
            route_lists, visited = self.recursive_traversal(
                rxn, [], 0, final_scale, [], 0
            )
            for route in route_lists:
                route_dict = {
                    "product_SMILES": product_smiles,
                    "route_name": route_name,
                    "routescore": route[0] / final_scale,
                    "visited_nodes": route[1],
                }
                route_data.append(route_dict)

        print(route_dict)

        """
        # column names
        headers = [
            "index",
            "product_SMILES",
            "route_name",
            "routescore",
            "visited_nodes",
        ]
        # read csv_data (check if empty)
        try:
            # file.csv is an empty csv file
            df = pd.read_csv(self.csv_path)
        except pd.errors.EmptyDataError:
            # print("The CSV file is empty")
            index_count = 0
            with open(self.csv_path, "w", newline="") as csvfile:
                for route in route_data:
                    route["index"] = index_count
                    index_count += 1
                writer = csv.DictWriter(csvfile, fieldnames=headers)
                writer.writeheader()
                writer.writerows(route_data)
        else:
            index_count = len(df.index)
            # write route_dict data to csv_data
            with open(self.csv_path, "a", newline="") as csvfile:
                for route in route_data:
                    route["index"] = index_count
                    index_count += 1
                writer = csv.DictWriter(csvfile, fieldnames=headers)
                writer.writerows(route_data)
            # print("Dataframe loaded successfully!!")
        """

    def recursive_traversal(
        self,
        rxn_node: int,
        visited: List[int],
        score: float,
        final_scale: float,
        route_list: List[List],
        count_neighbours: int,
    ) -> List[List]:
        """ Traverses through graph to find all routes that produce product_node

        Parameters
        ----------
        - rxn_node: current reaction node
        - visited: list of nodes that have already been visited in that specific route
        - score: sum of stepscores in the recursion search
        - final_scale: final scale of reaction in mols
        - route_list: main list of routescores and visited nodes for all routes
        - count_neighbours: number of neighbours in the current level

        Returns
        -------
        - route_list: listing containing all routescores and visited_nodes for each route
        format: [(routescore, [visited_nodes]), ...]
        """
        # NOTE: recursively traverse but calculate stepscore bottom up instead of
        # top down
        # calculate stepscore
        child_of_rxn = list(self.graph.predecessors(rxn_node))
        step_score, scale = Calculate().StepScore(
            self.graph, rxn_node, child_of_rxn, final_scale
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

        if len(not_last_mols) != 0:
            count_neighbours = len(
                not_last_mols
            )  # count the number of neighbour molecules in the same level
            for mol in not_last_mols:
                child_of_not_last_mol = list(self.graph.predecessors(mol))
                count_neighbours -= 1
                # visited one neighbour so the next node in the same level will have one less neighbour
                for rxn_node in child_of_not_last_mol:
                    intermediate_routes, new_visited = self.recursive_traversal(
                        rxn_node,
                        new_visited,
                        score,
                        scale,
                        route_list,
                        count_neighbours,
                    )
        elif count_neighbours != 0:
            pass
        else:
            route_list.append([score, new_visited])

        return route_list, new_visited

    def run(self, full_props_path):
        df = pd.read_csv(full_props_path)
        count = 0
        while count < len(df["smiles"]):
            self.route_search(df["smiles"][count], 0.000016666666)
            count += 1
            print(count)


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

    def molecule_classifier(self, mol_node: str, rxn_type: str, reaction_smiles: str):
        """Determines molecule type by substructure and gives appropriate equivalents
        
        Parameters
        ----------
        - mol_smiles: smiles of molecule
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
                    "name": "smilesB",
                    "substructure": [
                        Chem.MolFromSmarts("cBr"),
                        Chem.MolFromSmarts("cI"),
                    ],
                    "eq": 1,
                },
                {
                    "name": "smilesA",
                    "substructure": [Chem.MolFromSmarts("cB(O)O")],
                    "eq": 3,
                },
            ]
        elif rxn_type == "pentamerSuzuki":
            patts: list = [
                {
                    "name": "smilesAB",
                    "substructure": [
                        Chem.MolFromSmarts("c[B-]12OC(C[N+]1(C)CC(O2)=O)=O")
                    ],
                    "eq": 3,
                },
                {
                    "name": "smilesC",
                    "substructure": [
                        Chem.MolFromSmarts("cBr"),
                        Chem.MolFromSmarts("cI"),
                    ],
                    "eq": 1,
                },
            ]
        elif rxn_type == "deboc":
            patts: list = [
                {
                    "name": "smilesNBoc",
                    "substructure": [Chem.MolFromSmiles("O=C([N])OC(C)(C)C")],
                    "eq": 1,
                }
            ]
        elif rxn_type == "BHA":
            patts: list = [
                {
                    "name": "smilesNH",
                    "substructure": [Chem.MolFromSmiles("nH")],
                    "eq": 1,
                },
                {
                    "name": "smilesX",
                    "substructure": [Chem.MolFromSmiles("Br")],
                    "eq": 3,  # not sure if all halides
                },
            ]
        elif rxn_type == "SNAr":
            patts: list = [
                {
                    "name": "smilesAr",
                    "substructure": [Chem.MolFromSmarts("cF")],
                    "eq": 1,
                },
                {
                    "name": "smilesNu",
                    "substructure": [Chem.MolFromSmiles("c1ccc2c(c1)[nH]c1ccccc12")],
                    "eq": 2,
                },
            ]
        else:
            patts: list = []

        mol_info = []
        mol_smiles = mol_node["SMILES"]
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

    def StepScore(self, graph, rxn_node: int, mol_nodes: List[int], scale: float):
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
        # print(rxn_node)
        n_parr = float(rxn_json_data["n_parr"])
        tH = float(rxn_json_data["tH"]) / n_parr
        tM = float(rxn_json_data["tM"]) / n_parr
        yld = float(rxn_json_data["yield"])
        reaction_smiles = str(rxn_json_data["rxn_SMILES"])
        reaction_type = str(rxn_json_data["reaction_type"])
        # print(reaction_type)

        # molecule information
        # copy molecule nodes
        copy_mol_nodes = []
        for mol in mol_nodes:
            copy_mol_nodes.append(copy.copy(graph.nodes[mol]))

        # prepare molecules by adding smiles_type and equivalents to molecule_node info
        # NOTE: check if this actually modifies the graph's node: YES
        for mol_node in copy_mol_nodes:
            mol_node["smiles_type"] = Calculate().molecule_classifier(
                mol_node, reaction_type, reaction_smiles
            )[0]
            mol_node["eq_per_site"] = Calculate().molecule_classifier(
                mol_node, reaction_type, reaction_smiles
            )[1]

        # find reaction sites
        # print(reaction_type)
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
                if (
                    mol_node["smiles_type"] == "smilesB"
                    and reaction_type == "wingSuzuki"
                ):
                    reaction_sites = self.stoichiometry(mol_node["SMILES"], "Suzuki")
                elif (
                    mol_node["smiles_type"] == "smilesC"
                    and reaction_type == "pentamerSuzuki"
                ):
                    reaction_sites = self.stoichiometry(mol_node["SMILES"], "Suzuki")
                elif (
                    mol_node["smiles_type"] == "smilesNBoc" and reaction_type == "deboc"
                ):
                    reaction_sites = self.stoichiometry(mol_node["SMILES"], "deBoc")

        # process molecules by multiplying reaction site with equivalents
        # reagents and necessary reactants are multiplied
        # NOTE: molecule node is changed, so changes stay
        # print(reaction_sites)
        for mol_node in copy_mol_nodes:
            if reaction_type == "wingSuzuki":
                if mol_node["smiles_type"] != "smilesB":
                    mol_node["eq_per_site"] = reaction_sites * float(
                        mol_node["eq_per_site"]
                    )
            elif reaction_type == "pentamerSuzuki":
                if mol_node["smiles_type"] != "smilesC":
                    mol_node["eq_per_site"] = reaction_sites * float(
                        mol_node["eq_per_site"]
                    )
            elif reaction_type == "deboc":
                if mol_node["smiles_type"] != "smilesNBoc":
                    mol_node["eq_per_site"] = reaction_sites * float(
                        mol_node["eq_per_site"]
                    )
            elif reaction_type == "BHA":
                if mol_node["smiles_type"] != "smilesNH":
                    mol_node["eq_per_site"] = reaction_sites * float(
                        mol_node["eq_per_site"]
                    )
            elif reaction_type == "SNAr":
                if mol_node["smiles_type"] != "smilesAr":
                    mol_node["eq_per_site"] = reaction_sites * float(
                        mol_node["eq_per_site"]
                    )

        # calculate money and materials
        for mol_node in copy_mol_nodes:
            # print(mol_node)
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
        # print("time:", cost_time, "money:", cost_money, "materials:", cost_materials)
        # print("stepscore: ", cost)
        # process scale
        # (NOTE: inverse design because recursive search follows backward reaction steps)
        for mol_node in copy_mol_nodes:
            if reaction_type == "pentamerSuzuki":
                if mol_node["smiles_type"] == "smilesAB":
                    scale = scale * mol_node["eq_per_site"]
            elif reaction_type == "deboc":
                if mol_node["smiles_type"] == "smilesNBoc":
                    scale = scale * mol_node["eq_per_site"]
            elif reaction_type == "BHA":
                if mol_node["smiles_type"] == "smilesNH":
                    scale = scale * mol_node["eq_per_site"]
            elif reaction_type == "SNAr":
                if mol_node["smiles_type"] == "smilesAr":
                    scale = scale * mol_node["eq_per_site"]
        return cost, scale


tester = Search(GRAPH_PKL_PATH, JSON_PATH, CSV_PATH)

tester.route_search(
    "CC(C)(C)OC(=O)n1ccc2cc(-c3cccc(-c4cc(-c5cccc(-c6ccc7c(ccn7C(=O)OC(C)(C)C)c6)c5-n5c6ccccc6c6ccccc65)c5ccccc5c4)c3-n3c4ccccc4c4ccccc43)ccc21",
    0.0000166666666,
)


# tester.run(FULL_PROPS_PATH)
"""
substruct = Chem.MolFromSmarts("Br")
print(substruct)
mol = Chem.MolFromSmiles("C[N+]12CC(=O)O[B-]1(c1ccc(Br)s1)OC(=O)C2")
check = mol.HasSubstructMatch(substruct)
print(check)
"""

