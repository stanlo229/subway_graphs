import csv
import json
import pickle
from typing import Tuple, List

from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors, Draw
from IPython.display import display

IPythonConsole.ipython_useSVG = True
import pkg_resources

JSON_PATH = pkg_resources.resource_filename("subway", "basic/basic.json")
GRAPH_PKL_PATH = pkg_resources.resource_filename("subway", "basic/graph_basic.pkl")


class General:
    def load_pkl(graph_path):
        file = open(graph_path, "rb")
        data = pickle.load(file)
        file.close()
        return data


class Search:
    """
    Class for all functions needed to search through the graph to 
    produce routescores for a desired product.
    """

    def __init__(self, graph):
        """
        Parameters
        ----------
        - graph: generated graph from json file

        """
        self.graph = graph

    def route_search(self, product_smiles: str):
        """ Determines routescore of all routes used to produce product_smiles

        Parameters
        ----------
        - product_smiles: SMILES of desired product

        Returns
        -------
        - dictionary format: {route_name: (routescore, [visited nodes])}
        - columns of csv file:
            - route_name: Index of the route that get to the desired product
            - product_name: name of the product of each route
            - route_score: routescore of each route
            - scale: 
            - yield:
            - n_target:
        """
        # find product_node and product_node_index
        route_scores = {}
        product_node = None
        node_list = list(
            self.graph.nodes.data("SMILES")
        )  # compile all nodes with SMILES data
        for node in node_list:
            if node[1] == product_smiles:
                product_node = self.graph.nodes[node[0]]
                product_node_index = node[0]

        route_dict = {}
        child_of_product = list(self.graph.predecessors(product_node_index))
        route_name = 1
        for rxn in child_of_product:
            route_lists, visited = self.recursive_traversal(rxn, [], 0, 0.0001, [], 0)
            for route in route_lists:
                route_dict[route_name] = route
                route_name += 1

        print(route_dict)

    # scale going backwards will be different than going forward!
    def recursive_traversal(
        self,
        rxn_node: int,
        visited: List[int],
        score: float,
        scale: float,
        route_list: List[Tuple],
        count_neighbours: int,
    ) -> List[Tuple]:
        """ Traverses through graph to find all routes that produce molecule with node_index
        Parameters
        ----------
        - rxn_node: current reaction node
        - visited: list of visited nodes to keep track while recursive search
        - score: sum of stepscores in the recursion search
        - scale: scale of current reaction step
        - route_list: main list of routescores and visited nodes for each route
        - count_neighbours: number of neighbours in the current level

        Returns
        -------
        - route_list: list of tuples (routescores, [visited_nodes])
        """
        # calculate stepscore
        child_of_rxn = list(self.graph.predecessors(rxn_node))
        step_score, scale = Calculate.StepScore(rxn_node, child_of_rxn, scale)
        score += step_score
        scale = scale * scale

        # add current node to visited list
        new_visited = []
        new_visited.extend(visited)
        new_visited.append(rxn_node)

        # list of molecules that will continue the traversal
        not_last_mols = []

        # add molecules required for current reaction node to visited node list
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
                count_neighbours -= 1  # visited one neighbour so the next node in the same level will have one less neighbour
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
            route_list.append((score, new_visited))

        return route_list, new_visited


class Calculate:
    """
    Class containing all functions required to perform any calculations for StepScore
    """

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

    def molecule_classifier(self, mol_smiles: str):
        """Determines molecule type by substructure
        
        Parameters
        ----------
        - mol_smiles: smiles of molecule
        
        Returns
        -------
        - mol_type: the molecule's structure type
        """
        # figure out what substructures correspond to which smiles
        patts: dict = {
            "smilesA": [Chem.MolFromSmarts("cI")],
            "smilesB": [Chem.MolFromSmarts("cBr")],
            "smilesC": [Chem.MolFromSmarts("cBr")],
            "smilesNBoc": [Chem.MolFromSmiles("O=C(OC(C)(C)C)n1c2c(cc1)cccc2")],
            "smilesNH": [Chem.MolFromSmiles("nH")],
            "BHA-H": [Chem.MolFromSmiles("[H]n1c2c(cc1)cccc2")],
            "smilesBHA": [Chem.MolFromSmiles("n1(c2nccnc2)ccc3ccccc13")],
            "SNAr-F": [Chem.MolFromSmarts("cF")],
            "smilesSNAr": [Chem.MolFromSmarts("cn(c1c2cccc1)c3c2cccc3")],
        }

        mol = Chem.MolFromSmiles(mol_smiles)

    def add_on_eq(self, reaction_type: str, reaction_smiles: str):
        """Function that calculates the equivalents per reaction site for specific reactions
        
        Parameters
        ----------
        - reaction_type: reaction template for current step
        - reaction_smiles: reaction smile for current step

        Returns
        -------
        - eq: equivalents per reaction site corresponding to specific molecules
        """
        eq = 0
        return eq

    def StepScore(self, rxn_node, mol_nodes):
        """Perform calculations for the StepScore.

        Parameters
        ----------
        sm_list:        list of dictionaries corresponding to each starting material,
                        with SMILES and reaction equivalents
        product_smiles: SMILES of the desired product molecule
        target_smiles:  SMILES of the target molecule of the synthetic route
        rxn_type:       name of reaction to be carried out
        multiplier:     Multiplier (for reagents) based on number of reactive sites
        scale:          scale of the reaction in mols
        yld:            yield of the reaction
        manual:         whether the reaction if performed by a human (True) or a robot (False)

        Returns
        -------
        stepscore_results: dictionary containing relevant data for the stepscore
        stepscore_results = {
                             'StepScore': step_score,
                             'cost': cost,
                             'time': cost_time,
                             'money': cost_money,
                             'materials': cost_materials,
                             'yield': yld,
                             '# man steps': man_steps
                             }
        """
        return (1, 1)


graph = General.load_pkl(GRAPH_PKL_PATH)
search_graph = Search(graph)
a = search_graph.route_search("Brc1ccc2oc(Br)cc2c1")
