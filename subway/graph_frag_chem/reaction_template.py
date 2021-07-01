from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import Descriptors
import json
import pkg_resources
from typing import List
import pickle
import copy

from rdkit.Chem.Draw.IPythonConsole import _GetSubstructMatches
from subway.graph_frag_chem.graph_frag_chem import ProductMaker

JSON_PATH = pkg_resources.resource_filename(
    "subway", "data/reaction_template_nodes.json"
)

ADJ_PATH = pkg_resources.resource_filename("subway", "data/reaction_template_adj.pkl")


class ReactionTemplates:
    """ 
    Class containing reaction templates that will create reaction nodes and product/intermediate nodes from molecule nodes
    """

    def __init__(self, json_path, adjacency_path):
        self.json_path = json_path
        self.adj_path = adjacency_path
        file = open(json_path, "rb")
        self.json_data = json.load(file)
        file.close()
        file = open(adjacency_path, "rb")
        self.adj_list = pickle.load(file)
        file.close()

    def write_json(self, new_data, filename):
        """Function that writes new data to json file (list)
        """
        with open(filename, "r+") as file:
            # First we load existing data into a dict.
            file_data = json.load(file)
            # Join new_dat3a with file_data
            file_data.append(new_data)
            # Sets file's current position at offset.
            file.seek(0)
            # convert back to json.
            json.dump(file_data, file, indent=4)

    def reload_json(self):
        """Function that will re-load json data to ensure most updated json file is used.
        """
        file = open(self.json_path, "rb")
        self.json_data = json.load(file)
        file.close()

    def getNodes(self, smiles_list: List[str]):
        """Function that returns list of nodes from list of SMILES.
        Parameters
        ----------
        smiles_list: list of smiles

        Returns
        -------
        node_list: list of nodes from json

        """
        node_list = []
        for node in self.json_data:
            for smile in smiles_list:
                if node["type"] == "molecule":
                    if node["SMILES"] == smile:
                        node_list.append(node)

        return node_list

    def getNextId(self):
        """Function that returns (largest id + 1) in the json file.
        Parameters
        ----------
        None

        Returns
        -------
        next_id = id number of next node that needs to be created

        """
        id = 0
        for node in self.json_data:
            if node["id"] > id:
                id = node["id"]
        return id + 1

    def molecule_node_maker(
        self,
        mol_obj,
        cost_per_mol: float,
        cad: float,
        eq_per_site: float,
        name: str,
        manual: str,
    ):
        """ Function that will create molecule node and update json file.

        Parameters
        ----------
        mol_obj: rdkit molecule
        cost_per_mol: cost of molecule in CAD per mole
        cad: the cost in CAD of the molecule
        grams: the amount of the molecule in grams
        molar_mass: the molar mass of the molecule
        eq_per_site: number of equivalents per reaction site
        name: name of molecule (not consistent naming), intermediates/products don't have a name
        # rxn_sites: number of reaction sites on limiting reagent

        Returns
        -------
        None
        """
        molecule_node = {}
        molecule_node["Block_type"] = "-"
        molecule_node["Block_number"] = "0"
        molecule_node["Name"] = name
        molecule_node["SMILES"] = Chem.MolToSmiles(mol_obj)

        # compute g/mol with RDKit
        gmol = Descriptors.MolWt(mol_obj)
        molecule_node["g/mol"] = gmol

        molecule_node["CAD"] = cad
        molecule_node["$/mol"] = cost_per_mol
        molecule_node["eq_per_site"] = eq_per_site
        molecule_node["Manual?"] = manual
        molecule_node["id"] = self.getNextId()
        molecule_node["type"] = "molecule"

        # need to modify so it adds smoothly
        self.write_json(molecule_node, self.json_path)
        self.reload_json()

        return molecule_node["id"]

    def reaction_node_maker(
        self,
        reaction_type: str,
        reaction_smile: str,
        tH: float,
        tM: float,
        n_parr: int,
        yld: float,
        scale: float,
        # rxn_sites: int,
    ):
        """ Function that will create reaction node and update json file.

        Parameters
        ----------
        reaction_type: type of reaction
        reaction_smile: SMILE of reaction
        tH: human time
        tM: machine time
        n_parr: reactions performed in parallel
        yld: yield of reaction
        scale: scale of reaction
        # rxn_sites: number of reaction sites on limiting reagent

        Returns
        -------
        None
        """
        reaction_node = {}
        reaction_node["rxn_SMILES"] = reaction_smile
        reaction_node["t_H"] = tH
        reaction_node["t_M"] = tM
        reaction_node["n_parr"] = n_parr
        reaction_node["yield"] = yld
        reaction_node["scale"] = scale
        reaction_node["reaction_type"] = reaction_type
        reaction_node["type"] = "reaction"
        reaction_node["id"] = self.getNextId()

        self.write_json(reaction_node, self.json_path)
        self.reload_json()

        return reaction_node["id"]

    def adjacency_list_builder(
        self, reactants: List[str], reagents: List[str], reaction: int, product: int
    ):
        """ Function that will find ids of all necessary compents and add to the adjacency list.

        Parameters
        ----------
        reactants: list of reactant SMILES
        reagents: list of reagent SMILES
        reaction: id of reaction node
        product: id of product node

        Returns
        -------
        None
        """
        # get reactant and reagent ids and connect to reaction node
        for node in self.json_data:
            if node["type"] == "molecule":
                for rct in reactants:
                    if node["SMILES"] == rct:
                        self.adj_list.append([node["id"], reaction])
                for rgt in reagents:
                    if node["SMILES"] == rgt:
                        self.adj_list.append([node["id"], reaction])

        # add product node and reaction node connection
        file = open(self.adj_path, "wb")
        self.adj_list.append([reaction, product])
        pickle.dump(self.adj_list, file)
        file.close()

    def wingSuzuki(self, reactants_smiles: List[str]):
        """ Function that will contain reaction information for any wingSuzuki reaction.

        Parameters
        ----------
        reactants_smiles: list of the reactant's SMILES

        Returns
        -------
        reaction_dict = dictionary containing all information for reaction node
        """
        # reaction-specific information
        reagent_smiles = [
            "CC(C)C1=CC(=C(C(=C1)C(C)C)C2=CC(=CC=C2)P(C3CCCCC3)C4CCCCC4)C(C)C.C1=CC=C([C-]=C1)C2=CC=CC=C2N.Cl[Pd+]",
            "[O-]P(=O)([O-])[O-].[K+].[K+].[K+]",
        ]
        patts: list = [
            {
                "name": "smiles_halo-BMIDA",
                "substructure": [Chem.MolFromSmarts("cBr"), Chem.MolFromSmarts("cI"),],
                "eq": 1,
            },
            {
                "name": "smiles_boronic_acid",
                "substructure": [
                    Chem.MolFromSmarts("cB(O)O")
                ],  # explicitly define the H's to avoid MIDA boronate / BMIDA
                "eq": 3,
            },
        ]
        t_H = 2.5
        t_M = 28
        n_parr = 48
        yld = 1
        scale = 0.0001

        # get molecule nodes from json
        mol_node_list = self.getNodes(reactants_smiles)

        # process reactants to create reaction_node
        # find limiting reagent by calculating reaction sites
        # NOTE: for now we don't need this, since we're just trying to replicate
        rxn_sites = 0
        for rct in reactants_smiles:
            mol = Chem.MolFromSmiles(rct)
            for structure in patts:
                for substruct in structure["substructure"]:
                    if len(mol.GetSubstructMatches(substruct)) > rxn_sites:
                        rxn_sites = len(mol.GetSubstructMatches(substruct))
                        limiting_rct_smiles = rct

        # creating reaction node
        reaction_node = {}
        reaction_smile = ""
        # creating Reaction SMILE
        # add reactants
        index = 0
        while index < len(reactants_smiles):
            reaction_smile += reactants_smiles[index]
            index += 1
            if index != len(reactants_smiles):
                reaction_smile += "."
        reaction_smile += ">"

        # add reagents
        index = 0
        while index < len(reagent_smiles):
            reaction_smile += reagent_smiles[index]
            index += 1
            if index != len(reagent_smiles):
                reaction_smile += "."
        reaction_smile += ">"

        # add product (import graph_frag_chem / graph_ibm_rxn)
        product_list = ProductMaker().rct_to_pdt(mol_node_list, "wingSuzuki")

        product_smile = Chem.MolToSmiles(product_list[0])
        # add product molecule node and returns id of product node
        # check if node exists in json file
        duplicate_node = False
        for node in self.json_data:
            if node["type"] == "molecule":
                if node["SMILES"] == product_smile:
                    product_id = node["id"]
                    duplicate_node = True

        if not duplicate_node:
            product_id = self.molecule_node_maker(
                product_list[0], 0.0, 0.0, 0.0, "-", ""
            )

        # assume 1 product
        reaction_smile += product_smile

        # add reaction node to json and returns id of reaction node
        # check if reaction node exists in json file
        duplicate_rxn_node = False
        for node in self.json_data:
            if node["type"] == "reaction":
                if node["rxn_SMILES"] == reaction_smile:
                    reaction_id = node["id"]
                    duplicate_rxn_node = True

        if not duplicate_rxn_node:
            reaction_id = self.reaction_node_maker(
                "wingSuzuki", reaction_smile, t_H, t_M, n_parr, yld, scale
            )

        # add to adjacency list
        self.adjacency_list_builder(
            reactants_smiles, reagent_smiles, reaction_id, product_id
        )

    def pentamerSuzuki(self, reactants_smiles: List[str]):
        """ Function that will contain reaction information for any pentamerSuzuki reaction.
        
        Parameters
        ----------
        reactants_smiles: list of the reactant's SMILES

        Returns
        -------
        None
        """
        # reaction-specific information
        reagent_smiles = [
            "CC(C)C1=CC(=C(C(=C1)C(C)C)C2=CC(=CC=C2)P(C3CCCCC3)C4CCCCC4)C(C)C.C1=CC=C([C-]=C1)C2=CC=CC=C2N.Cl[Pd+]",
            "[O-]P(=O)([O-])[O-].[K+].[K+].[K+]",
        ]
        patts: list = [
            {
                "name": "smiles_BMIDA",
                "substructure": [Chem.MolFromSmiles("B1OC(=O)CN(C)CC(=O)O1")],
                "eq": 3,
            },
            {
                "name": "smiles_dihalide",
                "substructure": [Chem.MolFromSmarts("cBr"), Chem.MolFromSmarts("cI"),],
                "eq": 1,
            },
        ]
        t_H = 2.5
        t_M = 28
        n_parr = 48
        yld = 1
        scale = 0.0001

        # get molecule nodes from json
        mol_node_list = self.getNodes(reactants_smiles)

        # process reactants to create reaction_node
        # find limiting reagent by calculating reaction sites
        # NOTE: for now we don't need this, since we're just trying to replicate
        rxn_sites = 0
        for rct in reactants_smiles:
            mol = Chem.MolFromSmiles(rct)
            for structure in patts:
                for substruct in structure["substructure"]:
                    if len(mol.GetSubstructMatches(substruct)) > rxn_sites:
                        rxn_sites = len(mol.GetSubstructMatches(substruct))
                        limiting_rct_smiles = rct

        # creating reaction node
        reaction_node = {}
        reaction_smile = ""
        # creating Reaction SMILE
        # add reactants
        index = 0
        while index < len(reactants_smiles):
            reaction_smile += reactants_smiles[index]
            index += 1
            if index != len(reactants_smiles):
                reaction_smile += "."
        reaction_smile += ">"

        # add reagents
        index = 0
        while index < len(reagent_smiles):
            reaction_smile += reagent_smiles[index]
            index += 1
            if index != len(reagent_smiles):
                reaction_smile += "."
        reaction_smile += ">"

        # add product (import graph_frag_chem / graph_ibm_rxn)
        product_list = ProductMaker().rct_to_pdt(mol_node_list, "pentamerSuzuki")

        product_smile = Chem.MolToSmiles(product_list[0])
        # add product molecule node and returns id of product node
        # check if node exists in json file
        duplicate_node = False
        for node in self.json_data:
            if node["type"] == "molecule":
                if node["SMILES"] == product_smile:
                    product_id = node["id"]
                    duplicate_node = True

        if not duplicate_node:
            product_id = self.molecule_node_maker(
                product_list[0], 0.0, 0.0, 0.0, "-", ""
            )

        # assume 1 product
        reaction_smile += product_smile

        # add reaction node to json and returns id of reaction node
        # check if reaction node exists in json file
        duplicate_rxn_node = False
        for node in self.json_data:
            if node["type"] == "reaction":
                if node["rxn_SMILES"] == reaction_smile:
                    reaction_id = node["id"]
                    duplicate_rxn_node = True

        if not duplicate_rxn_node:
            reaction_id = self.reaction_node_maker(
                "pentamerSuzuki", reaction_smile, t_H, t_M, n_parr, yld, scale
            )

        # add to adjacency list
        self.adjacency_list_builder(
            reactants_smiles, reagent_smiles, reaction_id, product_id
        )

    def deBoc(self, reactants_smiles: List[str]):
        """ Function that will contain reaction information for any deBoc reaction.
        
        Parameters
        ----------
        reactants_smiles: list of the reactant's SMILES

        Returns
        -------
        None
        """
        # reaction-specific information
        reagent_smiles = ["C(=O)([O-])[O-].[K+].[K+]"]
        patts: list = [
            {
                "name": "smiles_NBoc",
                "substructure": [Chem.MolFromSmiles("CC(C)(C)OC(=O)")],
                "eq": 1,
            }
        ]
        t_H = 6.5
        t_M = 0
        n_parr = 1
        yld = 1
        scale = 0.0001

        # get molecule nodes from json
        mol_node_list = self.getNodes(reactants_smiles)

        # process reactants to create reaction_node
        # find limiting reagent by calculating reaction sites
        # NOTE: for now we don't need this, since we're just trying to replicate
        rxn_sites = 0
        for rct in reactants_smiles:
            mol = Chem.MolFromSmiles(rct)
            for structure in patts:
                for substruct in structure["substructure"]:
                    if len(mol.GetSubstructMatches(substruct)) > rxn_sites:
                        rxn_sites = len(mol.GetSubstructMatches(substruct))
                        limiting_rct_smiles = rct

        # creating reaction node
        reaction_node = {}
        reaction_smile = ""
        # creating Reaction SMILE
        # add reactants
        index = 0
        while index < len(reactants_smiles):
            reaction_smile += reactants_smiles[index]
            index += 1
            if index != len(reactants_smiles):
                reaction_smile += "."
        reaction_smile += ">"

        # add reagents
        index = 0
        while index < len(reagent_smiles):
            reaction_smile += reagent_smiles[index]
            index += 1
            if index != len(reagent_smiles):
                reaction_smile += "."
        reaction_smile += ">"

        # add product (import graph_frag_chem / graph_ibm_rxn)
        product_list = ProductMaker().rct_to_pdt(mol_node_list, "deBoc")

        product_smile = Chem.MolToSmiles(product_list[0])
        # add product molecule node and returns id of product node
        # check if node exists in json file
        duplicate_node = False
        for node in self.json_data:
            if node["type"] == "molecule":
                if node["SMILES"] == product_smile:
                    product_id = node["id"]
                    duplicate_node = True

        if not duplicate_node:
            product_id = self.molecule_node_maker(
                product_list[0], 0.0, 0.0, 0.0, "-", ""
            )

        # assume 1 product
        reaction_smile += product_smile

        # add reaction node to json and returns id of reaction node
        # check if reaction node exists in json file
        duplicate_rxn_node = False
        for node in self.json_data:
            if node["type"] == "reaction":
                if node["rxn_SMILES"] == reaction_smile:
                    reaction_id = node["id"]
                    duplicate_rxn_node = True

        if not duplicate_rxn_node:
            reaction_id = self.reaction_node_maker(
                "deBoc", reaction_smile, t_H, t_M, n_parr, yld, scale
            )

        # add to adjacency list
        self.adjacency_list_builder(
            reactants_smiles, reagent_smiles, reaction_id, product_id
        )

    def BHA(self, reactants_smiles: List[str]):
        """ Function that will contain reaction information for any BHA reaction.
        
        Parameters
        ----------
        reactants_smiles: list of the reactant's SMILES

        Returns
        -------
        None
        """
        # reaction-specific information
        reagent_smiles = [
            "C1=CC=C(C=C1)C=CC(=O)C=CC2=CC=CC=C2.C1=CC=C(C=C1)C=CC(=O)C=CC2=CC=CC=C2.C1=CC=C(C=C1)C=CC(=O)C=CC2=CC=CC=C2.[Pd].[Pd]",
            "CN(C)C1=CC=CC=C1C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4",
            "CC(C)(C)[O-].[Na+]",
        ]
        patts: list = [
            {
                "name": "smiles_NH",
                "substructure": [Chem.MolFromSmarts("[nH]")],
                "eq": 1,
            },
            {
                "name": "smiles_halide",
                "substructure": [Chem.MolFromSmiles("Br")],
                "eq": 3,
            },
        ]
        t_H = 6
        t_M = 0
        n_parr = 1
        yld = 1
        scale = 0.0001

        # get molecule nodes from json
        mol_node_list = self.getNodes(reactants_smiles)

        # process reactants to create reaction_node
        # find limiting reagent by calculating reaction sites
        # NOTE: for now we don't need this, since we're just trying to replicate
        rxn_sites = 0
        for rct in reactants_smiles:
            mol = Chem.MolFromSmiles(rct)
            for structure in patts:
                for substruct in structure["substructure"]:
                    if len(mol.GetSubstructMatches(substruct)) > rxn_sites:
                        rxn_sites = len(mol.GetSubstructMatches(substruct))
                        limiting_rct_smiles = rct

        # creating reaction node
        reaction_node = {}
        reaction_smile = ""
        # creating Reaction SMILE
        # add reactants
        index = 0
        while index < len(reactants_smiles):
            reaction_smile += reactants_smiles[index]
            index += 1
            if index != len(reactants_smiles):
                reaction_smile += "."
        reaction_smile += ">"

        # add reagents
        index = 0
        while index < len(reagent_smiles):
            reaction_smile += reagent_smiles[index]
            index += 1
            if index != len(reagent_smiles):
                reaction_smile += "."
        reaction_smile += ">"

        # add product (import graph_frag_chem / graph_ibm_rxn)
        product_list = ProductMaker().rct_to_pdt(mol_node_list, "BHA")

        product_smile = Chem.MolToSmiles(product_list[0])
        # add product molecule node and returns id of product node
        # check if node exists in json file
        duplicate_node = False
        for node in self.json_data:
            if node["type"] == "molecule":
                if node["SMILES"] == product_smile:
                    product_id = node["id"]
                    duplicate_node = True

        if not duplicate_node:
            product_id = self.molecule_node_maker(
                product_list[0], 0.0, 0.0, 0.0, "-", ""
            )

        # assume 1 product
        reaction_smile += product_smile

        # add reaction node to json and returns id of reaction node
        # check if reaction node exists in json file
        duplicate_rxn_node = False
        for node in self.json_data:
            if node["type"] == "reaction":
                if node["rxn_SMILES"] == reaction_smile:
                    reaction_id = node["id"]
                    duplicate_rxn_node = True

        if not duplicate_rxn_node:
            reaction_id = self.reaction_node_maker(
                "BHA", reaction_smile, t_H, t_M, n_parr, yld, scale
            )

        # add to adjacency list
        self.adjacency_list_builder(
            reactants_smiles, reagent_smiles, reaction_id, product_id
        )

    def SNAr(self, reactants_smiles: List[str]):
        """ Function that will contain reaction information for any SNAr reaction.
        
        Parameters
        ----------
        reactants_smiles: list of the reactant's SMILES

        Returns
        -------
        None
        """
        # reaction-specific information
        reagent_smiles = ["C(=O)([O-])[O-].[Cs+].[Cs+]"]
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
        t_H = 6.5
        t_M = 0
        n_parr = 1
        yld = 1
        scale = 0.0001

        # get molecule nodes from json
        mol_node_list = self.getNodes(reactants_smiles)

        # process reactants to create reaction_node
        # find limiting reagent by calculating reaction sites
        # NOTE: for now we don't need this, since we're just trying to replicate
        rxn_sites = 0
        for rct in reactants_smiles:
            mol = Chem.MolFromSmiles(rct)
            for structure in patts:
                for substruct in structure["substructure"]:
                    if len(mol.GetSubstructMatches(substruct)) > rxn_sites:
                        rxn_sites = len(mol.GetSubstructMatches(substruct))
                        limiting_rct_smiles = rct

        # creating reaction node
        reaction_node = {}
        reaction_smile = ""
        # creating Reaction SMILE
        # add reactants
        index = 0
        while index < len(reactants_smiles):
            reaction_smile += reactants_smiles[index]
            index += 1
            if index != len(reactants_smiles):
                reaction_smile += "."
        reaction_smile += ">"

        # add reagents
        index = 0
        while index < len(reagent_smiles):
            reaction_smile += reagent_smiles[index]
            index += 1
            if index != len(reagent_smiles):
                reaction_smile += "."
        reaction_smile += ">"

        # add product (import graph_frag_chem / graph_ibm_rxn)
        product_list = ProductMaker().rct_to_pdt(mol_node_list, "SNAr")
        # print(product_list)
        # multiple products, for loop is required
        for product in product_list:
            # copy of reaction smile so original doesn't get modified
            reaction_smile_copy = copy.copy(reaction_smile)
            product_smile = Chem.MolToSmiles(product)
            # print(product_smile)
            # add product molecule node and returns id of product node
            # check if node exists in json file
            duplicate_node = False
            for node in self.json_data:
                if node["type"] == "molecule":
                    if node["SMILES"] == product_smile:
                        product_id = node["id"]
                        duplicate_node = True
            if not duplicate_node:
                product_id = self.molecule_node_maker(product, 0.0, 0.0, 0.0, "-", "")

            # add product smile to reaction
            reaction_smile_copy += product_smile
            # add reaction node to json and returns id of reaction node
            # check if reaction node exists in json file
            duplicate_rxn_node = False
            for node in self.json_data:
                if node["type"] == "reaction":
                    if node["rxn_SMILES"] == reaction_smile_copy:
                        reaction_id = node["id"]
                        duplicate_rxn_node = True
            if not duplicate_rxn_node:
                reaction_id = self.reaction_node_maker(
                    "SNAr", reaction_smile_copy, t_H, t_M, n_parr, yld, scale
                )
            # add to adjacency list
            self.adjacency_list_builder(
                reactants_smiles, reagent_smiles, reaction_id, product_id
            )


# rxn_template = ReactionTemplates(JSON_PATH, ADJ_PATH)
# rxn_template.wingSuzuki(["OB(O)c1ccco1", "C[N+]12CC(=O)O[B-]1(c1cccnc1Br)OC(=O)C2"])


# tester for molecular weight
# mol = Chem.MolFromSmarts("CC(=CBr)c1([N+](C)(CC=O)CC=O)ccc2ncccc2c1")
# print(Descriptors.MolWt(mol))
