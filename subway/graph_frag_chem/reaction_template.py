from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import Descriptors
import json
import csv
import pkg_resources
from typing import List
import pickle
import copy

from rdkit.Chem.Draw.IPythonConsole import _GetSubstructMatches
from subway.graph_frag_chem.graph_frag_chem import ProductMaker
from subway.graph_frag_chem.IBM_RXN.graph_ibm_rxn import IBM_RXN

JSON_PATH = pkg_resources.resource_filename(
    "subway", "data/reaction_template/reaction_template_nodes.json"
)

ADJ_PATH = pkg_resources.resource_filename(
    "subway", "data/reaction_template/reaction_template_adj.pkl"
)

IBM_ERROR_CSV = pkg_resources.resource_filename(
    "subway", "data/reaction_template/IBM_RXN/ibm_error.csv"
)

# NEEDED FOR IBM_RXN. Generate this from your own account
API_KEY = "apk-ccbb3f6f74af6962119bc9f7461e784a263b14fef8080c320610c016b69057164d7317ae96eff2deda922d13b59ccc330e924ab679d8c7bdada2610ab690867df7277c559da08fc3ee10d269016ee747"
PROJECT_ID = "60e8605799348f0001561731"


class ReactionTemplates:
    """ 
    Class containing reaction templates that will create reaction nodes and product/intermediate nodes from molecule nodes
    """

    def __init__(self, json_path, adjacency_path, ibm_rxn_bool):
        self.json_path = json_path
        self.adj_path = adjacency_path
        file = open(json_path, "rb")
        self.json_data = json.load(file)
        file.close()
        file = open(adjacency_path, "rb")
        self.adj_list = pickle.load(file)
        file.close()
        self.ibm_rxn_bool = ibm_rxn_bool

    def write_json(self, new_data, filename):
        """Function that writes new data to json file (list)
        """
        with open(filename, "r+") as file:
            # First we load existing data into a dict.
            file_data = json.load(file)
            # Join new_data with file_data
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

    def update_error_csv(self, reaction_smile, ibm_error_path, reaction_type):
        """ Function that adds error reaction smile to csv because IBM RXN produces wrong product
        """
        with open(ibm_error_path, "a", newline="") as file:
            writer = csv.writer(file)
            writer.writerow([reaction_type, reaction_smile])

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
        mol_smiles,
        cost_per_mol: float,
        cad: float,
        eq_per_site: float,
        name: str,
        manual: str,
    ):
        """ Function that will create molecule node and update json file.

        Parameters
        ----------
        mol_smiles: SMILES of the molecule
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
        molecule_node["SMILES"] = mol_smiles

        mol_obj = Chem.MolFromSmiles(mol_smiles)

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

    def modify_reaction_new_product(self, old_reaction, new_product):
        """Function that modifies the old reaction by replacing the old product with new product using string manipulation.

        Return
        --------
        new_reaction: Reaction SMILE of new reaction
        """
        counter = 0
        product_string = ""
        for char in old_reaction:
            if counter == 2:
                product_string += char
            if char == ">":
                counter += 1

        new_reaction = old_reaction.replace(product_string, new_product)
        return new_reaction

    def graph_frag_Process(self, reactants_smiles, reagent_smiles, reaction_type):
        """Function that processes the incoming reactions and prepares for molecule, reaction node maker
        
        Parameters
        ----------
        reactants_smiles: list of reactants
        reagent_smiles: reagents for corresponding reaction type
        reaction_type: reaction type as input for graph_frag_chem

        Return
        -------
        product_list: list of products in rdkit mol
        reaction_list: list of reactions SMILES

        """
        product_list = []
        reaction_list = []
        mol_node_list = self.getNodes(reactants_smiles)
        # creating reaction node
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

        product_mol_list = ProductMaker().rct_to_pdt(mol_node_list, reaction_type)
        # replicate reaction for as many products
        for product in product_mol_list:
            product_smile = Chem.MolToSmiles(product)
            rxn = copy.copy(reaction_smile)
            rxn += product_smile
            product_list.append(product_smile)
            reaction_list.append(rxn)

        return product_list, reaction_list

    def IBM_RXN_Process(self, reactants_smiles, reagent_smiles, ibm_reagents):
        """Function that processes the incoming reactions and prepares for molecule, reaction node maker
        Parameters
        ----------
        reactants_smiles: list of reactant pairs
        reagent_smiles: reagents for corresponding reaction type
        ibm_reagents: list of reagents that count as reactants for ibm (with fragment notation "~")

        Return
        -------
        product_list: list of products in rdkit mol
        reaction_list: list of reactions SMILES

        """
        # need to carefully plan out how I want each reaction template to look like, and
        # how this function will reduce the clutter of code but also be integrated seamlessly and intuitively
        reaction_list = []
        # add ibm_reagents to reactants_smiles
        ibm_modified_reagents = []
        for rct_pair in reactants_smiles:
            for ibm_rgt in ibm_reagents:
                rct_pair += "."
                rct_pair += ibm_rgt
            ibm_modified_reagents.append(rct_pair)
        # get products from ibm_rxn
        product_list, ibm_reaction_list = IBM_RXN(
            API_KEY, PROJECT_ID
        ).get_product_batch(ibm_modified_reagents)
        index = 0
        # add reagents, syntax, and product to reactants to make Reaction SMILES
        for index in range(len(reactants_smiles)):
            rxn = reactants_smiles[index]
            rxn += ">"
            rgt_index = 0
            while rgt_index < len(reagent_smiles):
                rxn += reagent_smiles[rgt_index]
                # don't add "." to last reagent
                if rgt_index != len(reagent_smiles) - 1:
                    rxn += "."
                rgt_index += 1
            rxn += ">"
            if product_list[index] == None:
                reaction_list.append(ibm_reaction_list[index])
                continue
            else:
                product_smile = product_list[index]
            rxn += product_smile
            reaction_list.append(rxn)
            index += 1

        return product_list, reaction_list

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
        ibm_reagents = {
            "substruct": "B1OC(=O)CN(C)CC(=O)O1",
            "reagents": [
                "CC(C)C1=CC(=C(C(=C1)C(C)C)C2=CC(=CC=C2)P(C3CCCCC3)C4CCCCC4)C(C)C~C1=CC=C([C-]=C1)C2=CC=CC=C2N~Cl[Pd+]",
                "[O-]P(=O)([O-])[O-]~[K+]~[K+]~[K+]",
            ],
        }
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

        # create product_list and reaction_list from IBM_RXN or graph_frag_chem
        if self.ibm_rxn_bool:
            product_list, reaction_list = self.IBM_RXN_Process(
                reactants_smiles, reagent_smiles, ibm_reagents["reagents"]
            )
        else:
            product_list, reaction_list = self.graph_frag_Process(
                reactants_smiles, reagent_smiles, "wingSuzuki"
            )

        # multiple products, for loop is required
        for i in range(len(product_list)):
            product_smile = product_list[i]
            if product_smile == None:
                self.update_error_csv(reaction_list[i], IBM_ERROR_CSV, "wingSuzuki")
                continue
            else:
                # add product molecule node and returns id of product node
                # check if node exists in json file
                product_mol = Chem.MolFromSmiles(product_list[i])
                reaction_smile = reaction_list[i]
                # # if bmida is not present, then re-run with additional reagent. Can use single prediction
                # bmida_smiles = structures_for_ibm["substruct"]
                # reactant_list = []
                # if not product_mol.HasSubstructMatch(Chem.MolFromSmiles(bmida_smiles)):
                #     rxn = rdChemReactions.ReactionFromSmarts(reaction_smile)
                #     rcts = rxn.GetReactants()
                #     for rct_mol in rcts:
                #         reactant_list.append(Chem.MolToSmiles(rct_mol))
                #     reactant_list.append(structures_for_ibm["reagent"])
                #     product_smile = IBM_RXN(API_KEY, PROJECT_ID).get_product(
                #         reactant_list
                #     )[0]
                #     # modify the reaction with new product
                #     reaction_smile = self.modify_reaction_new_product(
                #         reaction_smile, product_smile
                #     )

                duplicate_node = False
                for node in self.json_data:
                    if node["type"] == "molecule":
                        if node["SMILES"] == product_smile:
                            product_id = node["id"]
                            duplicate_node = True
                if not duplicate_node:
                    product_id = self.molecule_node_maker(
                        product_smile, 0.0, 0.0, 0.0, "-", ""
                    )
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
        ibm_reagents = {
            "substruct": "B1OC(=O)CN(C)CC(=O)O1",
            "reagents": [
                "CC(C)C1=CC(=C(C(=C1)C(C)C)C2=CC(=CC=C2)P(C3CCCCC3)C4CCCCC4)C(C)C~C1=CC=C([C-]=C1)C2=CC=CC=C2N~Cl[Pd+]",
                "[O-]P(=O)([O-])[O-]~[K+]~[K+]~[K+]",
            ],
        }
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

        # create product_list and reaction_list from IBM_RXN or graph_frag_chem
        if self.ibm_rxn_bool:
            product_list, reaction_list = self.IBM_RXN_Process(
                reactants_smiles, reagent_smiles, ibm_reagents["reagents"]
            )
        else:
            product_list, reaction_list = self.graph_frag_Process(
                reactants_smiles, reagent_smiles, "pentamerSuzuki"
            )

        # multiple products, for loop is required
        for i in range(len(product_list)):
            product_smile = product_list[i]
            if product_smile == None:
                self.update_error_csv(reaction_list[i], IBM_ERROR_CSV, "pentamerSuzuki")
                continue
            else:
                # add product molecule node and returns id of product node
                # check if node exists in json file
                product_mol = Chem.MolFromSmiles(product_list[i])
                reaction_smile = reaction_list[i]
                # NOTE: check if single-substituted. If, yes. Then re-run with same B block and product.
                # check for both substructures
                reactant_list = []
                substructure_list = patts[1]["substructure"]
                if product_mol.HasSubstructMatch(
                    substructure_list[0]
                ) or product_mol.HasSubstructMatch(substructure_list[1]):
                    rxn = rdChemReactions.ReactionFromSmarts(reaction_smile)
                    rcts = rxn.GetReactants()
                    for rct_mol in rcts:
                        if rct_mol.HasSubstructMatch(
                            Chem.MolFromSmiles(ibm_reagents["substruct"])
                        ):
                            reactant_list.append(Chem.MolToSmiles(rct_mol))
                    reactant_list.append(product_smile)
                    product_smile = IBM_RXN(API_KEY, PROJECT_ID).get_product(
                        reactant_list
                    )[0]
                    # modify the reaction with new product
                    reaction_smile = self.modify_reaction_new_product(
                        reaction_smile, product_smile
                    )

                duplicate_node = False
                for node in self.json_data:
                    if node["type"] == "molecule":
                        if node["SMILES"] == product_smile:
                            product_id = node["id"]
                            duplicate_node = True
                if not duplicate_node:
                    product_id = self.molecule_node_maker(
                        product_smile, 0.0, 0.0, 0.0, "-", ""
                    )
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
        ibm_reagents = {"reagents": ["C(=O)([O-])[O-]~[K+]~[K+]"]}
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

        # create product_list and reaction_list from IBM_RXN or graph_frag_chem
        if self.ibm_rxn_bool:
            product_list, reaction_list = self.IBM_RXN_Process(
                reactants_smiles, reagent_smiles, ibm_reagents["reagents"]
            )
        else:
            product_list, reaction_list = self.graph_frag_Process(
                reactants_smiles, reagent_smiles, "deBoc"
            )

        # multiple products, for loop is required
        for i in range(len(product_list)):
            product_smile = product_list[i]
            if product_smile == None:
                self.update_error_csv(reaction_list[i], IBM_ERROR_CSV, "deBoc")
                continue
            else:
                # add product molecule node and returns id of product node
                # check if node exists in json file
                product_mol = Chem.MolFromSmiles(product_list[i])
                duplicate_node = False
                for node in self.json_data:
                    if node["type"] == "molecule":
                        if node["SMILES"] == product_smile:
                            product_id = node["id"]
                            duplicate_node = True
                if not duplicate_node:
                    product_id = self.molecule_node_maker(
                        product_smile, 0.0, 0.0, 0.0, "-", ""
                    )
                reaction_smile = reaction_list[i]
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
        ibm_reagents = {
            "reagents": [
                "C1=CC=C(C=C1)C=CC(=O)C=CC2=CC=CC=C2~C1=CC=C(C=C1)C=CC(=O)C=CC2=CC=CC=C2~C1=CC=C(C=C1)C=CC(=O)C=CC2=CC=CC=C2~[Pd]~[Pd]",
                "CN(C)C1=CC=CC=C1C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4",
                "CC(C)(C)[O-]~[Na+]",
            ]
        }
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

        # create product_list and reaction_list from IBM_RXN or graph_frag_chem
        if self.ibm_rxn_bool:
            product_list, reaction_list = self.IBM_RXN_Process(
                reactants_smiles, reagent_smiles, ibm_reagents["reagents"]
            )
        else:
            product_list, reaction_list = self.graph_frag_Process(
                reactants_smiles, reagent_smiles, "BHA"
            )

        # multiple products, for loop is required
        for i in range(len(product_list)):
            product_smile = product_list[i]
            if product_smile == None:
                self.update_error_csv(reaction_list[i], IBM_ERROR_CSV, "BHA")
                continue
            else:
                # add product molecule node and returns id of product node
                # check if node exists in json file
                product_mol = Chem.MolFromSmiles(product_list[i])
                duplicate_node = False
                for node in self.json_data:
                    if node["type"] == "molecule":
                        if node["SMILES"] == product_smile:
                            product_id = node["id"]
                            duplicate_node = True
                if not duplicate_node:
                    product_id = self.molecule_node_maker(
                        product_smile, 0.0, 0.0, 0.0, "-", ""
                    )
                reaction_smile = reaction_list[i]
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
        ibm_reagents = {"reagents": ["C(=O)([O-])[O-]~[Cs+]~[Cs+]"]}
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

        # create product_list and reaction_list from IBM_RXN or graph_frag_chem
        if self.ibm_rxn_bool:
            product_list, reaction_list = self.IBM_RXN_Process(
                reactants_smiles, reagent_smiles, ibm_reagents["reagents"]
            )
        else:
            product_list, reaction_list = self.graph_frag_Process(
                reactants_smiles, reagent_smiles, "SNAr"
            )

        # multiple products, for loop is required
        for i in range(len(product_list)):
            product_smile = product_list[i]
            if product_smile == None:
                self.update_error_csv(reaction_list[i], IBM_ERROR_CSV, "SNAr")
                continue
            else:
                # add product molecule node and returns id of product node
                # check if node exists in json file
                product_mol = Chem.MolFromSmiles(product_list[i])
                duplicate_node = False
                for node in self.json_data:
                    if node["type"] == "molecule":
                        if node["SMILES"] == product_smile:
                            product_id = node["id"]
                            duplicate_node = True
                if not duplicate_node:
                    product_id = self.molecule_node_maker(
                        product_smile, 0.0, 0.0, 0.0, "-", ""
                    )
                reaction_smile = reaction_list[i]
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
                        "SNAr", reaction_smile, t_H, t_M, n_parr, yld, scale
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
