import json
import pickle
import pkg_resources
from rdkit import Chem
from rdkit.Chem import rdChemReactions

JSON_PATH = pkg_resources.resource_filename(
    "subway", "data/reaction_template_nodes.json"
)

MAN_JSON_PATH = pkg_resources.resource_filename("subway", "data/manual_nodes.json")

ADJ_PATH = pkg_resources.resource_filename("subway", "data/reaction_template_adj.pkl")

# NOTE: reaction_template_nodes.json must contain A, B, C building blocks and any reagent and all manual molecules

# add manual reaction nodes and manual connections for adjacency list
class Manual_Reactions:
    """Class containing all functions and data necessary to create manual 
    reaction nodes and manual connections in adjacency list.
    """

    def __init__(self, json_path, adjacency_path, man_json_path):

        self.json_path = json_path
        self.adj_path = adjacency_path
        self.man_json_path = man_json_path

        file = open(man_json_path, "rb")
        self.man_json_data = json.load(file)
        file.close()

        file = open(json_path, "rb")
        self.json_data = json.load(file)
        file.close()

        file = open(adjacency_path, "rb")
        self.adj_list = pickle.load(file)
        file.close()

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

    def write_json(self, new_data, filename):
        """Function that writes new data (list) to json file (list)
        """
        with open(filename, "r+") as file:
            # First we load existing data into a dict.
            file_data = json.load(file)
            # Join new_dat3a with file_data
            file_data.extend(new_data)
            # Sets file's current position at offset.
            file.seek(0)
            # convert back to json.
            json.dump(file_data, file, indent=4)

    def add_manual_nodes(self):
        """Function that adds manual nodes to json file.
        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        man_id = self.getNextId()
        for man_node in self.man_json_data:
            man_node["id"] = man_id
            man_id += 1
        self.write_json(self.man_json_data, self.json_path)

    def add_manual_adj_list(self):
        """Function that adds manual connections to adjacency list.
        NOTE: must perform after manual reaction nodes are added.
        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        man_adj_list = []
        for man_node in self.man_json_data:
            rxn = rdChemReactions.ReactionFromSmarts(man_node["rxn_SMILES"])
            # find reaction node
            for node in self.json_data:
                if node["type"] == "reaction":
                    if node["rxn_SMILES"] == man_node["rxn_SMILES"]:
                        rxn_id = node["id"]

            # find reactant nodes
            reactant_list = rxn.GetReactants()
            for rct in reactant_list:
                rct_smile = Chem.MolToSmiles(rct)
                rct_smile = Chem.CanonSmiles(rct_smile)
                for node in self.json_data:
                    if node["type"] == "molecule":
                        if Chem.CanonSmiles(node["SMILES"]) == rct_smile:
                            man_adj_list.append([node["id"], rxn_id])

            # find reagent nodes
            reagent_list = rxn.GetAgents()
            for rgt in reagent_list:
                rgt_smile = Chem.MolToSmiles(rgt)
                rgt_smile = Chem.CanonSmiles(rgt_smile)
                for node in self.json_data:
                    if node["type"] == "molecule":
                        if Chem.CanonSmiles(node["SMILES"]) == rgt_smile:
                            man_adj_list.append([node["id"], rxn_id])

            # find product node
            product_list = rxn.GetProducts()
            for pdt in product_list:
                pdt_smile = Chem.MolToSmiles(pdt)
                pdt_smile = Chem.CanonSmiles(pdt_smile)
                for node in self.json_data:
                    if node["type"] == "molecule":
                        if Chem.CanonSmiles(node["SMILES"]) == pdt_smile:
                            man_adj_list.append([rxn_id, node["id"]])
        print(man_adj_list)

        file = open(self.adj_path, "wb")
        self.adj_list.extend(man_adj_list)
        pickle.dump(self.adj_list, file)
        file.close()


man = Manual_Reactions(JSON_PATH, ADJ_PATH, MAN_JSON_PATH)
man.add_manual_nodes()
man = Manual_Reactions(JSON_PATH, ADJ_PATH, MAN_JSON_PATH)
man.add_manual_adj_list()

# mol1 = Chem.MolFromSmiles("c1ccc(P(c2ccccc2)c2ccccc2)cc1")
# smile1 = Chem.CanonSmiles("c1ccc(P(c2ccccc2)c2ccccc2)cc1")
# mol2 = Chem.MolFromSmiles("P(c1ccccc1)(c2ccccc2)c3ccccc3")
# smile2 = Chem.CanonSmiles("P(c1ccccc1)(c2ccccc2)c3ccccc3")
# print(smile1, smile2)
# print(mol1 == mol2)
