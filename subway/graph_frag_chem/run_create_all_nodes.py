import json
import pickle
import pkg_resources
from rdkit import Chem

from subway.graph_frag_chem.reaction_template import ReactionTemplates
from subway.auto.graph_search_auto import Calculate

JSON_PATH = pkg_resources.resource_filename(
    "subway", "data/reaction_template_nodes.json"
)

ADJ_PATH = pkg_resources.resource_filename("subway", "data/reaction_template_adj.pkl")

# import building blocks
# import reaction_template.py
# figure out a way to combine the right building blocks for input to reaction template
# ex. which reaction template to use?

# clear adjacency list


class GraphPreBuild:
    """Class containing all the functions necessary to recreate all the nodes and adjacency list.
    """

    def __init__(self, json_path, adjacency_path):
        self.json_path = json_path
        self.adjacency_path = adjacency_path

        file = open(json_path, "rb")
        self.json_data = json.load(file)
        file.close()
        file = open(adjacency_path, "rb")
        self.adj_list = pickle.load(file)
        file.close()

    def clear_pkl(self, pkl_path):
        """Function that will clear the pickle file and replace it with an empty list.
        """
        empty = []
        file = open(pkl_path, "wb")
        pickle.dump(empty, file)
        file.close()

    def read_pkl(self, pkl_path):
        """Function that will read out pickle file.
        """
        file = open(pkl_path, "rb")
        data = pickle.load(file)
        file.close()
        print(len(data))

    def reload_json(self):
        """Function that will re-load json data to ensure most updated json file is used.
        """
        file = open(self.json_path, "rb")
        self.json_data = json.load(file)
        file.close()

    def route_Base(self):
        """ Function that will recreate all the intermediates/products and reactions for the Base routes.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        # wingSuzuki
        for first_node in self.json_data:
            if first_node["type"] == "molecule":
                if first_node["Block_type"] == "a":
                    a_node = first_node
                    for second_node in self.json_data:
                        if second_node["type"] == "molecule":
                            if second_node["Block_type"] == "b":
                                b_node = second_node
                                ReactionTemplates(
                                    self.json_path, self.adjacency_path
                                ).wingSuzuki([a_node["SMILES"], b_node["SMILES"]])
        # updates json data
        self.reload_json()

        # pentamerSuzuki
        # match smiles_BMIDA substructure to get intermediate
        for first_node in self.json_data:
            if first_node["type"] == "molecule" and first_node["Block_type"] == "-":
                smiles_type = Calculate().molecule_classifier(
                    first_node, "pentamerSuzuki"
                )[0]
                if smiles_type == "smiles_BMIDA":
                    ab_node = first_node
                    for second_node in self.json_data:
                        if second_node["type"] == "molecule":
                            if second_node["Block_type"] == "c":
                                c_node = second_node
                                ReactionTemplates(
                                    self.json_path, self.adjacency_path
                                ).pentamerSuzuki([ab_node["SMILES"], c_node["SMILES"]])
        # updates json data
        self.reload_json()

    def route_Buch(self):
        """ Function that will recreate all the intermediates/products and reactions for the Buch routes.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        # deBoc
        # match tert-Butyloxycarbonyl substructure
        for first_node in self.json_data:
            if first_node["type"] == "molecule" and first_node["Block_type"] == "-":
                mol = Chem.MolFromSmiles(first_node["SMILES"])
                bmida = Chem.MolFromSmiles("B1OC(=O)CN(C)CC(=O)O1")
                if not mol.HasSubstructMatch(bmida):
                    smiles_type_first = Calculate().molecule_classifier(
                        first_node, "deBoc"
                    )[0]
                    if smiles_type_first == "smiles_NBoc":
                        NBoc_node = first_node
                        ReactionTemplates(self.json_path, self.adjacency_path).deBoc(
                            [NBoc_node["SMILES"]]
                        )
        # updates json data
        self.reload_json()

        # BHA
        # match NH and halide substructure
        for first_node in self.json_data:
            if first_node["type"] == "molecule" and first_node["Block_type"] == "-":
                mol = Chem.MolFromSmiles(first_node["SMILES"])
                bmida = Chem.MolFromSmiles("B1OC(=O)CN(C)CC(=O)O1")
                if not mol.HasSubstructMatch(bmida):
                    smiles_type_first = Calculate().molecule_classifier(
                        first_node, "BHA"
                    )[0]
                    if smiles_type_first == "smiles_NH":
                        nh_node = first_node
                        for second_node in self.json_data:
                            if (
                                second_node["type"] == "molecule"
                                and second_node["Block_type"] == "buch"
                            ):
                                halide_node = second_node
                                ReactionTemplates(
                                    self.json_path, self.adjacency_path
                                ).BHA([nh_node["SMILES"], halide_node["SMILES"]])
        # updates json data
        self.reload_json()

    def route_SB(self):
        """ Function that will recreate all the intermediates/products and reactions for the SB routes.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        # SNAr
        # match aryl fluoride and carbazole substructures
        for first_node in self.json_data:
            if first_node["type"] == "molecule" and first_node["Block_type"] == "-":
                mol = Chem.MolFromSmiles(first_node["SMILES"])
                bmida = Chem.MolFromSmiles("B1OC(=O)CN(C)CC(=O)O1")
                if not mol.HasSubstructMatch(bmida):
                    smiles_type_first = Calculate().molecule_classifier(
                        first_node, "SNAr"
                    )[0]
                    if smiles_type_first == "smiles_ArF":
                        ar_node = first_node
                        for second_node in self.json_data:
                            if second_node["type"] == "molecule":
                                smiles_type_second = Calculate().molecule_classifier(
                                    second_node, "SNAr"
                                )[0]
                                if smiles_type_second == "smiles_Nu":
                                    nu_node = second_node
                                    ReactionTemplates(
                                        self.json_path, self.adjacency_path
                                    ).SNAr([ar_node["SMILES"], nu_node["SMILES"]])
        # updates json data
        self.reload_json()

        # deBoc
        # match tert-Butyloxycarbonyl substructure
        for first_node in self.json_data:
            if first_node["type"] == "molecule" and first_node["Block_type"] == "-":
                mol = Chem.MolFromSmiles(first_node["SMILES"])
                bmida = Chem.MolFromSmiles("B1OC(=O)CN(C)CC(=O)O1")
                if not mol.HasSubstructMatch(bmida):
                    smiles_type_first = Calculate().molecule_classifier(
                        first_node, "deBoc"
                    )[0]
                    if smiles_type_first == "smiles_NBoc":
                        NBoc_node = first_node
                        ReactionTemplates(self.json_path, self.adjacency_path).deBoc(
                            [NBoc_node["SMILES"]]
                        )
        # updates json data
        self.reload_json()

        # BHA
        # match NH and halide substructure
        for first_node in self.json_data:
            if first_node["type"] == "molecule" and first_node["Block_type"] == "-":
                mol = Chem.MolFromSmiles(first_node["SMILES"])
                bmida = Chem.MolFromSmiles("B1OC(=O)CN(C)CC(=O)O1")
                if not mol.HasSubstructMatch(bmida):
                    smiles_type_first = Calculate().molecule_classifier(
                        first_node, "BHA"
                    )[0]
                    if smiles_type_first == "smiles_NH":
                        nh_node = first_node
                        for second_node in self.json_data:
                            if second_node["type"] == "molecule":
                                smiles_type_second = Calculate().molecule_classifier(
                                    second_node, "BHA"
                                )[0]
                                if smiles_type_second == "smiles_halide":
                                    halide_node = second_node
                                    ReactionTemplates(
                                        self.json_path, self.adjacency_path
                                    ).BHA([nh_node["SMILES"], halide_node["SMILES"]])
        # updates json data
        self.reload_json()

    def route_BS(self):
        """ Function that will recreate all the intermediates/products and reactions for the BS routes.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        # deBoc
        # match tert-Butyloxycarbonyl substructure
        for first_node in self.json_data:
            if first_node["type"] == "molecule" and first_node["Block_type"] == "-":
                mol = Chem.MolFromSmiles(first_node["SMILES"])
                bmida = Chem.MolFromSmiles("B1OC(=O)CN(C)CC(=O)O1")
                if not mol.HasSubstructMatch(bmida):
                    smiles_type_first = Calculate().molecule_classifier(
                        first_node, "deBoc"
                    )[0]
                    if smiles_type_first == "smiles_NBoc":
                        NBoc_node = first_node
                        ReactionTemplates(self.json_path, self.adjacency_path).deBoc(
                            [NBoc_node["SMILES"]]
                        )
        # updates json data
        self.reload_json()

        # BHA
        # match NH and halide substructure
        for first_node in self.json_data:
            if first_node["type"] == "molecule" and first_node["Block_type"] == "-":
                mol = Chem.MolFromSmiles(first_node["SMILES"])
                bmida = Chem.MolFromSmiles("B1OC(=O)CN(C)CC(=O)O1")
                if not mol.HasSubstructMatch(bmida):
                    smiles_type_first = Calculate().molecule_classifier(
                        first_node, "BHA"
                    )[0]
                    if smiles_type_first == "smiles_NH":
                        nh_node = first_node
                        for second_node in self.json_data:
                            if second_node["type"] == "molecule":
                                smiles_type_second = Calculate().molecule_classifier(
                                    second_node, "BHA"
                                )[0]
                                if smiles_type_second == "smiles_halide":
                                    halide_node = second_node
                                    ReactionTemplates(
                                        self.json_path, self.adjacency_path
                                    ).BHA([nh_node["SMILES"], halide_node["SMILES"]])
        # updates json data
        self.reload_json()

        # SNAr
        # match aryl fluoride and carbazole substructures
        for first_node in self.json_data:
            if first_node["type"] == "molecule" and first_node["Block_type"] == "-":
                mol = Chem.MolFromSmiles(first_node["SMILES"])
                bmida = Chem.MolFromSmiles("B1OC(=O)CN(C)CC(=O)O1")
                if not mol.HasSubstructMatch(bmida):
                    smiles_type_first = Calculate().molecule_classifier(
                        first_node, "SNAr"
                    )[0]
                    if smiles_type_first == "smiles_ArF":
                        ar_node = first_node
                        for second_node in self.json_data:
                            if second_node["type"] == "molecule":
                                smiles_type_second = Calculate().molecule_classifier(
                                    second_node, "SNAr"
                                )[0]
                                if smiles_type_second == "smiles_Nu":
                                    nu_node = second_node
                                    ReactionTemplates(
                                        self.json_path, self.adjacency_path
                                    ).SNAr([ar_node["SMILES"], nu_node["SMILES"]])
        # updates json data
        self.reload_json()

    def route_SNAr(self):
        """ Function that will recreate all the intermediates/products and reactions for the SNAr routes.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        # SNAr
        # match aryl fluoride and carbazole substructures
        for first_node in self.json_data:
            if first_node["type"] == "molecule" and first_node["Block_type"] == "-":
                # checking if SMILES does not have bmida substructure
                # ensures that only pentamerSuzuki products are used
                mol = Chem.MolFromSmiles(first_node["SMILES"])
                bmida = Chem.MolFromSmiles("B1OC(=O)CN(C)CC(=O)O1")
                if not mol.HasSubstructMatch(bmida):
                    smiles_type_first = Calculate().molecule_classifier(
                        first_node, "SNAr"
                    )[0]
                    if smiles_type_first == "smiles_ArF":
                        ar_node = first_node
                        for second_node in self.json_data:
                            if (
                                second_node["type"] == "molecule"
                                and second_node["Block_type"] == "snar"
                            ):
                                nu_node = second_node
                                ReactionTemplates(
                                    self.json_path, self.adjacency_path
                                ).SNAr([ar_node["SMILES"], nu_node["SMILES"]])
        # updates json data
        self.reload_json()

    def route_All(self):
        """ Function that will recreate all the intermediates/products and reactions for All routes.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        # creates all wingSuzuki and pentamerSuzuki intermediates/products
        self.route_Base()
        # creates all SNAr intermediates/products (excluding routes that required deBoc+BHA products)
        self.route_SNAr()
        # creates all deBoc+BHA intermediates/products
        self.route_Buch()
        # creates remainder of routes that require deBoc+BHA products
        self.route_SNAr()


prebuild = GraphPreBuild(JSON_PATH, ADJ_PATH)
# prebuild.clear_pkl(ADJ_PATH)
# prebuild.read_pkl(ADJ_PATH)
# prebuild.route_Base()
# prebuild.route_SNAr()
# prebuild.route_Buch()
prebuild.route_All()
