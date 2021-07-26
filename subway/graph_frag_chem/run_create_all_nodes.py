import json
import pickle
import pkg_resources
from rdkit import Chem

from subway.graph_frag_chem.reaction_template import ReactionTemplates
from subway.auto.graph_search_auto import Calculate

JSON_PATH = pkg_resources.resource_filename(
    "subway", "data/reaction_template/reaction_template_nodes.json"
)

IBM_JSON_PATH = pkg_resources.resource_filename(
    "subway", "data/reaction_template/IBM_RXN/ibm_reaction_template_nodes.json"
)

ADJ_PATH = pkg_resources.resource_filename(
    "subway", "data/reaction_template/reaction_template_adj.pkl"
)

IBM_ADJ_PATH = pkg_resources.resource_filename(
    "subway", "data/reaction_template/IBM_RXN/ibm_reaction_template_adj.pkl"
)

SNAR_PATH = pkg_resources.resource_filename(
    "subway", "data/reaction_template/snar_reactants.txt"
)

BHA_PATH = pkg_resources.resource_filename(
    "subway", "data/reaction_template/BHA_reactants.txt"
)

DEBOC_PATH = pkg_resources.resource_filename(
    "subway", "data/reaction_template/deboc_reactants.txt"
)

# STEPS
# clear adj_list
# clear json up to id=63
# clear ibm_error.csv
# add manual
# run all


# {'ibm_mols': 3389, 'ibm_rxns': 3374, 'ibm_ws': 110, 'ibm_ps': 1833,
# 'ibm_deBoc': 204, 'ibm_BHA': 290, 'ibm_SNAr': 932, 'ibm_total': 6763}


class GraphPreBuild:
    """Class containing all the functions necessary to recreate all the nodes and adjacency list.
    """

    def __init__(self, json_path, adjacency_path, ibm_rxn_bool):
        self.json_path = json_path
        self.adjacency_path = adjacency_path

        file = open(json_path, "rb")
        self.json_data = json.load(file)
        file.close()
        file = open(adjacency_path, "rb")
        self.adj_list = pickle.load(file)
        file.close()
        self.ibm_rxn_bool = ibm_rxn_bool  # gives user an option to use graph_frag_chem (FALSE) or IBM_RXN (TRUE)

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
        wS_reactants = []
        for first_node in self.json_data:
            if first_node["type"] == "molecule":
                if first_node["Block_type"] == "a":
                    a_node = first_node
                    for second_node in self.json_data:
                        if second_node["type"] == "molecule":
                            if second_node["Block_type"] == "b":
                                b_node = second_node
                                if self.ibm_rxn_bool:
                                    wS_reactants.append(
                                        a_node["SMILES"] + "." + b_node["SMILES"]
                                    )
                                else:
                                    ReactionTemplates(
                                        self.json_path,
                                        self.adjacency_path,
                                        self.ibm_rxn_bool,
                                    ).wingSuzuki([a_node["SMILES"], b_node["SMILES"]])

        # if IBM_RXN is used, send the entire batch of reactions in tens
        if self.ibm_rxn_bool:
            batch_size = 10
            # separate list of reactants into batches
            for i in range(0, len(wS_reactants), batch_size):
                batch = wS_reactants[i : i + batch_size]
                ReactionTemplates(
                    self.json_path, self.adjacency_path, self.ibm_rxn_bool,
                ).wingSuzuki(batch)
        # updates json data
        self.reload_json()

        # pentamerSuzuki
        pS_reactants = []
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
                                if self.ibm_rxn_bool:
                                    pS_reactants.append(
                                        ab_node["SMILES"] + "." + c_node["SMILES"]
                                    )
                                else:
                                    ReactionTemplates(
                                        self.json_path,
                                        self.adjacency_path,
                                        self.ibm_rxn_bool,
                                    ).pentamerSuzuki(
                                        [ab_node["SMILES"], c_node["SMILES"]]
                                    )
        # # if IBM_RXN is used, send the entire batch of reactions in tens
        # if self.ibm_rxn_bool:
        #     batch_size = 10
        #     # separate list of reactants into batches
        #     for i in range(0, len(pS_reactants), batch_size):
        #         batch = pS_reactants[i : i + batch_size]
        #         ReactionTemplates(
        #             self.json_path, self.adjacency_path, self.ibm_rxn_bool,
        #         ).pentamerSuzuki(batch)
        # # updates json data
        # self.reload_json()

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
        deBoc_reactants = []
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
                        if self.ibm_rxn_bool:
                            deBoc_reactants.append(NBoc_node["SMILES"])
                        else:
                            ReactionTemplates(
                                self.json_path, self.adjacency_path, self.ibm_rxn_bool
                            ).deBoc([NBoc_node["SMILES"]])

        with open(DEBOC_PATH, "w") as f:
            for item in deBoc_reactants:
                f.write("%s\n" % item)

        # # used to skip the searching of reactant pairs (takes a long time)
        # deBoc_reactants = []
        # f = open(SNAR_PATH, "r")
        # for rcts in f:
        #     deBoc_reactants.append(rcts)

        # if IBM_RXN is used, send the entire batch of reactions in tens
        if self.ibm_rxn_bool:
            batch_size = 10
            # separate list of reactants into batches
            for i in range(0, len(deBoc_reactants), batch_size):
                batch = deBoc_reactants[i : i + batch_size]
                ReactionTemplates(
                    self.json_path, self.adjacency_path, self.ibm_rxn_bool,
                ).deBoc(batch)

        # updates json data
        self.reload_json()

        # BHA
        bha_reactants = []
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
                                if self.ibm_rxn_bool:
                                    bha_reactants.append(
                                        nh_node["SMILES"] + "." + halide_node["SMILES"]
                                    )
                                else:
                                    ReactionTemplates(
                                        self.json_path,
                                        self.adjacency_path,
                                        self.ibm_rxn_bool,
                                    ).BHA([nh_node["SMILES"], halide_node["SMILES"]])

        with open(BHA_PATH, "w") as f:
            for item in bha_reactants:
                f.write("%s\n" % item)

        # # used to skip the searching of reactant pairs (takes a long time)
        # bha_reactants = []
        # f = open(SNAR_PATH, "r")
        # for rcts in f:
        #     bha_reactants.append(rcts)

        # if IBM_RXN is used, send the entire batch of reactions in tens
        if self.ibm_rxn_bool:
            batch_size = 10
            # separate list of reactants into batches
            for i in range(0, len(bha_reactants), batch_size):
                batch = bha_reactants[i : i + batch_size]
                ReactionTemplates(
                    self.json_path, self.adjacency_path, self.ibm_rxn_bool,
                ).BHA(batch)
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
        snar_reactants = []
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
                                    if self.ibm_rxn_bool:
                                        snar_reactants.append(
                                            ar_node["SMILES"] + "." + nu_node["SMILES"]
                                        )
                                    else:
                                        ReactionTemplates(
                                            self.json_path,
                                            self.adjacency_path,
                                            self.ibm_rxn_bool,
                                        ).SNAr([ar_node["SMILES"], nu_node["SMILES"]])
        # if IBM_RXN is used, send the entire batch of reactions in tens
        if self.ibm_rxn_bool:
            batch_size = 10
            # separate list of reactants into batches
            for i in range(0, len(snar_reactants), batch_size):
                batch = snar_reactants[i : i + batch_size]
                ReactionTemplates(
                    self.json_path, self.adjacency_path, self.ibm_rxn_bool,
                ).SNAr(batch)

        # updates json data
        self.reload_json()

        # deBoc
        deBoc_reactants = []
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
                        if self.ibm_rxn_bool:
                            deBoc_reactants.append(NBoc_node["SMILES"])
                        else:
                            ReactionTemplates(
                                self.json_path, self.adjacency_path, self.ibm_rxn_bool
                            ).deBoc([NBoc_node["SMILES"]])

        # if IBM_RXN is used, send the entire batch of reactions in tens
        if self.ibm_rxn_bool:
            batch_size = 10
            # separate list of reactants into batches
            for i in range(0, len(deBoc_reactants), batch_size):
                batch = deBoc_reactants[i : i + batch_size]
                ReactionTemplates(
                    self.json_path, self.adjacency_path, self.ibm_rxn_bool,
                ).deBoc(batch)
        # updates json data
        self.reload_json()

        # BHA
        bha_reactants = []
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
                                if self.ibm_rxn_bool:
                                    bha_reactants.append(
                                        nh_node["SMILES"] + "." + halide_node["SMILES"]
                                    )
                                else:
                                    ReactionTemplates(
                                        self.json_path,
                                        self.adjacency_path,
                                        self.ibm_rxn_bool,
                                    ).BHA([nh_node["SMILES"], halide_node["SMILES"]])
        # if IBM_RXN is used, send the entire batch of reactions in tens
        if self.ibm_rxn_bool:
            batch_size = 10
            # separate list of reactants into batches
            for i in range(0, len(bha_reactants), batch_size):
                batch = bha_reactants[i : i + batch_size]
                ReactionTemplates(
                    self.json_path, self.adjacency_path, self.ibm_rxn_bool,
                ).BHA(batch)
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
        deBoc_reactants = []
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
                        if self.ibm_rxn_bool:
                            deBoc_reactants.append(NBoc_node["SMILES"])
                        else:
                            ReactionTemplates(
                                self.json_path, self.adjacency_path, self.ibm_rxn_bool
                            ).deBoc([NBoc_node["SMILES"]])

        # if IBM_RXN is used, send the entire batch of reactions in tens
        if self.ibm_rxn_bool:
            batch_size = 10
            # separate list of reactants into batches
            for i in range(0, len(deBoc_reactants), batch_size):
                batch = deBoc_reactants[i : i + batch_size]
                ReactionTemplates(
                    self.json_path, self.adjacency_path, self.ibm_rxn_bool,
                ).deBoc(batch)
        # updates json data
        self.reload_json()

        # BHA
        bha_reactants = []
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
                                if self.ibm_rxn_bool:
                                    bha_reactants.append(
                                        nh_node["SMILES"] + "." + halide_node["SMILES"]
                                    )
                                else:
                                    ReactionTemplates(
                                        self.json_path,
                                        self.adjacency_path,
                                        self.ibm_rxn_bool,
                                    ).BHA([nh_node["SMILES"], halide_node["SMILES"]])
        # if IBM_RXN is used, send the entire batch of reactions in tens
        if self.ibm_rxn_bool:
            batch_size = 10
            # separate list of reactants into batches
            for i in range(0, len(bha_reactants), batch_size):
                batch = bha_reactants[i : i + batch_size]
                ReactionTemplates(
                    self.json_path, self.adjacency_path, self.ibm_rxn_bool,
                ).BHA(batch)
        # updates json data
        self.reload_json()

        # SNAr
        snar_reactants = []
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
                                    if self.ibm_rxn_bool:
                                        snar_reactants.append(
                                            ar_node["SMILES"] + "." + nu_node["SMILES"]
                                        )
                                    else:
                                        ReactionTemplates(
                                            self.json_path,
                                            self.adjacency_path,
                                            self.ibm_rxn_bool,
                                        ).SNAr([ar_node["SMILES"], nu_node["SMILES"]])
        # if IBM_RXN is used, send the entire batch of reactions in tens
        if self.ibm_rxn_bool:
            batch_size = 10
            # separate list of reactants into batches
            for i in range(0, len(snar_reactants), batch_size):
                batch = snar_reactants[i : i + batch_size]
                ReactionTemplates(
                    self.json_path, self.adjacency_path, self.ibm_rxn_bool,
                ).SNAr(batch)

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
        snar_reactants = []
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
                                    if self.ibm_rxn_bool:
                                        snar_reactants.append(
                                            ar_node["SMILES"] + "." + nu_node["SMILES"]
                                        )
                                    else:
                                        ReactionTemplates(
                                            self.json_path,
                                            self.adjacency_path,
                                            self.ibm_rxn_bool,
                                        ).SNAr([ar_node["SMILES"], nu_node["SMILES"]])
        with open(SNAR_PATH, "w") as f:
            for item in snar_reactants:
                f.write("%s\n" % item)

        # # used to skip the searching of reactant pairs (takes a long time)
        # snar_reactants = []
        # f = open(SNAR_PATH, "r")
        # for rcts in f:
        #     snar_reactants.append(rcts)

        # if IBM_RXN is used, send the entire batch of reactions in tens
        if self.ibm_rxn_bool:
            batch_size = 10
            # separate list of reactants into batches
            for i in range(0, len(snar_reactants), batch_size):
                batch = snar_reactants[i : i + batch_size]
                ReactionTemplates(
                    self.json_path, self.adjacency_path, self.ibm_rxn_bool,
                ).SNAr(batch)

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
        # self.route_Base()
        # creates all deBoc+BHA intermediates/products (excluding routes that required SNAr products)
        self.route_Buch()
        # creates all SNAr intermediates/products
        self.route_SNAr()
        # creates remainder of routes that require SNAr products
        self.route_Buch()


# change to IBM_JSON_PATH and IBM_ADJ_PATH, if True
prebuild = GraphPreBuild(IBM_JSON_PATH, IBM_ADJ_PATH, True)
# prebuild.clear_pkl(IBM_ADJ_PATH)
# prebuild.read_pkl(IBM_ADJ_PATH)
# prebuild.route_Base()
# prebuild.route_SNAr()
# prebuild.route_Buch()
# prebuild.route_All()

# tester
# mol = Chem.MolFromSmiles(
#     "Fc1cc(F)cc(-c2cc(-c3ccc4oc(-c5cc(-c6cc(F)cc(F)c6)cc(C(F)(F)F)c5)cc4c3)cc(C(F)(F)F)c2)c1"
# )
# file = open(JSON_PATH, "rb")
# json_data = json.load(file)
# file.close()

# ReactionTemplates(JSON_PATH, ADJ_PATH).SNAr(
#     [
#         "Fc1cc(-c2cccc(F)c2-c2ccc3oc(-c4c(F)cccc4-c4ccnc(F)c4)cc3c2)ccn1",
#         "c1ccc2c(c1)[nH]c1ccccc12",
#     ]
# )

