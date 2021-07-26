import json
from subway.graph_frag_chem.run_create_all_nodes import IBM_JSON_PATH
import pandas as pd
import pkg_resources
import rdkit
from rdkit import Chem
from rdkit.Chem import rdmolops

ROUTESCORE_JSON = pkg_resources.resource_filename("subway", "data/auto_nodes.json")
FRAGCHEM_JSON = pkg_resources.resource_filename(
    "subway", "data/reaction_template/reaction_template_nodes.json"
)
IBM_JSON_PATH = pkg_resources.resource_filename(
    "subway", "data/reaction_template/IBM_RXN/ibm_reaction_template_nodes_new.json"
)
COMPARE_CSV = pkg_resources.resource_filename(
    "subway", "data/reaction_template/json_comparison.csv"
)
NOT_MATCH_CSV = pkg_resources.resource_filename(
    "subway", "data/reaction_template/not_matched_debug.csv"
)
FULL_PROPS_CSV = pkg_resources.resource_filename(
    "subway", "data/subway_maps/full_props.csv"
)


class Compare:
    """ Class that contains functions for comparing the routescore data to the new graph data
    """

    def __init__(
        self, routescore_json, fragchem_json, ibm_rxn_json, compare_csv, not_match_csv
    ):
        file = open(routescore_json, "rb")
        self.routescore_data = json.load(file)
        file.close()
        file = open(fragchem_json, "rb")
        self.fragchem_data = json.load(file)
        file.close()
        self.compare_path = compare_csv
        self.not_match_path = not_match_csv
        file = open(ibm_rxn_json, "rb")
        self.ibm_rxn_data = json.load(file)
        file.close()

    def check_duplicates(self):
        """Function that checks for duplicates in the fragchem json file"""
        for node in self.fragchem_data:
            if node["type"] == "molecule":
                current_id = node["id"]
                current_smile = node["SMILES"]
                for alt_node in self.fragchem_data:
                    if alt_node["type"] == "molecule":
                        if (
                            alt_node["id"] != current_id
                            and alt_node["SMILES"] == current_smile
                        ):
                            print(alt_node["id"], current_id)

    def full_props_comparison(self, full_props_path):
        """Function that compares number of matches in fragchem json file to full_props.csv"""
        matches = 0
        full_props_df = pd.read_csv(full_props_path)
        for product in full_props_df["smiles"]:
            for node in self.fragchem_data:
                if node["type"] == "molecule":
                    if product == node["SMILES"]:
                        matches += 1
        print(matches)

    def transforms_bmida_check(self, rs_smiles, ibm_smiles):
        """Function that converts the bmida structure in both SMILES to be the same, 
        and then returns whether they are the same (True) or not (False)
        """
        original_structure = "C[N+]12CC(=O)O[B-]1OC(=O)C2"
        bmida_modified_structure = "B1OC(=O)CN(C)CC(=O)O1"
        original_mol = Chem.MolFromSmiles(original_structure)
        bmida_mod_mol = Chem.MolFromSmiles(bmida_modified_structure)

        rs_mol = Chem.MolFromSmiles(rs_smiles)
        ibm_mol = Chem.MolFromSmiles(ibm_smiles)

        if rs_mol.HasSubstructMatch(bmida_mod_mol):
            rs_mol = rdmolops.ReplaceSubstructs(
                rs_mol, bmida_mod_mol, bmida_mod_mol, replaceAll=True
            )[0]

        rs_smiles = Chem.MolToSmiles(rs_mol)

        return rs_smiles == ibm_smiles

    def get_product_from_rxn(self, reaction_smile: str):
        """Function that does string concatenation to find product in reaction_smile
        """
        counter = 0
        product_smile = ""
        for char in reaction_smile:
            if counter == 2:
                product_smile += char
            if char == ">":
                counter += 1

        return product_smile

    def smiles_match(self):
        """ Function that counts the number of smiles matches in the routescore data to fragchem data

        Parameters
        ----------
        None

        Return (in .csv format)
        -------
        mol_matches: number of molecule matches
        rxn_matches: number of reaction matches
        total_matches: total number of matches
        rs_mols: total number of molecule nodes in routescore data
        rs_rxn: total number of reaction nodes in routescore data
        rs_xxx: total number of reaction nodes with specific type of reaction is routescore data
        rs_total: total number of nodes in routescore data
        fg_mols: total number of molecule nodes in fragchem data
        fg_rxn: total number of reaction nodes in fragchem data
        fg_xxx: total number of reaction nodes with specific type of reaction is fragchem data
        fg_total: total number of nodes in fragchem data
        """
        mol_matches = 0
        rxn_matches = 0
        total_matches = 0
        rs_mols = 0
        rs_rxns = 0
        rs_ws = 0
        rs_ps = 0
        rs_deBoc = 0
        rs_BHA = 0
        rs_SNAr = 0
        rs_total = 0
        fg_mols = 0
        fg_rxns = 0
        fg_ws = 0
        fg_ps = 0
        fg_deBoc = 0
        fg_BHA = 0
        fg_SNAr = 0
        fg_total = 0
        # finds number of matches between two datasets
        # find list of smiles that are not matched in routescore data
        matched_list = []
        for rs_node in self.routescore_data:
            for fg_node in self.fragchem_data:
                if rs_node["type"] == "molecule" and fg_node["type"] == "molecule":
                    if rs_node["SMILES"] == fg_node["SMILES"]:
                        mol_matches += 1
                        matched_list.append(rs_node["SMILES"])
                elif rs_node["type"] == "reaction" and fg_node["type"] == "reaction":
                    if rs_node["rxn_SMILES"] == fg_node["rxn_SMILES"]:
                        rxn_matches += 1
        total_matches = mol_matches + rxn_matches

        # use matched_list to find the not matched smiles in routescore data
        # list of routescore smiles
        rs_mol_smiles_list = []
        for rs_node in self.routescore_data:
            if rs_node["type"] == "molecule":
                rs_mol_smiles_list.append(rs_node["SMILES"])

        matched_list = set(matched_list)
        rs_mol_smiles_list = set(rs_mol_smiles_list)
        not_matched_list = rs_mol_smiles_list - matched_list
        not_matched_dict = {}
        not_matched_dict["not_matched"] = []
        for not_match in not_matched_list:
            not_matched_dict["not_matched"].append(not_match)

        not_match_df = pd.DataFrame(data=not_matched_dict)
        not_match_df.to_csv(self.not_match_path)

        # finds number of molecules and reactions in routescore data
        for rs_node in self.routescore_data:
            if rs_node["type"] == "molecule":
                rs_mols += 1
            elif rs_node["type"] == "reaction":
                if rs_node["reaction_type"] == "wingSuzuki":
                    rs_ws += 1
                elif rs_node["reaction_type"] == "pentamerSuzuki":
                    rs_ps += 1
                elif rs_node["reaction_type"] == "deBoc":
                    rs_deBoc += 1
                elif rs_node["reaction_type"] == "BHA":
                    rs_BHA += 1
                elif rs_node["reaction_type"] == "SNAr":
                    rs_SNAr += 1
                rs_rxns += 1
            rs_total += 1

        # finds number of molecules and reactions in fragchem data
        for fg_node in self.fragchem_data:
            if fg_node["type"] == "molecule":
                fg_mols += 1
            elif fg_node["type"] == "reaction":
                if fg_node["reaction_type"] == "wingSuzuki":
                    fg_ws += 1
                elif fg_node["reaction_type"] == "pentamerSuzuki":
                    fg_ps += 1
                elif fg_node["reaction_type"] == "deBoc":
                    fg_deBoc += 1
                elif fg_node["reaction_type"] == "BHA":
                    fg_BHA += 1
                elif fg_node["reaction_type"] == "SNAr":
                    fg_SNAr += 1
                fg_rxns += 1
            fg_total += 1

        matches_dict = {
            "mol_matches": mol_matches,
            "rxn_matches": rxn_matches,
            "total_matches": total_matches,
            "rs_mols": rs_mols,
            "rs_rxns": [rs_rxns, rs_ws, rs_ps, rs_deBoc, rs_BHA, rs_SNAr],
            "rs_total": rs_total,
            "fg_mols": fg_mols,
            "fg_rxns": [fg_rxns, fg_ws, fg_ps, fg_deBoc, fg_BHA, fg_SNAr],
            "fg_total": fg_total,
        }

        df = pd.DataFrame(
            data=matches_dict,
            index=["number", "wingSuzuki", "pentamerSuzuki", "deBoc", "BHA", "SNAr"],
        )
        df.to_csv(self.compare_path)

    def mol_match(self):
        """ Function that counts the number of smiles matches in the routescore data to IBM_RXN data

        Parameters
        ----------
        None

        Return (in .csv format)
        -------
        mol_matches: number of molecule matches
        rxn_matches: number of reaction matches
        total_matches: total number of matches
        rs_mols: total number of molecule nodes in routescore data
        rs_rxn: total number of reaction nodes in routescore data
        rs_xxx: total number of reaction nodes with specific type of reaction is routescore data
        rs_total: total number of nodes in routescore data
        ibm_mols: total number of molecule nodes in IBM_RXN data
        ibm_rxn: total number of reaction nodes in IBM_RXN data
        ibm_xxx: total number of reaction nodes with specific type of reaction is IBM_RXN data
        ibm_total: total number of nodes in IBM_RXN data
        """
        mol_matches = 0
        rxn_matches = 0
        total_matches = 0
        rs_mols = 0
        rs_rxns = 0
        rs_ws = 0
        rs_ps = 0
        rs_deBoc = 0
        rs_BHA = 0
        rs_SNAr = 0
        rs_total = 0
        ibm_mols = 0
        ibm_rxns = 0
        ibm_ws = 0
        ibm_ps = 0
        ibm_deBoc = 0
        ibm_BHA = 0
        ibm_SNAr = 0
        ibm_total = 0
        # finds number of matches between two datasets
        # find list of mols that are not matched in routescore data
        matched_list = []
        print("START")
        index = 0
        for rs_node in self.routescore_data:
            print("count: ", index)
            for ibm_node in self.ibm_rxn_data:
                if rs_node["type"] == "molecule" and ibm_node["type"] == "molecule":
                    if self.transforms_bmida_check(
                        rs_node["SMILES"], ibm_node["SMILES"]
                    ):
                        mol_matches += 1
                        matched_list.append(rs_node["SMILES"])
            index += 1
            # elif rs_node["type"] == "reaction" and ibm_node["type"] == "reaction":
            #     if rs_node["rxn_SMILES"] == ibm_node["rxn_SMILES"]:
            #         rxn_matches += 1
        print("DONE")
        total_matches = mol_matches + rxn_matches

        # use matched_list to find the not matched smiles in routescore data
        # list of routescore smiles
        rs_mol_smiles_list = []
        for rs_node in self.routescore_data:
            if rs_node["type"] == "molecule":
                rs_mol_smiles_list.append(rs_node["SMILES"])

        matched_list = set(matched_list)
        rs_mol_smiles_list = set(rs_mol_smiles_list)
        not_matched_list = rs_mol_smiles_list - matched_list
        not_matched_dict = {}
        not_matched_dict["not_matched"] = []
        for not_match in not_matched_list:
            not_matched_dict["not_matched"].append(not_match)

        not_match_df = pd.DataFrame(data=not_matched_dict)
        not_match_df.to_csv(self.not_match_path)

        # finds number of molecules and reactions in routescore data
        for rs_node in self.routescore_data:
            if rs_node["type"] == "molecule":
                rs_mols += 1
            elif rs_node["type"] == "reaction":
                if rs_node["reaction_type"] == "wingSuzuki":
                    rs_ws += 1
                elif rs_node["reaction_type"] == "pentamerSuzuki":
                    rs_ps += 1
                elif rs_node["reaction_type"] == "deBoc":
                    rs_deBoc += 1
                elif rs_node["reaction_type"] == "BHA":
                    rs_BHA += 1
                elif rs_node["reaction_type"] == "SNAr":
                    rs_SNAr += 1
                rs_rxns += 1
            rs_total += 1

        # finds number of molecules and reactions in fragchem data
        for ibm_node in self.ibm_rxn_data:
            if ibm_node["type"] == "molecule":
                ibm_mols += 1
            elif ibm_node["type"] == "reaction":
                if ibm_node["reaction_type"] == "wingSuzuki":
                    ibm_ws += 1
                elif ibm_node["reaction_type"] == "pentamerSuzuki":
                    ibm_ps += 1
                elif ibm_node["reaction_type"] == "deBoc":
                    ibm_deBoc += 1
                elif ibm_node["reaction_type"] == "BHA":
                    ibm_BHA += 1
                elif ibm_node["reaction_type"] == "SNAr":
                    ibm_SNAr += 1
                ibm_rxns += 1
            ibm_total += 1

        matches_dict = {
            "mol_matches": mol_matches,
            "rxn_matches": rxn_matches,
            "total_matches": total_matches,
            "rs_mols": rs_mols,
            "rs_rxns": [rs_rxns, rs_ws, rs_ps, rs_deBoc, rs_BHA, rs_SNAr],
            "rs_total": rs_total,
            "ibm_mols": ibm_mols,
            "ibm_rxns": [ibm_rxns, ibm_ws, ibm_ps, ibm_deBoc, ibm_BHA, ibm_SNAr],
            "ibm_total": ibm_total,
        }

        df = pd.DataFrame(
            data=matches_dict,
            index=["number", "wingSuzuki", "pentamerSuzuki", "deBoc", "BHA", "SNAr"],
        )
        df.to_csv(self.compare_path)

    def count_mol_rxn(self):
        ibm_mols = 0
        ibm_rxns = 0
        ibm_ws = 0
        ibm_ps = 0
        ibm_deBoc = 0
        ibm_BHA = 0
        ibm_SNAr = 0
        ibm_total = 0
        # finds number of molecules and reactions in fragchem data
        for ibm_node in self.ibm_rxn_data:
            if ibm_node["type"] == "molecule":
                ibm_mols += 1
            elif ibm_node["type"] == "reaction":
                if ibm_node["reaction_type"] == "wingSuzuki":
                    ibm_ws += 1
                elif ibm_node["reaction_type"] == "pentamerSuzuki":
                    ibm_ps += 1
                elif ibm_node["reaction_type"] == "deBoc":
                    ibm_deBoc += 1
                elif ibm_node["reaction_type"] == "BHA":
                    ibm_BHA += 1
                elif ibm_node["reaction_type"] == "SNAr":
                    ibm_SNAr += 1
                ibm_rxns += 1
            ibm_total += 1

        ibm_dict = {
            "ibm_mols": ibm_mols,
            "ibm_rxns": ibm_rxns,
            "ibm_ws": ibm_ws,
            "ibm_ps": ibm_ps,
            "ibm_deBoc": ibm_deBoc,
            "ibm_BHA": ibm_BHA,
            "ibm_SNAr": ibm_SNAr,
            "ibm_total": ibm_total,
        }
        print(ibm_dict)

    def count_tertbutyl_ethanoate(self):
        count = 0
        for ibm_node in self.ibm_rxn_data:
            if ibm_node["type"] == "molecule":
                if Chem.MolFromSmiles(ibm_node["SMILES"]).HasSubstructMatch(
                    Chem.MolFromSmiles("O=C(OC(C)(C)C)")
                ):
                    count += 1
            print(count)

    def compare_wS(self):
        ws_matches = 0
        ws_match_list = []
        for rs_node in self.routescore_data:
            if (
                rs_node["type"] == "reaction"
                and rs_node["reaction_type"] == "wingSuzuki"
            ):
                rs_product_smile = self.get_product_from_rxn(rs_node["rxn_SMILES"])
                for ibm_node in self.ibm_rxn_data:
                    if (
                        ibm_node["type"] == "reaction"
                        and ibm_node["reaction_type"] == "wingSuzuki"
                    ):
                        ibm_product_smile = self.get_product_from_rxn(
                            ibm_node["rxn_SMILES"]
                        )
                        if self.transforms_bmida_check(
                            rs_product_smile, ibm_product_smile
                        ):
                            ws_matches += 1
                            ws_match_list.append(ibm_product_smile)

        for ibm_node in self.ibm_rxn_data:
            if (
                ibm_node["type"] == "reaction"
                and ibm_node["reaction_type"] == "wingSuzuki"
            ):
                ibm_product_smile = self.get_product_from_rxn(ibm_node["rxn_SMILES"])
                if ibm_product_smile not in ws_match_list:
                    print(ibm_node["rxn_SMILES"])

        print(ws_matches)

    def compare_pS(self):
        ps_matches = 0
        ps_match_list = []
        count = 0
        for rs_node in self.routescore_data:
            if (
                rs_node["type"] == "reaction"
                and rs_node["reaction_type"] == "pentamerSuzuki"
            ):
                rs_product_smile = self.get_product_from_rxn(rs_node["rxn_SMILES"])
                for ibm_node in self.ibm_rxn_data:
                    if (
                        ibm_node["type"] == "reaction"
                        and ibm_node["reaction_type"] == "pentamerSuzuki"
                    ):
                        ibm_product_smile = self.get_product_from_rxn(
                            ibm_node["rxn_SMILES"]
                        )
                        if self.transforms_bmida_check(
                            rs_product_smile, ibm_product_smile
                        ):
                            ps_matches += 1
                            ps_match_list.append(ibm_product_smile)
                            print(ps_matches)

        for ibm_node in self.ibm_rxn_data:
            if (
                ibm_node["type"] == "reaction"
                and ibm_node["reaction_type"] == "pentamerSuzuki"
            ):
                ibm_product_smile = self.get_product_from_rxn(ibm_node["rxn_SMILES"])
                if ibm_product_smile not in ps_match_list:
                    print(ibm_node["rxn_SMILES"])

        print(ps_matches)

    def compare_deboc(self):
        deboc_matches = 0
        deboc_match_list = []
        count = 0
        for rs_node in self.routescore_data:
            if rs_node["type"] == "reaction" and rs_node["reaction_type"] == "deBoc":
                rs_product_smile = self.get_product_from_rxn(rs_node["rxn_SMILES"])
                for ibm_node in self.ibm_rxn_data:
                    if (
                        ibm_node["type"] == "reaction"
                        and ibm_node["reaction_type"] == "deBoc"
                    ):
                        ibm_product_smile = self.get_product_from_rxn(
                            ibm_node["rxn_SMILES"]
                        )
                        if self.transforms_bmida_check(
                            rs_product_smile, ibm_product_smile
                        ):
                            deboc_matches += 1
                            deboc_match_list.append(ibm_product_smile)
                            print(deboc_matches)

        for ibm_node in self.ibm_rxn_data:
            if ibm_node["type"] == "reaction" and ibm_node["reaction_type"] == "deBoc":
                ibm_product_smile = self.get_product_from_rxn(ibm_node["rxn_SMILES"])
                if ibm_product_smile not in deboc_match_list:
                    print(ibm_node["rxn_SMILES"])

        print(deboc_matches)

    def compare_BHA(self):
        BHA_matches = 0
        BHA_match_list = []
        count = 0
        for rs_node in self.routescore_data:
            if rs_node["type"] == "reaction" and rs_node["reaction_type"] == "BHA":
                rs_product_smile = self.get_product_from_rxn(rs_node["rxn_SMILES"])
                for ibm_node in self.ibm_rxn_data:
                    if (
                        ibm_node["type"] == "reaction"
                        and ibm_node["reaction_type"] == "BHA"
                    ):
                        ibm_product_smile = self.get_product_from_rxn(
                            ibm_node["rxn_SMILES"]
                        )
                        if self.transforms_bmida_check(
                            rs_product_smile, ibm_product_smile
                        ):
                            BHA_matches += 1
                            BHA_match_list.append(ibm_product_smile)
                            print(BHA_matches)

        for ibm_node in self.ibm_rxn_data:
            if ibm_node["type"] == "reaction" and ibm_node["reaction_type"] == "BHA":
                ibm_product_smile = self.get_product_from_rxn(ibm_node["rxn_SMILES"])
                if ibm_product_smile not in BHA_match_list:
                    print(ibm_node["rxn_SMILES"])

        print(BHA_matches)

    def compare_SNAr(self):
        SNAr_matches = 0
        SNAr_match_list = []
        count = 0
        for rs_node in self.routescore_data:
            if rs_node["type"] == "reaction" and rs_node["reaction_type"] == "SNAr":
                rs_product_smile = self.get_product_from_rxn(rs_node["rxn_SMILES"])
                for ibm_node in self.ibm_rxn_data:
                    if (
                        ibm_node["type"] == "reaction"
                        and ibm_node["reaction_type"] == "SNAr"
                    ):
                        ibm_product_smile = self.get_product_from_rxn(
                            ibm_node["rxn_SMILES"]
                        )
                        if self.transforms_bmida_check(
                            rs_product_smile, ibm_product_smile
                        ):
                            SNAr_matches += 1
                            SNAr_match_list.append(ibm_product_smile)
                            print(SNAr_matches)

        for ibm_node in self.ibm_rxn_data:
            if ibm_node["type"] == "reaction" and ibm_node["reaction_type"] == "SNAr":
                ibm_product_smile = self.get_product_from_rxn(ibm_node["rxn_SMILES"])
                if ibm_product_smile not in SNAr_match_list:
                    print(ibm_node["rxn_SMILES"])

        print(SNAr_matches)


compare = Compare(
    ROUTESCORE_JSON, FRAGCHEM_JSON, IBM_JSON_PATH, COMPARE_CSV, NOT_MATCH_CSV
)
# compare.mol_match()
# compare.count_mol_rxn()
# compare.compare_wS()
# compare.compare_pS()
# compare.compare_deboc()
# compare.compare_BHA()
compare.compare_SNAr()
# compare.count_tertbutyl_ethanoate()
# compare.check_duplicates()
# compare.full_props_comparison(FULL_PROPS_CSV)
