import json
import pandas as pd
import pkg_resources

ROUTESCORE_JSON = pkg_resources.resource_filename("subway", "data/auto_nodes.json")
FRAGCHEM_JSON = pkg_resources.resource_filename(
    "subway", "data/reaction_template_nodes.json"
)
COMPARE_CSV = pkg_resources.resource_filename("subway", "data/json_comparison.csv")
NOT_MATCH_CSV = pkg_resources.resource_filename("subway", "data/not_matched_debug.csv")
FULL_PROPS_CSV = pkg_resources.resource_filename(
    "subway", "data/subway_maps/full_props.csv"
)


class Compare:
    """ Class that contains functions for comparing the routescore data to the new graph data
    """

    def __init__(self, routescore_json, fragchem_json, compare_csv, not_match_csv):
        file = open(routescore_json, "rb")
        self.routescore_data = json.load(file)
        file.close()
        file = open(fragchem_json, "rb")
        self.fragchem_data = json.load(file)
        file.close()
        self.compare_path = compare_csv
        self.not_match_path = not_match_csv

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

    def check_duplicates(self):
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
        matches = 0
        full_props_df = pd.read_csv(full_props_path)
        for product in full_props_df["smiles"]:
            for node in self.fragchem_data:
                if node["type"] == "molecule":
                    if product == node["SMILES"]:
                        matches += 1
        print(matches)


compare = Compare(ROUTESCORE_JSON, FRAGCHEM_JSON, COMPARE_CSV, NOT_MATCH_CSV)
compare.smiles_match()
# compare.check_duplicates()
# compare.full_props_comparison(FULL_PROPS_CSV)

