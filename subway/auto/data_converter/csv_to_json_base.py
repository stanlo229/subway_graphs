import csv
import json
import pickle
from rdkit import Chem

import pkg_resources

JSON_PATH = pkg_resources.resource_filename("subway", "data/auto_nodes.json")

SUM_SUZUKI_FILEPATH = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Rxn_Summs\Suzuki_summary.pkl"
SUM_BHA_FILEPATH = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Rxn_Summs\Buchwald_summary.pkl"
SUM_DEBOC_FILEPATH = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Rxn_Summs\Buchwald_deprotection_summary.pkl"
SUM_SNAR_FILEPATH = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Rxn_Summs\SNAr_summary.pkl"

TGT_BASE_FILEPATH = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Targets\targets_Base_BLANK.csv"
TGT_BS_FILEPATH = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Targets\targets_B-S_BLANK.csv"
TGT_BUCH_FILEPATH = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Targets\targets_Buch_BLANK.csv"
TGT_SB_FILEPATH = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Targets\targets_S-B_BLANK.csv"
TGT_SNAR_FILEPATH = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Targets\targets_SNAr_BLANK.csv"


def check_json(json_path, reaction_smile):
    with open(json_path, encoding="utf-8") as file:
        jsonReader = json.load(file)
        for node in jsonReader:
            if node["type"] == "reaction":
                if node["rxn_SMILES"] == reaction_smile:
                    return False
                else:
                    continue
        return True


def sum_properties(sum_path):
    file = open(sum_path, "rb")
    sum_data = pickle.load(file)
    file.close()
    return (
        sum_data["t_H"],
        sum_data["t_M"],
        sum_data["yield"],
        sum_data["scale"],
        sum_data["n_parr"],
    )


# add reaction json to json
def add_rxn_json(
    sum_suzuki_path, sum_bha_path, sum_deboc_path, sum_snar_path, tgt_path, json_path
):
    # count number of nodes in json
    with open(json_path, encoding="utf-8") as file:
        jsonReader = json.load(file)
        rxn_id = len(jsonReader)  # 3861

    data = []
    # wingSuzuki
    rxn_repo = []
    tH, tM, yld, scale, n_parr = sum_properties(sum_suzuki_path)
    with open(tgt_path, encoding="utf-8") as csvf:
        csvReader = csv.DictReader(csvf)
        for rows in csvReader:
            temp_data = {}
            temp_data["tH"] = tH
            temp_data["tM"] = tM
            temp_data["yield"] = yld
            temp_data["scale"] = scale
            temp_data["n_parr"] = n_parr
            rxn_SMILES = (
                Chem.CanonSmiles(rows["a"])
                + "."
                + Chem.CanonSmiles(rows["b"])
                + ">"
                + "CC(C)C1=CC(=C(C(=C1)C(C)C)C2=CC(=CC=C2)P(C3CCCCC3)C4CCCCC4)C(C)C.C1=CC=C([C-]=C1)C2=CC=CC=C2N.Cl[Pd+]"
                + "."
                + "[O-]P(=O)([O-])[O-].[K+].[K+].[K+]"
                + ">"
                + Chem.CanonSmiles(rows["ab"])
            )
            temp_data["rxn_SMILES"] = rxn_SMILES
            if rxn_SMILES not in rxn_repo:
                rxn_repo.append(rxn_SMILES)
                temp_data["id"] = rxn_id
                temp_data["type"] = "reaction"
                temp_data["reaction_type"] = "wingSuzuki"
                rxn_id += 1
                data.append(temp_data)

    # pentamerSuzuki
    rxn_repo = []
    tH, tM, yld, scale, n_parr = sum_properties(sum_suzuki_path)
    with open(tgt_path, encoding="utf-8") as csvf2:
        csvReader2 = csv.DictReader(csvf2)
        for rows in csvReader2:
            temp_data = {}
            temp_data["tH"] = tH
            temp_data["tM"] = tM
            temp_data["yld"] = yld
            temp_data["scale"] = scale
            temp_data["n_parr"] = n_parr
            rxn_SMILES = (
                Chem.CanonSmiles(rows["ab"])
                + "."
                + Chem.CanonSmiles(rows["c"])
                + ">"
                + "CC(C)C1=CC(=C(C(=C1)C(C)C)C2=CC(=CC=C2)P(C3CCCCC3)C4CCCCC4)C(C)C.C1=CC=C([C-]=C1)C2=CC=CC=C2N.Cl[Pd+]"
                + "."
                + "[O-]P(=O)([O-])[O-].[K+].[K+].[K+]"
                + ">"
                + Chem.CanonSmiles(rows["pentamer"])
            )
            temp_data["rxn_SMILES"] = rxn_SMILES
            if rxn_SMILES not in rxn_repo:
                rxn_repo.append(rxn_SMILES)
                temp_data["id"] = rxn_id
                temp_data["type"] = "reaction"
                temp_data["rxn_type"] = "pentamerSuzuki"
                rxn_id += 1
                data.append(temp_data)

    with open(json_path, "a", encoding="utf-8") as jsonf:
        jsonf.write(json.dumps(data, indent=4))
        jsonf.close()


add_rxn_json(
    SUM_SUZUKI_FILEPATH,
    SUM_BHA_FILEPATH,
    SUM_DEBOC_FILEPATH,
    SUM_SNAR_FILEPATH,
    TGT_BASE_FILEPATH,
    JSON_PATH,
)

