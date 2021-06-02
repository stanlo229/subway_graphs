import csv
import json
import pickle

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
                rows["a"]
                + "."
                + rows["b"]
                + ">"
                + "CC(C)C1=CC(=C(C(=C1)C(C)C)C2=CC(=CC=C2)P(C3CCCCC3)C4CCCCC4)C(C)C.C1=CC=C([C-]=C1)C2=CC=CC=C2N.Cl[Pd+]"
                + "."
                + "[O-]P(=O)([O-])[O-].[K+].[K+].[K+]"
                + ">"
                + rows["ab"]
            )
            temp_data["rxn_SMILES"] = rxn_SMILES
            temp_data["id"] = rxn_id
            temp_data["type"] = "reaction"
            temp_data["target"] = "Buch"
            temp_data["reaction_type"] = "wingSuzuki"
            rxn_id += 1
            data.append(temp_data)

    # pentamerSuzuki
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
                rows["ab"]
                + "."
                + rows["c"]
                + ">"
                + "CC(C)C1=CC(=C(C(=C1)C(C)C)C2=CC(=CC=C2)P(C3CCCCC3)C4CCCCC4)C(C)C.C1=CC=C([C-]=C1)C2=CC=CC=C2N.Cl[Pd+]"
                + "."
                + "[O-]P(=O)([O-])[O-].[K+].[K+].[K+]"
                + ">"
                + rows["N-Boc"]
            )
            temp_data["rxn_SMILES"] = rxn_SMILES
            temp_data["id"] = rxn_id
            temp_data["type"] = "reaction"
            temp_data["target"] = "Buch"
            temp_data["rxn_type"] = "pentamerSuzuki"
            rxn_id += 1
            data.append(temp_data)

    # deboc
    tH, tM, yld, scale, n_parr = sum_properties(sum_deboc_path)
    with open(tgt_path, encoding="utf-8") as csvf3:
        csvReader3 = csv.DictReader(csvf3)
        for rows in csvReader3:
            temp_data = {}
            temp_data["tH"] = tH
            temp_data["tM"] = tM
            temp_data["yld"] = yld
            temp_data["scale"] = scale
            temp_data["n_parr"] = n_parr
            rxn_SMILES = (
                rows["N-Boc"] + ">" + "C(=O)([O-])[O-].[K+].[K+]" + ">" + rows["N-H"]
            )
            temp_data["rxn_SMILES"] = rxn_SMILES
            temp_data["id"] = rxn_id
            temp_data["type"] = "reaction"
            temp_data["target"] = "Buch"
            temp_data["rxn_type"] = "deboc"
            rxn_id += 1
            data.append(temp_data)

    # BHA
    tH, tM, yld, scale, n_parr = sum_properties(sum_bha_path)
    with open(tgt_path, encoding="utf-8") as csvf4:
        csvReader4 = csv.DictReader(csvf4)
        for rows in csvReader4:
            temp_data = {}
            temp_data["id"] = rxn_id
            temp_data["tH"] = tH
            temp_data["tM"] = tM
            temp_data["yld"] = yld
            temp_data["scale"] = scale
            temp_data["n_parr"] = n_parr
            rxn_SMILES = (
                rows["N-H"]
                + "."
                + rows["halide"]
                + ">"
                + "C1=CC=C(C=C1)C=CC(=O)C=CC2=CC=CC=C2.C1=CC=C(C=C1)C=CC(=O)C=CC2=CC=CC=C2.C1=CC=C(C=C1)C=CC(=O)C=CC2=CC=CC=C2.[Pd].[Pd]"
                + "."
                + "CN(C)C1=CC=CC=C1C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4"
                + "."
                + "CC(C)(C)[O-].[Na+]"
                + ">"
                + rows["pentamer"]
            )
            temp_data["rxn_SMILES"] = rxn_SMILES
            temp_data["id"] = rxn_id
            temp_data["type"] = "reaction"
            temp_data["target"] = "Buch"
            temp_data["rxn_type"] = "BHA"
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
    TGT_BUCH_FILEPATH,
    JSON_PATH,
)

