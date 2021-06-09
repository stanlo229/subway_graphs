import pickle
import csv
import json
from rdkit import Chem

import pkg_resources

"""
Format: [node_from, node_to]
"""

JSON_PATH = pkg_resources.resource_filename("subway", "data/auto_nodes.json")
PKL_PATH = pkg_resources.resource_filename("subway", "data/adj_list.pkl")

TGT_BASE_FILEPATH = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Targets\targets_Base_BLANK.csv"
TGT_BS_FILEPATH = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Targets\targets_B-S_BLANK.csv"
TGT_BUCH_FILEPATH = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Targets\targets_Buch_BLANK.csv"
TGT_SB_FILEPATH = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Targets\targets_S-B_BLANK.csv"
TGT_SNAR_FILEPATH = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Targets\targets_SNAr_BLANK.csv"


def clear_pkl(pkl_path):
    empty = []
    file = open(pkl_path, "wb")
    pickle.dump(empty, file)
    file.close()


def read_pkl(pkl_path):
    file = open(pkl_path, "rb")
    data = pickle.load(file)
    file.close()
    print(len(data))


def target_to_adj(pkl_path, tgt_path, json_path):
    # open file
    file = open(pkl_path, "rb")
    data = pickle.load(file)
    file.close()

    # check if file is empty
    if len(data) != 0:
        adj_list = data
    else:
        adj_list = []

    with open(json_path, encoding="utf-8") as jsonf:
        jsonReader = json.load(jsonf)

    with open(tgt_path, encoding="utf-8") as csvf:
        csvReader = csv.DictReader(csvf)
        for rows in csvReader:
            # get SMILES of all molecules in each route
            # remove ones not in target file
            a_SMILES = Chem.CanonSmiles(rows["a"])
            b_SMILES = Chem.CanonSmiles(rows["b"])
            ab_SMILES = Chem.CanonSmiles(rows["ab"])
            c_SMILES = Chem.CanonSmiles(rows["c"])
            pent_SMILES = Chem.CanonSmiles(rows["pentamer"])
            xphos_SMILES = "CC(C)C1=CC(=C(C(=C1)C(C)C)C2=CC(=CC=C2)P(C3CCCCC3)C4CCCCC4)C(C)C.C1=CC=C([C-]=C1)C2=CC=CC=C2N.Cl[Pd+]"
            kpo_SMILES = "[O-]P(=O)([O-])[O-].[K+].[K+].[K+]"

            # find appropriate nodes
            for item in jsonReader:
                if item["type"] == "molecule":
                    if item["SMILES"] == a_SMILES:
                        a_m_node = item
                        a_mol_id = a_m_node["id"]
                    elif item["SMILES"] == b_SMILES:
                        b_m_node = item
                        b_mol_id = b_m_node["id"]
                    elif item["SMILES"] == ab_SMILES:
                        ab_m_node = item
                        ab_mol_id = ab_m_node["id"]
                    elif item["SMILES"] == xphos_SMILES:
                        xphos_m_node = item
                        xphos_mol_id = xphos_m_node["id"]
                    elif item["SMILES"] == kpo_SMILES:
                        kpo_m_node = item
                        kpo_mol_id = kpo_m_node["id"]
                    elif item["SMILES"] == c_SMILES:
                        c_m_node = item
                        c_mol_id = c_m_node["id"]
                    elif item["SMILES"] == pent_SMILES:
                        pent_m_node = item
                        pent_mol_id = pent_m_node["id"]

            # wingSuzuki
            rxn_SMILES = (
                a_SMILES
                + "."
                + b_SMILES
                + ">"
                + xphos_SMILES
                + "."
                + kpo_SMILES
                + ">"
                + ab_SMILES
            )
            for item in jsonReader:
                if item["type"] == "reaction":
                    if item["rxn_SMILES"] == rxn_SMILES:
                        rxn_node = item
                        rxn_id = rxn_node["id"]
                        adj_list.append([a_mol_id, rxn_id])
                        adj_list.append([b_mol_id, rxn_id])
                        adj_list.append([xphos_mol_id, rxn_id])
                        adj_list.append([kpo_mol_id, rxn_id])
                        adj_list.append([rxn_id, ab_mol_id])

            # pentamerSuzuki
            rxn_SMILES = (
                ab_SMILES
                + "."
                + c_SMILES
                + ">"
                + xphos_SMILES
                + "."
                + kpo_SMILES
                + ">"
                + pent_SMILES
            )
            for item in jsonReader:
                if item["type"] == "reaction":
                    if item["rxn_SMILES"] == rxn_SMILES:
                        rxn_node = item
                        rxn_id = rxn_node["id"]
                        adj_list.append([ab_mol_id, rxn_id])
                        adj_list.append([c_mol_id, rxn_id])
                        adj_list.append([xphos_mol_id, rxn_id])
                        adj_list.append([kpo_mol_id, rxn_id])
                        adj_list.append([rxn_id, pent_mol_id])
            # deboc
            # BHA
            # SNAr
    file = open(pkl_path, "wb")
    # dump information to that file
    pickle.dump(adj_list, file)
    # close the file
    file.close()


target_to_adj(PKL_PATH, TGT_BASE_FILEPATH, JSON_PATH)
# read_pkl(PKL_PATH)
# clear_pkl(PKL_PATH)
