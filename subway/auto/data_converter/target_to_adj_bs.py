import pickle
import csv
import json

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
    print(data)


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
            a_SMILES = rows["a"]
            b_SMILES = rows["b"]
            ab_SMILES = rows["ab"]
            c_SMILES = rows["c"]
            pent_SMILES = rows["pentamer"]
            F_SMILES = rows["F"]
            carb_SMILES = rows["carbazole"]
            xphos_SMILES = "CC(C)C1=CC(=C(C(=C1)C(C)C)C2=CC(=CC=C2)P(C3CCCCC3)C4CCCCC4)C(C)C.C1=CC=C([C-]=C1)C2=CC=CC=C2N.Cl[Pd+]"
            kpo_SMILES = "[O-]P(=O)([O-])[O-].[K+].[K+].[K+]"
            pdba_SMILES = "C1=CC=C(C=C1)C=CC(=O)C=CC2=CC=CC=C2.C1=CC=C(C=C1)C=CC(=O)C=CC2=CC=CC=C2.C1=CC=C(C=C1)C=CC(=O)C=CC2=CC=CC=C2.[Pd].[Pd]"
            dave_SMILES = "CN(C)C1=CC=CC=C1C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4"
            naot_SMILES = "CC(C)(C)[O-].[Na+]"
            cs_SMILES = "C(=O)([O-])[O-].[Cs+].[Cs+]"
            kco_SMILES = "C(=O)([O-])[O-].[K+].[K+]"
            NH_SMILES = rows["N-H"]
            NBoc_SMILES = rows["N-Boc"]
            hal_SMILES = rows["halide"]

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

                    elif item["SMILES"] == NH_SMILES:
                        nh_m_node = item
                        nh_mol_id = nh_m_node["id"]

                    elif item["SMILES"] == hal_SMILES:
                        hal_m_node = item
                        hal_mol_id = hal_m_node["id"]

                    elif item["SMILES"] == F_SMILES:
                        F_m_node = item
                        F_mol_id = F_m_node["id"]

                    elif item["SMILES"] == cs_SMILES:
                        cs_m_node = item
                        cs_mol_id = cs_m_node["id"]

                    elif item["SMILES"] == kco_SMILES:
                        kco_m_node = item
                        kco_mol_id = kco_m_node["id"]

                    elif item["SMILES"] == carb_SMILES:
                        carb_m_node = item
                        carb_mol_id = carb_m_node["id"]

                    elif item["SMILES"] == naot_SMILES:
                        naot_m_node = item
                        naot_mol_id = naot_m_node["id"]

                    elif item["SMILES"] == pdba_SMILES:
                        pdba_m_node = item
                        pdba_mol_id = pdba_m_node["id"]

                    elif item["SMILES"] == dave_SMILES:
                        dave_m_node = item
                        dave_mol_id = dave_m_node["id"]

                    elif item["SMILES"] == NBoc_SMILES:
                        nboc_m_node = item
                        nboc_mol_id = nboc_m_node["id"]

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
                + NBoc_SMILES
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
            adj_list.append([rxn_id, nboc_mol_id])

            # deboc
            rxn_SMILES = NBoc_SMILES + ">" + kco_SMILES + ">" + NH_SMILES
            for item in jsonReader:
                if item["type"] == "reaction":
                    if item["rxn_SMILES"] == rxn_SMILES:
                        rxn_node = item
            rxn_id = rxn_node["id"]
            adj_list.append([nboc_mol_id, rxn_id])
            adj_list.append([kco_mol_id, rxn_id])
            adj_list.append([rxn_id, nh_mol_id])

            # BHA
            rxn_SMILES = (
                NH_SMILES
                + "."
                + hal_SMILES
                + ">"
                + pdba_SMILES
                + "."
                + dave_SMILES
                + "."
                + naot_SMILES
                + ">"
                + F_SMILES
            )
            for item in jsonReader:
                if item["type"] == "reaction":
                    if item["rxn_SMILES"] == rxn_SMILES:
                        rxn_node = item

            rxn_id = rxn_node["id"]
            adj_list.append([nh_mol_id, rxn_id])
            adj_list.append([hal_mol_id, rxn_id])
            adj_list.append([dave_mol_id, rxn_id])
            adj_list.append([naot_mol_id, rxn_id])
            adj_list.append([pdba_mol_id, rxn_id])
            adj_list.append([rxn_id, F_mol_id])

            # SNAr
            rxn_SMILES = (
                F_SMILES + "." + carb_SMILES + ">" + cs_SMILES + ">" + pent_SMILES
            )
            for item in jsonReader:
                if item["type"] == "reaction":
                    if item["rxn_SMILES"] == rxn_SMILES:
                        rxn_node = item

            rxn_id = rxn_node["id"]
            adj_list.append([F_mol_id, rxn_id])
            adj_list.append([carb_mol_id, rxn_id])
            adj_list.append([cs_mol_id, rxn_id])
            adj_list.append([rxn_id, pent_mol_id])

    file = open(pkl_path, "wb")
    # dump information to that file
    pickle.dump(adj_list, file)
    # close the file
    file.close()


target_to_adj(PKL_PATH, TGT_BS_FILEPATH, JSON_PATH)
# read_pkl(PKL_PATH)

