import pickle
import csv
import json


tgt_SNAr_FilePath = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Targets\target_SNAr_BLANK.csv"
tgt_Base_FilePath = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Targets\target_Base_BLANK.csv"
tgt_BS_FilePath = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Targets\target_B-S_BLANK.csv"
tgt_SB_FilePath = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Targets\target_S-B_BLANK.csv"
tgt_Buch_FilePath = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Targets\target_Buch_BLANK.csv"

mol_jsonFilePath = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_graphs\subway\data\molecule_node_auto.json"
rxn_jsonFilePath = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_graphs\subway\data\reaction_node_auto.json"

PIK = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_graphs\subway\data\auto_adj_list.pkl"


def target_to_adj(PIK, tgtFile, mol_jsonFile, rxn_jsonFile):
    # open file
    file = open(PIK, "rb")
    data = pickle.load(file)
    file.close()
    if len(data) != 0:
        adj_from = data[0]
        adj_to = data[1]
    else:
        adj_from = []
        adj_to = []
    with open(mol_jsonFile, encoding="utf-8") as csvf:
        csvReader = json.load(csvf)

    with open(rxn_jsonFile, encoding="utf-8") as csvf3:
        csvReader3 = json.load(csvf3)

    with open(tgtFile, encoding="utf-8") as csvf2:
        csvReader2 = csv.DictReader(csvf2)
        for rows in csvReader2:
            # get SMILES of all molecules in each route
            a_SMILES = rows["a"]
            b_SMILES = rows["b"]
            ab_SMILES = rows["ab"]
            c_SMILES = rows["c"]
            pent_SMILES = rows["pentamer"]
            F_SMILES = rows["F"]
            carb_SMILES = rows["carbazole"]
            # NH_SMILES = rows["N-H"]
            # NBoc_SMILES = rows["N-Boc"]
            # hal_SMILES = rows["halide"]

            # match SMILES to SMILES in .json for each rxn
            # Step 1
            a_m_node = next(
                (item for item in csvReader if item["SMILES"] == a_SMILES), None
            )
            a_mol_id = a_m_node["id"]
            b_m_node = next(
                (item for item in csvReader if item["SMILES"] == b_SMILES), None
            )
            b_mol_id = b_m_node["id"]
            ab_m_node = next(
                (item for item in csvReader if item["SMILES"] == ab_SMILES), None
            )
            ab_mol_id = ab_m_node["id"]
            rxn_SMILES = a_SMILES + "." + b_SMILES + ">>" + ab_SMILES
            rxn_node = next(
                (item for item in csvReader3 if item["rxn_SMILES"] == rxn_SMILES), None
            )
            rxn_id = str(rxn_node["id"])
            rxn_id = rxn_id.replace("r", "")
            rxn_id = int(rxn_id)
            adj_from.append(a_mol_id)
            adj_from.append(b_mol_id)
            adj_from.append(ab_mol_id)
            for i in range(3):
                adj_to.append(rxn_id)

            # Step 2
            ab_m_node = next(
                (item for item in csvReader if item["SMILES"] == ab_SMILES), None
            )
            ab_mol_id = ab_m_node["id"]
            c_m_node = next(
                (item for item in csvReader if item["SMILES"] == c_SMILES), None
            )
            c_mol_id = c_m_node["id"]
            F_m_node = next(
                (item for item in csvReader if item["SMILES"] == F_SMILES), None
            )
            F_mol_id = F_m_node["id"]
            rxn_SMILES = ab_SMILES + "." + c_SMILES + ">>" + F_SMILES
            rxn_node = next(
                (item for item in csvReader3 if item["rxn_SMILES"] == rxn_SMILES), None
            )
            rxn_id = str(rxn_node["id"])
            rxn_id = rxn_id.replace("r", "")
            rxn_id = int(rxn_id)
            adj_from.append(ab_mol_id)
            adj_from.append(c_mol_id)
            adj_from.append(F_mol_id)
            for i in range(3):
                adj_to.append(rxn_id)

            # Step 3
            F_m_node = next(
                (item for item in csvReader if item["SMILES"] == F_SMILES), None
            )
            F_mol_id = F_m_node["id"]
            carb_node = next(
                (item for item in csvReader if item["SMILES"] == carb_SMILES), None
            )
            carb_mol_id = carb_node["id"]
            pent_node = next(
                (item for item in csvReader if item["SMILES"] == pent_SMILES), None
            )
            pent_mol_id = pent_node["id"]
            rxn_SMILES = F_SMILES + "." + carb_SMILES + ">>" + pent_SMILES
            rxn_node = next(
                (item for item in csvReader3 if item["rxn_SMILES"] == rxn_SMILES), None
            )
            rxn_id = str(rxn_node["id"])
            rxn_id = rxn_id.replace("r", "")
            rxn_id = int(rxn_id)
            adj_from.append(F_mol_id)
            adj_from.append(carb_mol_id)
            adj_from.append(pent_mol_id)
            for i in range(3):
                adj_to.append(rxn_id)

            """
            # Step 4
            nboc_m_node = next((item for item in csvReader if item["SMILES"] == NBoc_SMILES), None)
            nboc_mol_id = nboc_m_node["id"]
            nh_node = next((item for item in csvReader if item["SMILES"] == NH_SMILES), None)
            nh_mol_id = nh_node["id"]
            rxn_SMILES = NBoc_SMILES + ">>" + NH_SMILES
            rxn_node = next((item for item in csvReader3 if item["rxn_SMILES"] == rxn_SMILES), None)
            rxn_id = str(rxn_node["id"])
            rxn_id = rxn_id.replace("r", "")
            rxn_id = int(rxn_id)
            adj_from.append(nboc_mol_id)
            adj_from.append(nh_mol_id)
            for i in range(2):
                adj_to.append(rxn_id)
            """
            """
            # Step 5
            nh_node = next((item for item in csvReader if item["SMILES"] == NH_SMILES), None)
            nh_mol_id = nh_node["id"]
            hal_node = next((item for item in csvReader if item["SMILES"] == hal_SMILES), None)
            hal_mol_id = hal_node["id"]
            pent_node = next((item for item in csvReader if item["SMILES"] == pent_SMILES), None)
            pent_mol_id = pent_node["id"]
            rxn_SMILES = NH_SMILES + "." + hal_SMILES + ">>" + pent_SMILES
            rxn_node = next((item for item in csvReader3 if item["rxn_SMILES"] == rxn_SMILES), None)
            rxn_id = str(rxn_node["id"])
            rxn_id = rxn_id.replace("r", "")
            rxn_id = int(rxn_id)
            adj_from.append(nh_mol_id)
            adj_from.append(hal_mol_id)
            adj_from.append(pent_mol_id)
            for i in range(3):
                adj_to.append(rxn_id)
            """
    # append data
    data = []
    data.append(adj_from)
    data.append(adj_to)
    # open a file, where you ant to store the data
    file = open(PIK, "wb")
    # dump information to that file
    pickle.dump(data, file)
    # close the file
    file.close()


def clear_pkl(PIK):
    empty = []
    file = open(PIK, "wb")
    pickle.dump(empty, file)
    file.close()


def read_pkl(PIK):
    file = open(PIK, "rb")
    data = pickle.load(file)
    for i in data:
        print(i)
        # print("Number of connections:", len(i))
    file.close()


# target_to_adj(PIK, tgt_SNAr_FilePath, mol_jsonFilePath, rxn_jsonFilePath)
# clear_pkl(PIK)
read_pkl(PIK)

