import csv
import json

# might need an add json file


def make_json(csvFilePath, jsonFilePath):
    data = []
    with open(csvFilePath, encoding="utf-8") as csvf:
        csvReader = csv.DictReader(csvf)
        # Convert each row into a dictionary
        # and add it to data
        id_count = 0
        for rows in csvReader:
            rows["id"] = id_count
            id_count += 1
            if rows["Manual?"] == "Yes":
                rows["manual"] = True
            else:
                rows["manual"] = False
            data.append(rows)
    # Open a json writer, and use the json.dumps()
    # function to dump data
    with open(jsonFilePath, "w", encoding="utf-8") as jsonf:
        jsonf.write(json.dumps(data, indent=4))


# action
# csvFilePath = r'C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Inventory\Inventory.csv'
# csvFilePath = r'C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Rxn_Type\Buchwald_deprotection.csv'
# jsonFilePath = r'C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_graphs\subway\data\molecule_node_auto.json'


def add_json(csvFilePath, jsonFilePath):
    data = []
    with open(jsonFilePath, encoding="utf-8") as csvf:
        csvReader = json.load(csvf)
        # Convert each row into a dictionary
        # and add it to data
        id_count = 0
        for rows in csvReader:
            id_count += 1

    with open(csvFilePath, encoding="utf-8") as csvf:
        csvReader = csv.DictReader(csvf)
        for rows in csvReader:
            rows["id"] = id_count
            id_count += 1
            data.append(rows)

    with open(jsonFilePath, "a", encoding="utf-8") as jsonf:
        jsonf.write(json.dumps(data, indent=4))
        jsonf.close()


sum_suzuki_FilePath = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Rxn_Summs\Suzuki_summary.csv"
sum_BHA_FilePath = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Rxn_Summs\Buchwald_summary.csv"
sum_deboc_FilePath = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Rxn_Summs\Buchwald_deprotection_summary.csv"
sum_SNAr_FilePath = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Rxn_Summs\SNAr_summary.csv"
tgtFilePath = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Targets\target_Base_BLANK.csv"
jsonFilePath = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_graphs\subway\data\reaction_node_auto.json"


def sum_properties(sum_FilePath):
    with open(sum_FilePath, encoding="utf-8") as csvf:
        csvReader = csv.DictReader(csvf)
        for rows in csvReader:
            tH = rows["t_H"]
            tM = rows["t_M"]
            yld = rows["yield"]
            scale = rows["scale"]
            n_parr = rows["n_parr"]
    return tH, tM, yld, scale, n_parr


def add_rxn_json(
    sum_suzuki_FilePath,
    sum_BHA_FilePath,
    sum_deboc_FilePath,
    sum_SNAr_FilePath,
    jsonFilePath,
    tgtFilePath,
):
    data = []
    with open(jsonFilePath, encoding="utf-8") as csvf:
        csvReader = json.load(csvf)
        # Convert each row into a dictionary
        # and add it to data
        id_count = 0
        for rows in csvReader:
            id_count += 1
        # print(id_count)

    # RS_Base
    # Step 1
    tH, tM, yld, scale, n_parr = sum_properties(sum_suzuki_FilePath)
    with open(tgtFilePath, encoding="utf-8") as csvf:
        csvReader = csv.DictReader(csvf)
        for rows in csvReader:
            temp_data = {}
            temp_data["id"] = id_count
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
            temp_data["target"] = "Base"
            temp_data["type"] = "reaction"
            temp_data["reaction_type"] = "wingSuzuki"
            id_count += 1
            data.append(temp_data)

    # Step 2
    tH, tM, yld, scale, n_parr = sum_properties(sum_suzuki_FilePath)
    with open(tgtFilePath, encoding="utf-8") as csvf2:
        csvReader2 = csv.DictReader(csvf2)
        for rows in csvReader2:
            temp_data = {}
            temp_data["id"] = str(id_count) + "r"
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
                + rows["pentamer"]
            )
            temp_data["rxn_SMILES"] = rxn_SMILES
            temp_data["type"] = "Base"
            temp_data["rxn_type"] = "pentamerSuzuki"
            id_count += 1
            data.append(temp_data)
    """
    # Step 3
    tH, tM, yld, scale, n_parr = sum_properties(sum_deboc_FilePath)
    with open(tgtFilePath, encoding="utf-8") as csvf3:
        csvReader3 = csv.DictReader(csvf3)
        for rows in csvReader3:
            temp_data = {}
            temp_data["id"] = str(id_count) + "r"
            temp_data["tH"] = tH
            temp_data["tM"] = tM
            temp_data["yld"] = yld
            temp_data["scale"] = scale
            temp_data["n_parr"] = n_parr
            rxn_SMILES = (
                rows["N-Boc"] + ">" + "C(=O)([O-])[O-].[K+].[K+]" + ">" + rows["N-H"]
            )
            temp_data["rxn_SMILES"] = rxn_SMILES
            temp_data["type"] = "BS"
            temp_data["rxn_type"] = "deboc"
            id_count += 1
            data.append(temp_data)

    # Step 4
    tH, tM, yld, scale, n_parr = sum_properties(sum_BHA_FilePath)
    with open(tgtFilePath, encoding="utf-8") as csvf4:
        csvReader4 = csv.DictReader(csvf4)
        for rows in csvReader4:
            temp_data = {}
            temp_data["id"] = str(id_count) + "r"
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
                + rows["F"]
            )
            temp_data["rxn_SMILES"] = rxn_SMILES
            temp_data["type"] = "BS"
            temp_data["rxn_type"] = "BHA"
            id_count += 1
            data.append(temp_data)

    # Step 5
    tH, tM, yld, scale, n_parr = sum_properties(sum_SNAr_FilePath)
    with open(tgtFilePath, encoding="utf-8") as csvf5:
        csvReader5 = csv.DictReader(csvf5)
        for rows in csvReader5:
            temp_data = {}
            temp_data["id"] = str(id_count) + "r"
            temp_data["tH"] = tH
            temp_data["tM"] = tM
            temp_data["yld"] = yld
            temp_data["scale"] = scale
            temp_data["n_parr"] = n_parr
            rxn_SMILES = (
                rows["F"]
                + "."
                + rows["carbazole"]
                + ">"
                + "C(=O)([O-])[O-].[K+].[K+]"
                + ">"
                + rows["pentamer"]
            )
            temp_data["rxn_SMILES"] = rxn_SMILES
            temp_data["type"] = "BS"
            temp_data["rxn_type"] = "SNAr"
            id_count += 1
            data.append(temp_data)
    """

    with open(jsonFilePath, "a", encoding="utf-8") as jsonf:
        jsonf.write(json.dumps(data, indent=4))
        jsonf.close()


# add_json(csvFilePath, jsonFilePath)
add_rxn_json(
    sum_suzuki_FilePath,
    sum_BHA_FilePath,
    sum_deboc_FilePath,
    sum_SNAr_FilePath,
    jsonFilePath,
    tgtFilePath,
)

