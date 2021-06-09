import csv
import json
import pickle
from rdkit import Chem

import pkg_resources

"""
csv_path = pkg_resources.resource_filename(
    "subway_maps_repo", "Inventory/Inventory.csv"
)
"""
CSV_PATH = r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_maps_repo\Inventory\Inventory.csv"
JSON_PATH = pkg_resources.resource_filename("subway", "data/auto_nodes.json")

# add molecule from inventory to json
def add_mol_json(csv_path, json_path):
    data = []
    with open(csv_path, encoding="utf-8") as csvf:
        csvReader = csv.DictReader(csvf)
        id = 0
        for rows in csvReader:
            del rows["Unnamed: 0"]
            del rows["Unnamed: 0.1"]
            del rows["Unnamed: 0.1.1"]
            rows["id"] = id
            rows["type"] = "molecule"
            rows["eq_per_site"] = rows["Quantity"]
            rows["SMILES"] = Chem.CanonSmiles(rows["SMILES"])
            del rows["Quantity"]
            id += 1
            data.append(rows)

    with open(json_path, "w", encoding="utf-8") as jsonf:
        jsonf.write(json.dumps(data, indent=4))


add_mol_json(CSV_PATH, JSON_PATH)

