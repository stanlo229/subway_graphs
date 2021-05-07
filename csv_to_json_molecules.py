import csv
import json

def make_json(csvFilePath, jsonFilePath):
    data = {}
    with open(csvFilePath, encoding='utf-8') as csvf:
        csvReader = csv.DictReader(csvf)
        # Convert each row into a dictionary
        # and add it to data
        for rows in csvReader:
            key = rows['Frag_type']
            data[key] = rows
    # Open a json writer, and use the json.dumps()
    # function to dump data
    with open(jsonFilePath, 'w', encoding='utf-8') as jsonf:
        jsonf.write(json.dumps(data, indent=4))

# action
csvFilePath = r'Inventory\Inventory_MASTER.csv'
jsonFilePath = r'molecule_node_auto.json'

make_json(csvFilePath, jsonFilePath)