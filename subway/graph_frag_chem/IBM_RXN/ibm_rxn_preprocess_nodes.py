import rdkit
from rdkit import Chem
import json
import pkg_resources

IBM_JSON_PATH = pkg_resources.resource_filename(
    "subway", "data/reaction_template/ibm_reaction_template_nodes.json"
)

bmida_substructure = "B1OC(=O)CN(C)CC(=O)O1"
bmida_mol = Chem.MolFromSmiles(bmida_substructure)


def change_bmida(json_path):
    """Function that changes the smiles of the molecules in the json files 
    by converting the existing bmida structure to a SMILES RDKIT-readable format
    """
    file = open(json_path, "rb")
    json_data = json.load(file)
    file.close()
    new_json_data = []

    for node in json_data:
        if node["Block_type"] == "b":
            old_bmida_smiles = node["SMILES"]
            old_bmida_mol = Chem.MolFromSmiles(old_bmida_smiles)
            new_bmida_mol = Chem.ReplaceSubstructs(
                old_bmida_mol, bmida_mol, bmida_mol, replaceAll=True
            )[0]
            node["SMILES"] = Chem.MolToSmiles(new_bmida_mol)
        new_json_data.append(node)

    new_json_file = open(json_path, "w")
    json.dump(new_json_data, new_json_file)
    new_json_file.close()


change_bmida(IBM_JSON_PATH)

