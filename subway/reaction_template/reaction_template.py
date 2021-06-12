from rdkit import Chem
from rdkit.Chem import rdChemReactions
import json
import pkg_resources
from typing import List

from rdkit.Chem.Draw.IPythonConsole import _GetSubstructMatches

JSON_PATH = pkg_resources.resource_filename(
    "subway", "data/reaction_template_nodes.json"
)


class ReactionTemplates:
    """ 
    Class containing reaction templates that will create reaction nodes from molecule nodes
    """

    def __init__(self, json_path):
        file = open(json_path, "rb")
        self.json_data = json.load(file)
        file.close()

    def wingSuzuki(self, reactants_smiles: List[str]):
        """ Function that will contain reaction information for any wingSuzuki reaction.
        """
        # reaction-specific information
        reagent_smiles = [
            "[O-]P(=O)([O-])[O-].[K+].[K+].[K+]",
            "CC(C)C1=CC(=C(C(=C1)C(C)C)C2=CC(=CC=C2)P(C3CCCCC3)C4CCCCC4)C(C)C.C1=CC=C([C-]=C1)C2=CC=CC=C2N.Cl[Pd+]",
        ]
        patts: list = [
            {
                "name": "smilesB",
                "substructure": [Chem.MolFromSmarts("cBr"), Chem.MolFromSmarts("cI"),],
                "eq": 1,
            },
            {
                "name": "smilesA",
                "substructure": [
                    Chem.MolFromSmarts("cB(O)O")
                ],  # explicitly define the H's to avoid MIDA boronate / BMIDA
                "eq": 3,
            },
        ]
        t_H = 2.5
        t_M = 28
        n_parr = 48
        yld = 1
        scale = 0.0001

        # process reactants to create reaction_node
        # find limiting reagent by calculating reaction sites
        rxn_sites = 0
        for rct in reactants_smiles:
            mol = Chem.MolFromSmiles(rct)
            for structure in patts:
                for substruct in structure["substructure"]:
                    if len(mol.GetSubstructMatches(substruct)) > rxn_sites:
                        rxn_sites = len(mol.GetSubstructMatches(substruct))
                        limiting_rct_smiles = rct

        # creating reaction node
        reaction_node = {}
        reaction_smile = ""
        # creating Reaction SMILE
        # add reactants
        index = 0
        while index < len(reactants_smiles):
            reaction_smile += reactants_smiles[index]
            index += 1
            if index != len(reactants_smiles):
                reaction_smile += "."
        reaction_smile += ">"

        # add reagents
        index = 0
        while index < len(reagent_smiles):
            reaction_smile += reagent_smiles[index]
            index += 1
            if index != len(reagent_smiles):
                reaction_smile += "."
        reaction_smile += ">"

        # add product

        # setup reaction node
        reaction_node["rxn_SMILES"] = reaction_smile
        reaction_node["t_H"] = t_H
        reaction_node["t_M"] = t_M
        reaction_node["n_parr"] = n_parr
        reaction_node["yield"] = yld
        reaction_node["scale"] = scale
        reaction_node["type"] = "reaction"
        reaction_node["reaction_type"] = "wingSuzuki"

        # add reaction node to json

    def pentamerSuzuki(self, reactants_smiles: List[str]):
        """ Function that will contain reaction information for any pentamerSuzuki reaction.
        """
        # reaction-specific information
        reagent_smiles = [
            "[O-]P(=O)([O-])[O-].[K+].[K+].[K+]",
            "CC(C)C1=CC(=C(C(=C1)C(C)C)C2=CC(=CC=C2)P(C3CCCCC3)C4CCCCC4)C(C)C.C1=CC=C([C-]=C1)C2=CC=CC=C2N.Cl[Pd+]",
        ]
        patts: list = [
            {
                "name": "smilesAB",
                "substructure": [Chem.MolFromSmarts("c[B-]12OC(C[N+]1(C)CC(O2)=O)=O")],
                "eq": 3,
            },
            {
                "name": "smilesC",
                "substructure": [Chem.MolFromSmarts("cBr"), Chem.MolFromSmarts("cI"),],
                "eq": 1,
            },
        ]
        t_H = 2.5
        t_M = 28
        n_parr = 48

        # process reactants to create reaction_node

        reaction_node = {}

    def deBoc(self, reactants_smiles: List[str]):
        """ Function that will contain reaction information for any pentamerSuzuki reaction.
        """
        # reaction-specific information
        reagent_smiles = ["C(=O)([O-])[O-].[K+].[K+]"]
        patts: list = [
            {
                "name": "smilesNBoc",
                "substructure": [Chem.MolFromSmiles("CC(C)(C)OC(=O)")],
                "eq": 1,
            }
        ]
        t_H = 6.5
        t_M = 0
        n_parr = 1

        # process reactants to create reaction_node

        reaction_node = {}

    def BHA(self, reactants_smiles: List[str]):
        """ Function that will contain reaction information for any pentamerSuzuki reaction.
        """
        # reaction-specific information
        reagent_smiles = [
            "C1=CC=C(C=C1)C=CC(=O)C=CC2=CC=CC=C2.C1=CC=C(C=C1)C=CC(=O)C=CC2=CC=CC=C2.C1=CC=C(C=C1)C=CC(=O)C=CC2=CC=CC=C2.[Pd].[Pd]",
            "CN(C)C1=CC=CC=C1C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4",
            "CC(C)(C)[O-].[Na+]",
        ]
        patts: list = [
            {
                "name": "smilesNH",
                "substructure": [Chem.MolFromSmarts("[nH]")],
                "eq": 1,
            },
            {"name": "smilesX", "substructure": [Chem.MolFromSmiles("Br")], "eq": 3,},
        ]
        t_H = 6
        t_M = 0
        n_parr = 1

        # process reactants to create reaction_node

        reaction_node = {}

    def SNAr(self, reactants_smiles: List[str]):
        """ Function that will contain reaction information for any pentamerSuzuki reaction.
        """
        # reaction-specific information
        reagent_smiles = ["C(=O)([O-])[O-].[Cs+].[Cs+]"]
        patts: list = [
            {"name": "smilesAr", "substructure": [Chem.MolFromSmarts("cF")], "eq": 1,},
            {
                "name": "smilesNu",
                "substructure": [Chem.MolFromSmiles("c1ccc2c(c1)[nH]c1ccccc12")],
                "eq": 2,
            },
        ]
        t_H = 6.5
        t_M = 0
        n_parr = 1

        # process reactants to create reaction_node
        reaction_node = {}


rxn_template = ReactionTemplates(JSON_PATH)
rxn_template.wingSuzuki(["OB(O)c1ccco1", "C[N+]12CC(=O)O[B-]1(c1cccnc1Br)OC(=O)C2"])
