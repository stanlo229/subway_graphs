from typing import List

import pkg_resources
import copy
from itertools import combinations
from IPython.display import display
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, rdChemReactions
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw.IPythonConsole import _GetSubstructMatches
from rdkit.Chem.MolStandardize import rdMolStandardize

IPythonConsole.ipython_useSVG = True
from subway.auto.graph_search_auto import Calculate

# smiles list
# input X number of molecules smiles and output the product smiles
# make the output of this and ibm rxn the same. Easily switchable and testable.


class ProductMaker:
    """Class containing substructure patterns and functions required to create product from reactant molecules

    """

    def __init__(self):
        # substructures that are designed to be replaced and are
        # categorized by their reaction type and smiles identifier
        self.patts = [
            {
                "rxn_type": "wingSuzuki",
                "smiles_boronic_acid": [
                    Chem.MolFromSmiles(
                        "B1OC(C)(C)C(C)(C)O1"
                    ),  # for the tetramethyl dioxaborolane
                    Chem.MolFromSmiles("B(O)O"),
                ],
                "smiles_halo-BMIDA": [
                    Chem.MolFromSmarts("Br"),
                    Chem.MolFromSmarts("I"),
                ],
            },
            {
                "rxn_type": "pentamerSuzuki",
                "smiles_BMIDA": [
                    Chem.MolFromSmiles("B1OC(=O)CN(C)CC(=O)O1")
                ],  # same as Chem.MolFromSmarts("c[B-]12OC(C[N+]1(C)CC(O2)=O)=O") but different in rdkit
                "smiles_dihalide": [Chem.MolFromSmarts("Br"), Chem.MolFromSmarts("I"),],
            },
            {
                "rxn_type": "BHA",
                "smiles_NH": [Chem.MolFromSmarts("[nH]")],
                "smiles_halide": [Chem.MolFromSmarts("Br"),],
            },
            {
                "rxn_type": "deBoc",
                "smiles_NBoc": [Chem.MolFromSmiles("CC(C)(C)OC(=O)")],
            },
            {
                "rxn_type": "SNAr",
                "smiles_Nu": [Chem.MolFromSmarts("[nH]")],
                "smiles_ArF": [Chem.MolFromSmarts("F")],
            },
        ]

    def reactant_stoichiometry(self, reactant_smile: str, rxn_type: str):
        """Determines number of reaction sites for a reactant
        
        Parameters
        ----------
        - reactant_smile: SMILE of reactant
        - rxn_type: type of reaction
        
        Returns
        -------
        - rxn_sites: number of reaction sites
        """
        mol = Chem.MolFromSmiles(reactant_smile)
        rxn_sites = 0
        for patt in self.patts:
            if patt["rxn_type"] == rxn_type:
                keys = list(patt.keys())
                keys = keys[1:]
                for key in keys:
                    for substruct in patt[key]:
                        rxn_sites += len(mol.GetSubstructMatches(substruct))
        return rxn_sites

    def rct_to_pdt(self, reactant_nodes: List, rxn_type: str):
        """Determines product when given list of reactant smiles and reaction type
        
        Parameters
        ----------
        - reactant_nodes: list of reactant nodes
        - rxn_type: type of reaction
        
        Returns
        -------
        - product_list: list of SMILES of product molecule
        """
        # reactant smiles pre-process
        reactant_smiles = []
        for rct_node in reactant_nodes:
            reactant_smiles.append(rct_node["SMILES"])

        # identify the reactant smiles by rxn_type and smiles matching
        # use molecule classifier
        mol_info = []
        rxn_site_info = []
        for rct in reactant_nodes:
            mol_dict = {}
            mol_dict["smiles"] = rct["SMILES"]
            mol_dict["smilestype"] = Calculate().molecule_classifier(rct, rxn_type)[0]
            # use stoichiometry to figure out which has more reaction sites
            rxn_site_info.append(self.reactant_stoichiometry(rct["SMILES"], rxn_type))
            mol_dict["rdkit_mol"] = Chem.MolFromSmiles(rct["SMILES"])
            mol_info.append(mol_dict)

        # NOTE: works best with only 2 reactants.
        # find limiting reactant smiles with lowest reaction sites
        lim_rct_smile = reactant_smiles[rxn_site_info.index(min(rxn_site_info))]
        # use self.patts and delete smiles of limiting reagent (less reaction sites)

        # create product by deleting substructure and replacing it
        product_list = []  # list of combined molecules (products)
        # NOTE: I think the only solution is to make it reaction-specific, but talk to Martin!
        # NOTE: when replacing substructure, where does it go? Ensure that connections make sense. CHECK.
        # NOTE: double check above by making sure all nodes are matched in auto_nodes.json and reaction_template_nodes.json
        if rxn_type == "wingSuzuki" or rxn_type == "pentamerSuzuki":
            for patt in self.patts:
                if patt["rxn_type"] == rxn_type:
                    # find index of replacement and delete the substructure with less reaction sites
                    for mol in mol_info:
                        try:
                            substruct = patt[mol["smilestype"]]
                        except:
                            print("error: wrong smilestype")
                        if mol["smiles"] == lim_rct_smile:
                            for sub in substruct:
                                if mol["rdkit_mol"].HasSubstructMatch(
                                    sub
                                ):  # check if molecule has substructure in matched smiles type from patts
                                    # find index of substructure in molecule
                                    substruct_index_tuple = mol[
                                        "rdkit_mol"
                                    ].GetSubstructMatch(sub)
                                    if 0 in substruct_index_tuple:
                                        substruct_index = 0
                                    else:
                                        substruct_index = min(substruct_index_tuple) - 1
                                    mol["rdkit_mol"] = Chem.DeleteSubstructs(
                                        mol["rdkit_mol"], sub
                                    )
                                    substitute_grp = mol["rdkit_mol"]
                                    continue
                    for mol in mol_info:
                        try:
                            substruct = patt[mol["smilestype"]]
                        except:
                            print("error: wrong smilestype")
                        # replace smiles on the other but with the rest of the previously modified molecule
                        if mol["smiles"] != lim_rct_smile:
                            for sub in substruct:
                                if mol["rdkit_mol"].HasSubstructMatch(
                                    sub
                                ):  # check if molecule has substructure in matched smiles type from patts
                                    mol["rdkit_mol"] = Chem.ReplaceSubstructs(
                                        mol["rdkit_mol"],
                                        sub,
                                        substitute_grp,
                                        replaceAll=True,
                                        replacementConnectionPoint=substruct_index,
                                    )[0]
                                    product_list.append(mol["rdkit_mol"])
                                    continue
        elif rxn_type == "SNAr":
            # two types of structure in SNAr reactants
            substruct_paradifluorobenzene = Chem.MolFromSmiles("Fc1ccc(F)cc1")
            substruct_thiadiazole = Chem.MolFromSmiles("c2ccc1nsnc1c2")
            # needs better name (ask Martin)
            substruct_thiadiazole_dimethyl = Chem.MolFromSmiles("Cc3nc2cc1nsnc1cc2nc3C")
            substruct_furan = Chem.MolFromSmiles("c1ccoc1")
            substruct_diphenyl_ether = Chem.MolFromSmiles("c2ccc(Oc1ccccc1)cc2")

            for patt in self.patts:
                if patt["rxn_type"] == rxn_type:
                    # find index of replacement and delete the substructure with less reaction sites
                    for mol in mol_info:
                        # replace substructure with the entire nH containing group
                        nH_grp = Chem.MolFromSmarts("[nH]")
                        if mol["rdkit_mol"].HasSubstructMatch(nH_grp):
                            substitute_grp = mol["rdkit_mol"]
                            for atom in mol["rdkit_mol"].GetAtoms():
                                if (
                                    atom.GetAtomicNum() == 7
                                ):  # atomic number 7 corresponds to Nitrogen
                                    nH_index = atom.GetIdx()
                    for mol in mol_info:
                        try:
                            substruct = patt[mol["smilestype"]]
                        except:
                            print("error: wrong smilestype")
                        if not mol["rdkit_mol"].HasSubstructMatch(nH_grp):
                            for sub in substruct:
                                if mol["rdkit_mol"].HasSubstructMatch(
                                    sub
                                ):  # check if molecule has substructure in matched smiles type from patts
                                    # find index of nH group for proper attachment
                                    # replace only aromatic cF
                                    # arbitrary atomic element that is not in dataset
                                    temp_replace_grp = Chem.MolFromSmiles("[Rn]")
                                    atom_index_list = []
                                    for atom in mol["rdkit_mol"].GetAtoms():
                                        # check if F
                                        if atom.GetAtomicNum() == 9:
                                            # look at neighbor groups
                                            for neighbor in atom.GetNeighbors():
                                                # only aromatic C
                                                if neighbor.GetIsAromatic():
                                                    atom_index_list.append(
                                                        atom.GetIdx()
                                                    )
                                    # creates all the necessary combinations
                                    if len(atom_index_list) == 2:
                                        for index in atom_index_list:
                                            atom = mol["rdkit_mol"].GetAtomWithIdx(
                                                index
                                            )
                                            atom.SetAtomicNum(
                                                86
                                            )  # change to a new structure that is not present in dataset
                                        mol["rdkit_mol"] = Chem.ReplaceSubstructs(
                                            mol["rdkit_mol"],
                                            temp_replace_grp,
                                            substitute_grp,
                                            replaceAll=True,
                                            replacementConnectionPoint=nH_index,
                                        )[0]
                                        for atom in mol["rdkit_mol"].GetAtoms():
                                            if atom.GetAtomicNum() == 7:
                                                atom.SetNumExplicitHs(0)
                                        product_list.append(mol["rdkit_mol"])

                                    elif mol["rdkit_mol"].HasSubstructMatch(
                                        substruct_furan
                                    ):
                                        pairings = []
                                        index = 0
                                        while index < len(atom_index_list):
                                            pairings.append(
                                                [
                                                    atom_index_list[index],
                                                    atom_index_list[index + 1],
                                                ]
                                            )
                                            index += 2

                                        for j in range(len(pairings)):
                                            comb = list(combinations(pairings, j + 1))
                                            for comb_pairings in comb:
                                                temp_pairings = []
                                                for pair in comb_pairings:
                                                    temp_pairings.extend(pair)
                                                # new copy of mol so original is not modified
                                                new_mol_ArF = copy.copy(
                                                    mol["rdkit_mol"]
                                                )
                                                for atom_index in temp_pairings:
                                                    new_mol_ArF.GetAtomWithIdx(
                                                        atom_index
                                                    ).SetAtomicNum(86)
                                                new_mol_ArF = Chem.ReplaceSubstructs(
                                                    new_mol_ArF,
                                                    temp_replace_grp,
                                                    substitute_grp,
                                                    replaceAll=True,
                                                    replacementConnectionPoint=nH_index,
                                                )[0]
                                                for atom in new_mol_ArF.GetAtoms():
                                                    if atom.GetAtomicNum() == 7:
                                                        atom.SetNumExplicitHs(0)
                                                product_list.append(new_mol_ArF)

                                    elif (
                                        mol["rdkit_mol"].HasSubstructMatch(
                                            substruct_paradifluorobenzene
                                        )
                                        or mol["rdkit_mol"].HasSubstructMatch(
                                            substruct_thiadiazole_dimethyl
                                        )
                                        or mol["rdkit_mol"].HasSubstructMatch(
                                            substruct_diphenyl_ether
                                        )
                                    ):
                                        # atom index pairings go around the structure.
                                        # So, if there are 6 total F's. The 1st and 4th F's are paired.
                                        divider = len(atom_index_list) // 2
                                        # pairings of F groups
                                        pairings = []
                                        for i in range(int(len(atom_index_list) / 2)):
                                            pairings.append(
                                                [
                                                    atom_index_list[i],
                                                    atom_index_list[i + divider],
                                                ]
                                            )
                                        for j in range(len(pairings)):
                                            comb = list(combinations(pairings, j + 1))
                                            for comb_pairings in comb:
                                                temp_pairings = []
                                                for pair in comb_pairings:
                                                    temp_pairings.extend(pair)
                                                # new copy of mol so original is not modified
                                                new_mol_ArF = copy.copy(
                                                    mol["rdkit_mol"]
                                                )
                                                for atom_index in temp_pairings:
                                                    new_mol_ArF.GetAtomWithIdx(
                                                        atom_index
                                                    ).SetAtomicNum(86)
                                                new_mol_ArF = Chem.ReplaceSubstructs(
                                                    new_mol_ArF,
                                                    temp_replace_grp,
                                                    substitute_grp,
                                                    replaceAll=True,
                                                    replacementConnectionPoint=nH_index,
                                                )[0]
                                                for atom in new_mol_ArF.GetAtoms():
                                                    if atom.GetAtomicNum() == 7:
                                                        atom.SetNumExplicitHs(0)
                                                product_list.append(new_mol_ArF)

                                    elif mol["rdkit_mol"].HasSubstructMatch(
                                        substruct_thiadiazole
                                    ):
                                        pairings = []
                                        first_index = 0
                                        last_index = len(atom_index_list) - 1
                                        for i in range(int(len(atom_index_list) / 2)):
                                            pairings.append(
                                                [
                                                    atom_index_list[first_index],
                                                    atom_index_list[last_index],
                                                ]
                                            )
                                            first_index += 1
                                            last_index -= 1
                                        for j in range(len(pairings)):
                                            comb = list(combinations(pairings, j + 1))
                                            for comb_pairings in comb:
                                                temp_pairings = []
                                                for pair in comb_pairings:
                                                    temp_pairings.extend(pair)
                                                # new copy of mol so original is not modified
                                                new_mol_ArF = copy.copy(
                                                    mol["rdkit_mol"]
                                                )
                                                for atom_index in temp_pairings:
                                                    new_mol_ArF.GetAtomWithIdx(
                                                        atom_index
                                                    ).SetAtomicNum(86)
                                                new_mol_ArF = Chem.ReplaceSubstructs(
                                                    new_mol_ArF,
                                                    temp_replace_grp,
                                                    substitute_grp,
                                                    replaceAll=True,
                                                    replacementConnectionPoint=nH_index,
                                                )[0]
                                                for atom in new_mol_ArF.GetAtoms():
                                                    if atom.GetAtomicNum() == 7:
                                                        atom.SetNumExplicitHs(0)
                                                product_list.append(new_mol_ArF)

        elif rxn_type == "BHA":
            for patt in self.patts:
                if patt["rxn_type"] == rxn_type:
                    # find index of replacement and delete the substructure with less reaction sites
                    for mol in mol_info:
                        # replace indole substructure with combined structure (indole + 2-bromopyrazine)
                        nH_grp = Chem.MolFromSmarts("[nH]")
                        if mol["rdkit_mol"].HasSubstructMatch(nH_grp):
                            indole_grp = Chem.MolFromSmiles("c2ccc1[nH]ccc1c2")
                            substitute_grp = Chem.MolFromSmiles(
                                "c1ccc2c(c1)ccn2c3cnccn3"
                            )
                            mol["rdkit_mol"] = Chem.ReplaceSubstructs(
                                mol["rdkit_mol"],
                                indole_grp,
                                substitute_grp,
                                replaceAll=True,
                            )[0]
                            for atom in mol["rdkit_mol"].GetAtoms():
                                if atom.GetAtomicNum() == 7:
                                    atom.SetNumExplicitHs(0)
                            product_list.append(mol["rdkit_mol"])

        elif rxn_type == "deBoc":
            # there should only be one reactant for this reaction
            # find the necessary substructure and remove it
            for patt in self.patts:
                if patt["rxn_type"] == rxn_type:
                    for mol in mol_info:
                        try:
                            substruct = patt[mol["smilestype"]]
                        except:
                            print("error: wrong smilestype")
                        for sub in substruct:
                            if mol["rdkit_mol"].HasSubstructMatch(sub):
                                mol["rdkit_mol"] = Chem.DeleteSubstructs(
                                    mol["rdkit_mol"], sub
                                )
                                # get index of indole substructure
                                indole_grp = Chem.MolFromSmiles("c2ccc1[nH]ccc1c2")
                                carbazole_grp = Chem.MolFromSmiles(
                                    "c1ccc3c(c1)[nH]c2ccccc23"
                                )
                                indole_index_tuple = mol[
                                    "rdkit_mol"
                                ].GetSubstructMatches(indole_grp)
                                carbazole_index_tuple = mol[
                                    "rdkit_mol"
                                ].GetSubstructMatches(carbazole_grp)
                                indole_index_list = []
                                carbazole_index_list = []
                                for tuple in indole_index_tuple:
                                    indole_index_list.extend(list(tuple))
                                for tuple in carbazole_index_tuple:
                                    carbazole_index_list.extend(list(tuple))
                                substruct_index_list = [
                                    item
                                    for item in indole_index_list
                                    if item not in carbazole_index_list
                                ]

                            for atom in mol["rdkit_mol"].GetAtoms():
                                if (
                                    atom.GetIsAromatic() == True  # aromatic nitrogen
                                    and atom.GetAtomicNum() == 7
                                    and atom.GetIdx()
                                    in substruct_index_list  # nitrogens in the indole substructure
                                ):  # atomic number 7 corresponds to Nitrogen

                                    atom.SetNumExplicitHs(
                                        1
                                    )  # gives 1 hydrogen to the N group
                            product_list.append(mol["rdkit_mol"])

        return product_list


# NOTE: tester
# NOTE: testing the replace algorithm
# molF = Chem.MolFromSmiles("Fc1c(F)c(Br)c2nsnc2c1Br")
# molCz = Chem.MolFromSmiles()
# print(molCz)
# substructF = Chem.MolFromSmiles("F")
# substructCz = Chem.MolFromSmiles("c1ccc2c(c1)[nH]c1ccccc12")
# print(molCz.HasSubstructMatch(substructCz))
# molCz2 = Chem.DeleteSubstructs(molCz, substructCz)
# # img = Draw.MolToImageFile(molA2, "new_mol2.png")

# NOTE: testing the atom index algorithm
# for atom in molCz.GetAtoms():
#     if atom.GetAtomicNum() == 7:
#         idx = atom.GetIdx()

# NOTE: testing the substruct index algorithm
# molA = Chem.MolFromSmiles("CS(=O)(=O)c1cncc(B(O)O)c1")
# substructA = Chem.MolFromSmiles("B(O)O")
# new_molA = Chem.DeleteSubstructs(molA, substructA)
# molB = Chem.MolFromSmiles("C[N+]12CC(=O)O[B-]1(c1cc(I)cc(C(F)(F)F)c1)OC(=O)C2")
# substructB = Chem.MolFromSmiles("I")
# molA = Chem.MolFromSmiles("OB(O)c1ccc2ncccc2c1")
# substructA = Chem.MolFromSmiles("B(O)O")
# new_molA = Chem.DeleteSubstructs(molA, substructA)
# molB = Chem.MolFromSmiles("C[N+]12CC(=O)O[B-]1(c1c(F)cccc1Br)OC(=O)C2")
# substructB = Chem.MolFromSmiles("Br")
# substruct_index_tuple = molA.GetSubstructMatch(substructA)
# print(substruct_index_tuple)
# if 0 in substruct_index_tuple:
#     substruct_index = 0
# else:
#     substruct_index = min(substruct_index_tuple) - 1
# new_mol = Chem.ReplaceSubstructs(
#     molB,
#     substructB,
#     new_molA,
#     replaceAll=True,
#     replacementConnectionPoint=substruct_index,
# )
# new_mol_smiles = Chem.MolToSmiles(new_mol[0])
# print(new_mol_smiles)
# img = Draw.MolToImageFile(new_mol[0], "new_mol.png")
# mol = Chem.MolFromSmarts("Brccccc-c1cc(-cccccBr)c2ccccc2c1")
# print(mol)
# img = Draw.MolToImageFile(mol, "smarts_mol.png")

# NOTE: testing removing NBoc group
# molA = Chem.MolFromSmarts(
#     "CC(C)(C)OC(=O)n1ccc2cc(-c3ccc4oc(-c5cc(-[nH]6c7ccccc7c7ccccc76)c(-c6cc7cc(-c8ccc9c(ccn9C(=O)OC(C)(C)C)c8)ccc7o6)cc5-[nH]5c6ccccc6c6ccccc65)cc4c3)ccc21"
# )
# print(molA)
# print(molA)
# substructA = Chem.MolFromSmiles("CC(C)(C)OC=O")
# new_molA = Chem.DeleteSubstructs(molA, substructA)
# atoms = new_molA.GetAtoms()
# for atom in atoms:
#     if atom.GetIsAromatic() == True and atom.GetAtomicNum() == 7:
#         atom.SetNumExplicitHs(1)

# print(Chem.MolToSmiles(new_molA))
# img = Draw.MolToImageFile(molA, "snar_mol.png")

# NOTE: Testing SNAr reaction
# 3 cases: replace all Fs, replace only aromatic Fs, replace specific F's
# NOTE: Experiment
from rdkit.Chem.Draw import rdMolDraw2D

d = rdMolDraw2D.MolDraw2DCairo(1000, 800)
d.drawOptions().addAtomIndices = True
product_list = []
mol_ArF = Chem.MolFromSmiles(
    "Cc1nc2c(-c3cccc(-c4cc(-n5c6ccccc6c6ccccc65)cc(-n5c6ccccc6c6ccccc65)c4)c3-n3c4ccccc4c4ccccc43)c3nsnc3c(-c3cccc(-c4cc(-n5c6ccccc6c6ccccc65)cc(-n5c6ccccc6c6ccccc65)c4)c3-n3c4ccccc4c4ccccc43)c2nc1C"
)
sym_list = list(Chem.rdmolfiles.CanonicalRankAtoms(mol_ArF, breakTies=False))
print(sym_list)
sym_dict = {}
for sym_class in sym_list:
    if str(sym_class) not in list(sym_dict.keys()):
        sym_dict[str(sym_class)] = 1
    elif str(sym_class) in list(sym_dict.keys()):
        sym_dict[str(sym_class)] += 1
print(sym_dict)
print(mol_ArF.GetAtomWithIdx(18).GetAtomicNum())

d.DrawMolecule(mol_ArF)
d.FinishDrawing()
d.WriteDrawingText("atom_annotation_5.png")
# mol_Cz = Chem.MolFromSmiles("c1ccc2c(c1)[nH]c1ccccc12")
# nH_grp = Chem.MolFromSmarts("[nH]")
# substitute_grp = mol_Cz
# for atom in mol_Cz.GetAtoms():
#     if atom.GetAtomicNum() == 7:  # atomic number 7 corresponds to Nitrogen
#         nH_index = atom.GetIdx()
# sub = Chem.MolFromSmiles("[Rn]")

# for atom in mol_ArF.GetAtoms():
#     atom.SetAtomMapNum(atom.GetIdx())
# img = Draw.MolToImageFile(mol_ArF, "snar_mol_indexed.png")
# atom_index_list = []
# for atom in mol_ArF.GetAtoms():
#     if atom.GetAtomicNum() == 9:
#         for neighbor in atom.GetNeighbors():
#             if neighbor.GetIsAromatic():
#                 atom_index_list.append(atom.GetIdx())

# substruct_paradifluorobenzene = Chem.MolFromSmiles("Fc1ccc(F)cc1")
# substruct_thiadiazole = Chem.MolFromSmiles("c2ccc1nsnc1c2")
# # print(mol_ArF.HasSubstructMatch(substruct_thiadiazole))
# # print(mol_ArF.HasSubstructMatch(substruct_paradifluorobenzene))
# if len(atom_index_list) == 2:
#     for index in atom_index_list:
#         atom = mol_ArF.GetAtomWithIdx(index)
#         atom.SetAtomicNum(
#             86
#         )  # change to a new structure that is not present in dataset
#         mol_ArF = Chem.ReplaceSubstructs(
#             mol_ArF,
#             sub,
#             substitute_grp,
#             replaceAll=True,
#             replacementConnectionPoint=nH_index,
#         )[0]
#     product_list.append(mol_ArF)
# elif mol_ArF.HasSubstructMatch(substruct_paradifluorobenzene):
#     # atom index pairings go around the structure.
#     # So, if there are 6 total F's. The 1st and 4th F's are paired.
#     divider = len(atom_index_list) // 2
#     # pairings of F groups
#     pairings = []
#     for i in range(int(len(atom_index_list) / 2)):
#         pairings.append([atom_index_list[i], atom_index_list[i + divider]])
#     for j in range(len(pairings)):
#         comb = list(combinations(pairings, j + 1))
#         for comb_pairings in comb:
#             temp_pairings = []
#             for pair in comb_pairings:
#                 temp_pairings.extend(pair)
#             # new copy of mol so original is not modified
#             new_mol_ArF = copy.copy(mol_ArF)
#             for atom_index in temp_pairings:
#                 new_mol_ArF.GetAtomWithIdx(atom_index).SetAtomicNum(86)
#             new_mol_ArF = Chem.ReplaceSubstructs(
#                 new_mol_ArF,
#                 sub,
#                 substitute_grp,
#                 replaceAll=True,
#                 replacementConnectionPoint=nH_index,
#             )[0]
#             for atom in new_mol_ArF.GetAtoms():
#                 if atom.GetAtomicNum() == 7:
#                     atom.SetNumExplicitHs(0)
#             product_list.append(new_mol_ArF)

# elif mol_ArF.HasSubstructMatch(substruct_thiadiazole):
#     pairings = []
#     first_index = 0
#     last_index = len(atom_index_list) - 1
#     for i in range(int(len(atom_index_list) / 2)):
#         pairings.append([atom_index_list[first_index], atom_index_list[last_index]])
#         first_index += 1
#         last_index -= 1
#     for j in range(len(pairings)):
#         comb = list(combinations(pairings, j + 1))
#         for comb_pairings in comb:
#             temp_pairings = []
#             for pair in comb_pairings:
#                 temp_pairings.extend(pair)
#             # new copy of mol so original is not modified
#             new_mol_ArF = copy.copy(mol_ArF)
#             for atom_index in temp_pairings:
#                 new_mol_ArF.GetAtomWithIdx(atom_index).SetAtomicNum(86)
#             new_mol_ArF = Chem.ReplaceSubstructs(
#                 new_mol_ArF,
#                 sub,
#                 substitute_grp,
#                 replaceAll=True,
#                 replacementConnectionPoint=nH_index,
#             )[0]
#             for atom in new_mol_ArF.GetAtoms():
#                 if atom.GetAtomicNum() == 7:
#                     atom.SetNumExplicitHs(0)
#             product_list.append(new_mol_ArF)

# for product in product_list:
#     print(Chem.MolToSmiles(product))
#     print(Chem.MolFromSmiles(Chem.MolToSmiles(product)))

# mol_ArF = Chem.ReplaceSubstructs(
#     mol_ArF, sub, substitute_grp, replaceAll=True, replacementConnectionPoint=nH_index,
# )[0]
# product_list.append(mol_ArF)
# mol = product_list[0]
# for atom in mol.GetAtoms():
#     if atom.GetAtomicNum() == 7:
#         atom.SetNumExplicitHs(0)

# img = Draw.MolToImageFile(mol, "snar_mol.png")


# NOTE: testing indole grp substructure for deBoc
# indole_grp = Chem.MolFromSmiles("c2ccc1[nH]ccc1c2")
# carbazole_grp = Chem.MolFromSmiles("c1ccc3c(c1)[nH]c2ccccc23")
# substruct_index_list = []
# substruct_index_list2 = []
# molagain = Chem.MolFromSmiles(
#     "Cc1nc2c(-c3c(-c4ccc5[nH]ccc5c4)cccc3-n3c4ccccc4c4ccccc43)c3nsnc3c(-c3c(-c4ccc5[nH]ccc5c4)cccc3-n3c4ccccc4c4ccccc43)c2nc1C"
# )
# substruct_index_tuples = molagain.GetSubstructMatches(indole_grp)
# substruct_index_tuples2 = molagain.GetSubstructMatches(carbazole_grp)
# for tuple in substruct_index_tuples:
#     substruct_index_list.extend(list(tuple))
# for tuple in substruct_index_tuples2:
#     substruct_index_list2.extend(list(tuple))
# print(substruct_index_list, substruct_index_list2)
# duplicate_index_list = list(set(substruct_index_list) & set(substruct_index_list2))
# result_list = [
#     item for item in substruct_index_list if item not in substruct_index_list2
# ]

# print(result_list)

# for atom in molagain.GetAtoms():
#     atom.SetAtomMapNum(atom.GetIdx())
# img = Draw.MolToImageFile(molagain, "deBoc_mol_index.png")

# find index of substructure deletion #NOTE: DOESN'T work because deleteSubstruct rearranges molecule
# if mol["rdkit_mol"].HasSubstructMatch(sub):
# substruct_index_list = []
# substruct_index_tuples = mol[
# "rdkit_mol"
# ].GetSubstructMatches(sub)
# print("tuples: ", substruct_index_tuples)
# for tuple in substruct_index_tuples:
# if 0 in tuple:
# substruct_index_list.append(0)
# else:
# substruct_index_list.append(min(tuple) - 1)
# print("list: ", substruct_index_list)
# print("before: ", Chem.MolToSmiles(mol["rdkit_mol"]))
# mol["rdkit_mol"] = Chem.DeleteSubstructs(
# mol["rdkit_mol"], sub
# )
# print("after: ", Chem.MolToSmiles(mol["rdkit_mol"]))

# NOTE: bmida has boronic acid substructure
# bmida = Chem.MolFromSmiles("B1OC(=O)CN(C)CC(=O)O1")
# substruct = Chem.MolFromSmiles("B(O)O")
# print(bmida.HasSubstructMatch(substruct))

# NOTE: BHA TESTING
# bha_mol = Chem.MolFromSmiles(
#     "Fc1cc(-c2cccnc2-c2ccc3[nH]ccc3c2)c(F)cc1-c1cccnc1-c1ccc2[nH]ccc2c1"
# )
# indole_grp = Chem.MolFromSmiles("c2ccc1[nH]ccc1c2")
# substitute_grp = Chem.MolFromSmiles("c1ccc2c(c1)ccn2c3cnccn3")
# new_bha_mol = Chem.ReplaceSubstructs(
#     bha_mol, indole_grp, substitute_grp, replaceAll=True
# )
# print(Chem.MolToSmiles(new_bha_mol[0]))
# img = Draw.MolToImageFile(new_bha_mol[0], "new_BHA_mol.png")
