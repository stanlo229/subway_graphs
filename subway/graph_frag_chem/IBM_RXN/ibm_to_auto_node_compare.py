import json
import pkg_resources
import rdkit
from rdkit import Chem
from rdkit.Chem import rdmolops, rdChemReactions


def transforms_bmida_check(rs_smiles, ibm_smiles):
    """Function that converts the bmida structure in both SMILES to be the same, 
    and then returns whether they are the same (True) or not (False)
    """
    original_structure = "C[N+]12CC(=O)O[B-]1OC(=O)C2"
    bmida_modified_structure = "B1OC(=O)CN(C)CC(=O)O1"
    original_mol = Chem.MolFromSmiles(original_structure)
    bmida_mod_mol = Chem.MolFromSmiles(bmida_modified_structure)

    rs_mol = Chem.MolFromSmiles(rs_smiles)
    ibm_mol = Chem.MolFromSmiles(ibm_smiles)

    if rs_mol.HasSubstructMatch(bmida_mod_mol):
        rs_mol = rdmolops.ReplaceSubstructs(
            rs_mol, bmida_mod_mol, bmida_mod_mol, replaceAll=True
        )[0]

    rs_smiles = Chem.MolToSmiles(rs_mol)

    return rs_smiles == ibm_smiles


def reaction_check(rs_reaction, ibm_reaction):
    """Function that converts Reaction SMILES into reaction and recovers the 
    products and reactants to determine if routescore data and ibm_rxn data is the same
    """
    # NOTE: SB and BS reactions will have same product, so reactants must be included as well
    rs_rxn = rdChemReactions.ReactionFromSmarts(rs_reaction)
    ibm_rxn = rdChemReactions.ReactionFromSmarts(ibm_reaction)
    rs_reactants = rs_rxn.GetReactants()
    rs_product = rs_rxn.GetProductTemplate(0)
    ibm_reactants = ibm_rxn.GetReactants()
    ibm_product = ibm_rxn.GetProductTemplate(0)

    reaction_bool = {"reactants": False, "product": False}

    # check if reactants are the same
    # NOTE: deboc reaction is different since the reagent is included in reactants
    reactant_bool = False
    for rs_rct in rs_reactants:
        for ibm_rct in ibm_reactants:
            if transforms_bmida_check(Chem.MolToSmiles(rs_rct), Chem.MolToSmiles(ibm_rct)):


    print(rs_rxn == ibm_rxn)


# transforms_bmida_check(
#     "C[N+]12CC(=O)O[B-]1(c1sccc1-c1ccc3ncccc3c1)OC(=O)C2",
#     "CN1CC(=O)OB(c2sccc2-c2ccc3ncccc3c2)OC(=O)C1",
# )

reaction_check(
    "OB(O)c1ccc2ncccc2c1.CN1CC(=O)OB(c2sccc2Br)OC(=O)C1>CC(C)C1=CC(=C(C(=C1)C(C)C)C2=CC(=CC=C2)P(C3CCCCC3)C4CCCCC4)C(C)C.C1=CC=C([C-]=C1)C2=CC=CC=C2N.Cl[Pd+].[O-]P(=O)([O-])[O-].[K+].[K+].[K+]>CN1CC(=O)OB(c2sccc2-c2ccc3ncccc3c2)OC(=O)C1",
    "OB(O)c1ccc2ncccc2c1.CN1CC(=O)OB(c2sccc2Br)OC(=O)C1>CC(C)C1=CC(=C(C(=C1)C(C)C)C2=CC(=CC=C2)P(C3CCCCC3)C4CCCCC4)C(C)C.C1=CC=C([C-]=C1)C2=CC=CC=C2N.Cl[Pd+].[O-]P(=O)([O-])[O-].[K+].[K+].[K+]>CN1CC(=O)OB(c2sccc2-c2ccc3ncccc3c2)OC(=O)C1",
)

