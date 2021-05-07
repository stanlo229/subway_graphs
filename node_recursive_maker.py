import routescore
import pickle
import json
from rdkit import Chem
from rdkit.Chem import rdChemReactions
import timeit

path_to_molecule = 'molecule_node_basic.json'
with open(path_to_molecule) as f:
    mol_data = json.load(f)

path_to_reaction = 'reaction_node_basic.json'
with open(path_to_reaction) as f2:
    rxn_data = json.load(f2)

class Node:
    """"Class for creating the tree which are all of the routes from a desired product"""

    def __init__(self, SMILES: str, rxn_binary: int): #rxn binary string
        self.rxn_binary = rxn_binary
        if rxn_binary == 0: #molecular node
            self.parent = []
            self.child = []
            self.data = list(filter(lambda x:x["SMILES"] == SMILES, mol_data)) #not the fastest approach, can optimize speed here
            for r in rxn_data:
                rxn = rdChemReactions.ReactionFromSmarts(r["rxn_SMILES"])
                product = rxn.GetProductTemplate(0)
                data_SMILES = Chem.MolToSmiles(product)
                a = Chem.CanonSmiles(data_SMILES)
                b = Chem.CanonSmiles(SMILES)
                if a == b:
                    self.child.append(Node(r["rxn_SMILES"], 1))
                '''reactant = rxn.GetReactants()
                for i in reactant:
                    r_SMILES = Chem.MolToSmiles(i)
                    y = Chem.CanonSmiles(r_SMILES)
                    if b == y:
                        self.parent.append(Node(r["rxn_SMILES"], 1))'''
        elif rxn_binary == 1:
            #searches for specific reaction smiles in .json file
            self.data = list(filter(lambda x:x["rxn_SMILES"]== SMILES, rxn_data)) #not the fastest approach, can optimize speed here
            rxn_SMILES = self.data[0]["rxn_SMILES"]
            rxn = rdChemReactions.ReactionFromSmarts(rxn_SMILES)
            #product = rxn.GetProductTemplate(0)
            #product_SMILES = Chem.MolToSmiles(product)
            #self.parent = []
            self.child = []
            for m in mol_data:
                x = Chem.CanonSmiles(m["SMILES"])
                reactant = rxn.GetReactants()
                agents = rxn.GetAgents()
                for i in reactant:
                    r_SMILES = Chem.MolToSmiles(i)
                    y = Chem.CanonSmiles(r_SMILES)
                    if x == y:
                        self.child.append(Node(m["SMILES"], 0))
                for j in agents:
                    a_SMILES = Chem.MolToSmiles(j)
                    z = Chem.CanonSmiles(a_SMILES)
                    if x == z:
                        self.child.append(Node(m["SMILES"], 0))


def RouteSearch(root_node: Node):
    scores = {}
    if root_node.rxn_binary == 0 and len(root_node.child) != 0:
        for route in root_node.child:
            route_name = route.data["rxn_SMILES"]
            scores[route_name] += RouteSearch(route)
    elif root_node.rxn_binary == 1 and len(root_node.child) != 0: #checks if molecule and is not last molecule
        cH = 52.97
        cM = 17.13
        tH = root_node.data["tH"]
        tM = root_node.data["tM"]
        yld = root_node.data["yield"]
        mw_stepscore = 0
        cost_stepscore = 0
        for mol in root_node.child:
            quantity = mol.data["quantity"]
            cost = mol.data["cost"]
            mw = mol.data["g/mol"]
            mw_stepscore += mw*quantity
            cost_stepscore += cost*quantity
        stepscore = tH * ((tH*cH) + cost_stepscore) * (mw_stepscore)
        return stepscore
    elif root_node.rxn_binary == 0 and len(root_node.child) == 0:
        return scores
            

trial = Node("Brc1ccc2oc(Br)cc2c1", 0)
print(trial.child[0].child)
