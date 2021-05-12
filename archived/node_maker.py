import json


path_to_molecule = 'molecule_node_basic.json'
with open(path_to_molecule) as f:
    mol_data = json.load(f)

path_to_reaction = 'reaction_node_basic.json'
with open(path_to_reaction) as f2:
    rxn_data = json.load(f2)

class MoleculeNode:

    def __init__(self, SMILES: str):
        '''
        parent_r_nodes: list of next possible reaction nodes
        child_r_nodes: list of previous possible reaction nodes
        NOTE: 
        starting material: no child node 
        product: no parent node
        '''
        self.parent_r_nodes = []
        self.child_r_nodes = []
        self.data = list(filter(lambda x:x["SMILES"]== SMILES, mol_data)) #not the fastest approach, can optimize speed here
    
    def insert(self, parent_r_nodes, child_r_nodes):
        # needs editing to update multiple parent/child nodes like an empty list with appending
        self.parent_r_nodes.append(parent_r_nodes)
        self.child_r_nodes.append(child_r_nodes)


class ReactionNode:

    def __init__(self, rxn_SMILES: str):
        '''
        parent_m_nodes: list of product molecule nodes (never null)
        child_m_nodes: list of molecule nodes used (never null)
        '''
        self.parent_m_nodes = []
        self.child_m_nodes = []
        #searches for specific reaction smiles in .json file
        self.data = list(filter(lambda x:x["rxn_SMILES"]== rxn_SMILES, rxn_data)) #not the fastest approach, can optimize speed here
    
    def insert(self, parent_m_nodes, child_m_nodes):
        self.parent_m_nodes.append(parent_m_nodes)
        self.child_m_nodes.append(child_m_nodes)

print(ReactionNode("Brc1cc(C=O)c(O)cc1.BrC(Br)(Br)Br>P(c1ccccc1)(c2ccccc2)c3ccccc3.CCN(CC)CC>Brc1cc(/C=C(Br)/Br)c(O)cc1").data[0]["rxn_SMILES"])