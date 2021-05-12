'''
mol_to_rxn_adj_matrix = [[1,0],[1,0],[1,0],[1,0],[1,1],[0,1],[0,1],[0,1]]

rxn_to_mol_adj_matrix = [[1,1,1,1,1,0,0,0],[0,0,0,0,1,1,1,1]]

mol_to_rxn_adj_list = [[0],[0],[0],[0],[0,1],[1],[1],[1]]

rxn_to_mol_adj_list = [[0,1,2,3,4],[4,5,6,7]]

rxn_to_mol_from = [0,0,0,0,0,1,1,1,1]
rxn_to_mol_to   = [0,1,2,3,4,4,5,6,7]
'''

mol_to_rxn_from = [0,1,2,3,4,4,5,6,7]
mol_to_rxn_to   = [0,0,0,0,0,1,1,1,1]
#NOTE: next step is to add more reactions and molecules to the mix 
# and see if it works with creating part of the graph