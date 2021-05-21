"""
mol_to_rxn_adj_matrix = [[1,0],[1,0],[1,0],[1,0],[1,1],[0,1],[0,1],[0,1]]

rxn_to_mol_adj_matrix = [[1,1,1,1,1,0,0,0],[0,0,0,0,1,1,1,1]]

mol_to_rxn_adj_list = [[0],[0],[0],[0],[0,1],[1],[1],[1]]

rxn_to_mol_adj_list = [[0,1,2,3,4],[4,5,6,7]]

rxn_to_mol_from = [0,0,0,0,0,1,1,1,1]
rxn_to_mol_to   = [0,1,2,3,4,4,5,6,7]
"""

# case 1
mol_to_rxn_from = [0, 1, 2, 3, 4, 4, 5, 6, 7]
mol_to_rxn_to = [0, 0, 0, 0, 0, 1, 1, 1, 1]

"""
# case 2
mol_to_rxn_from = [0,1,2,3,4,4,4,5,6,7,8,9]
mol_to_rxn_to   = [0,0,0,0,0,1,2,1,1,1,2,2]
"""
"""
# case 3
mol_to_rxn_from = [0,1,2,3,4,4,4,5,6,7,7,8,9,10,11]
mol_to_rxn_to   = [0,0,0,0,0,1,2,1,1,1,3,2,2,3,3]
"""
"""
# case 4
mol_to_rxn_from = [0,1,2,3,4,4,4,5,6,7,7,7,8,9,10,11,11,12]
mol_to_rxn_to   = [0,0,0,0,0,1,2,1,1,1,3,4,2,2,3,3,4,4]
"""
"""
# case 5
mol_to_rxn_from = [0,1,2,2,3,4,4,4,4,5,5,6,6,7,7,7,8,9,10,11,11,12]
mol_to_rxn_to   = [0,0,0,1,0,0,1,2,4,1,2,1,2,1,3,4,2,2,3,3,4,4]
"""

# case 6
# mol_to_rxn_from = [0, 1, 2, 2, 3, 4, 4, 4, 4, 5, 6, 7, 7, 8, 9, 10, 11, 11, 12, 12]
# mol_to_rxn_to = [0, 0, 0, 1, 0, 0, 1, 2, 4, 1, 1, 1, 3, 2, 2, 3, 3, 4, 4, 5]

# case 7
# mol_to_rxn_from = [0,1,2,3,3,4,4,4,4,5,6,7,8,9]
# mol_to_rxn_to   = [0,0,0,0,3,0,1,2,3,1,1,1,2,2]

# case 8
# mol_to_rxn_from = [0,1,2,2,3,4,4,4,4,4,5,5,6,6,7,7,7,8,9,10,11,11,12]
# mol_to_rxn_to   = [0,0,0,1,0,0,1,2,4,5,1,2,1,2,1,3,4,2,2,3,3,4,4]

# NOTE: next step is to add more reactions and molecules to the mix
# and see if it works with creating part of the graph
