import pkg_resources
import pickle

pkl_path = pkg_resources.resource_filename("subway", "basic/adj_list_basic.pkl")

"""
Format: [[node_from, node_to], ...]
"""
adj_list = [[2, 0], [3, 0], [4, 0], [5, 0], [0, 6], [6, 1], [7, 1], [8, 1], [1, 9]]

for (i, j) in adj_list:
    print(i, j)


def make_pkl(pkl_path, data):
    file = open(pkl_path, "wb")
    pickle.dump(data, file)
    file.close()


make_pkl(pkl_path, adj_list)
