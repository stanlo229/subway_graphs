import pkg_resources
import pickle

pkl_path = pkg_resources.resource_filename("subway", "basic/adj_list_basic.pkl")

"""
Format: [[node_from, node_to], ...]
"""
adj_list = [
    [2, 0],
    [3, 0],
    [4, 0],
    [5, 0],
    [0, 6],
    [6, 1],
    [7, 1],
    [8, 1],
    [1, 9],
    [10, 6],
    [11, 10],
    [12, 10],
    [13, 9],
    [14, 13],
    [15, 13],
    [4, 13],
    [16, 11],
    [17, 12],
    [18, 16],
    [19, 17],
]


def make_pkl(pkl_path, data):
    file = open(pkl_path, "wb")
    pickle.dump(data, file)
    file.close()


make_pkl(pkl_path, adj_list)
