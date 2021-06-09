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
    [1, 11],
    [11, 10],
    [10, 9],
    [12, 10],
    [14, 10],
    [15, 10],
    [14, 13],
    [15, 13],
    [18, 13],
    [19, 13],
    [13, 12],
]
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
    [1, 11],
    [11, 10],
    [10, 9],
    [12, 10],
    [14, 10],
    [15, 10],
    [14, 13],
    [15, 13],
    [18, 13],
    [19, 13],
    [13, 12],
]
"""


def make_pkl(pkl_path, data):
    file = open(pkl_path, "wb")
    pickle.dump(data, file)
    file.close()


make_pkl(pkl_path, adj_list)
