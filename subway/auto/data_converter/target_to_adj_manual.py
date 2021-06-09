import pickle
import csv
import json

import pkg_resources

"""
Format: [node_from, node_to]
"""

JSON_PATH = pkg_resources.resource_filename("subway", "data/auto_nodes.json")
PKL_PATH = pkg_resources.resource_filename("subway", "data/adj_list.pkl")


def clear_pkl(pkl_path):
    empty = []
    file = open(pkl_path, "wb")
    pickle.dump(empty, file)
    file.close()


def read_pkl(pkl_path):
    file = open(pkl_path, "rb")
    data = pickle.load(file)
    file.close()
    print(len(data))


def pkl_to_txt(pkl_path):
    file = open(pkl_path, "rb")
    data = pickle.load(file)
    file.close()
    f = open("adj_list.txt", "w+")
    f.write(str(data))
    f.close()


def manual_to_adj(pkl_path, json_path):
    # open file
    file = open(pkl_path, "rb")
    data = pickle.load(file)
    file.close()

    # check if file is empty
    if len(data) != 0:
        adj_list = data
    else:
        adj_list = []
    file = open(pkl_path, "wb")

    manual_data = [
        [7750, 32],
        [3872, 7750],
        [3873, 7750],
        [3874, 7750],
        [3868, 7749],
        [7749, 3872],
        [3869, 7749],
        [3870, 7749],
        [3871, 7749],
        [7752, 3884],
        [3878, 7752],
        [3879, 7752],
        [3880, 7752],
        [42, 7752],
        [7751, 42],
        [3875, 7751],
        [3876, 7751],
        [3877, 7751],
        [7753, 40],
        [3881, 7753],
        [3882, 7753],
        [3883, 7753],
    ]
    adj_list.extend(manual_data)

    # dump information to that file
    pickle.dump(adj_list, file)
    # close the file
    file.close()


# manual_to_adj(PKL_PATH, JSON_PATH)
# pkl_to_txt(PKL_PATH)
# clear_pkl(PKL_PATH)
# read_pkl(PKL_PATH)
