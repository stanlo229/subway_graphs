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
        [12821, 32],
        [3872, 12821],
        [3873, 12821],
        [3874, 12821],
        [3868, 12820],
        [12820, 3872],
        [3869, 12820],
        [3870, 12820],
        [3871, 12820],
        [12823, 42],
        [3878, 12823],
        [3879, 12823],
        [3880, 12823],
        [3884, 12823],
        [12822, 3884],
        [3875, 12822],
        [3876, 12822],
        [3877, 12822],
        [12824, 40],
        [3881, 12824],
        [3872, 12824],
        [3883, 12824],
    ]
    adj_list.extend(manual_data)

    # dump information to that file
    pickle.dump(adj_list, file)
    # close the file
    file.close()


# manual_to_adj(PKL_PATH, JSON_PATH)
# pkl_to_txt(PKL_PATH)
# clear_pkl(PKL_PATH)
read_pkl(PKL_PATH)
