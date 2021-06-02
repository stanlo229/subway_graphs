import pickle
import pkg_resources
import itertools

PKL_PATH = pkg_resources.resource_filename("subway", "data/adj_list.pkl")


def remove_dup_pkl(pkl_path):
    file = open(pkl_path, "rb")
    data = pickle.load(file)
    file.close()
    data.sort()
    data = list(data for data, _ in itertools.groupby(data))

    file = open(pkl_path, "wb")
    # dump information to that file
    pickle.dump(data, file)
    # close the file
    file.close()


remove_dup_pkl(PKL_PATH)
