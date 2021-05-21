import pickle

data = pickle.load(
    open(
        r"C:\Users\Stanley Lo\Documents\Summer 2021 NSERC\subway_graphs\subway\data\auto_adj_list.pkl",
        "rb",
    )
)
output = open("write.txt", "w")
output.write(str(data))
output.flush()
output.close()
