import csv
import pkg_resources
import pandas as pd
import math

FULL_PROPS_CSV_PATH = pkg_resources.resource_filename(
    "subway", "data/subway_maps/full_props.csv"
)

RESULTS_CSV_PATH = pkg_resources.resource_filename("subway", "data/results.csv")


def accuracy_checker(full_props_csv, results_csv):
    full_props_df = pd.read_csv(full_props_csv)
    results_df = pd.read_csv(results_csv)

    results_index = 0
    matches = 0
    total = 0
    for pdt_smile in results_df["product_SMILES"]:
        full_props_index = 0
        for smile in full_props_df["smiles"]:
            if pdt_smile == smile:
                if round(results_df["routescore"][results_index], -1) == round(
                    full_props_df["route_score"][full_props_index], -1
                ):
                    matches += 1
                elif (
                    results_df["routescore"][results_index]
                    != full_props_df["route_score"][full_props_index]
                ):
                    print(full_props_index)
            full_props_index += 1
        total += 1
        results_index += 1

    print(matches / total)


accuracy_checker(FULL_PROPS_CSV_PATH, RESULTS_CSV_PATH)
