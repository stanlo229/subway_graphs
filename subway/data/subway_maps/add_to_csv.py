import pandas as pd
import pkg_resources

FULL_PROPS_PATH = pkg_resources.resource_filename(
    "subway", "data/subway_maps/full_props.csv"
)
SB_PATH = pkg_resources.resource_filename(
    "subway", "data/subway_maps/S-B_man_routes.csv"
)


def add_to_csv(full_props_path, s_b_path):
    full_props_df = pd.read_csv(full_props_path)
    s_b_df = pd.read_csv(s_b_path)
    index = 0
    while index < len(s_b_df):
        smiles = s_b_df["pentamer"][index]
        routescore = s_b_df["RouteScore"][index]
        full_props_df = full_props_df.append(
            {"smiles": smiles, "route_score": routescore}, ignore_index=True
        )
        index += 1
    full_props_df.to_csv(full_props_path)


add_to_csv(FULL_PROPS_PATH, SB_PATH)
