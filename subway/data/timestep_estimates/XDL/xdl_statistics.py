import pandas as pd
import pkg_resources

TIMESTEP_CSV = pkg_resources.resource_filename(
    "subway", "data/timestep_estimates/XDL/Timestep Estimates Manual Annotation.csv"
)
ERROR_CSV = pkg_resources.resource_filename(
    "subway", "data/timestep_estimates/XDL/xdl_errors.csv"
)
GLOBAL_STATS_CSV = pkg_resources.resource_filename(
    "subway", "data/timestep_estimates/XDL/global_stats.csv"
)


class xdl_stats:
    """
    Class that contains functions for summarizing statistics of XDL outputs, and creates csv list of errors (w/ procedure)
    """

    def __init__(self, timestep_csv):
        self.data = pd.read_csv(timestep_csv)

    def create_error_csv(self, error_csv):
        """
        Function that compiles all "Error: " msgs
        """
        error_df = pd.DataFrame()
        error_df["Input(Procedure Step)"] = " "
        error_df["Error Message"] = " "
        index = 0
        for i in range(self.data.shape[0]):  # iterate through rows
            for j in range(self.data.shape[1]):  # iterate through columns
                value = self.data.iloc[i, j]
                if isinstance(value, str):
                    if value[0:5] == "Error":
                        error_df.at[index, "Input(Procedure Step)"] = self.data.iloc[
                            i, j - 2
                        ]
                        error_df.at[index, "Error Message"] = self.data.iloc[i, j]
                        index += 1
        error_df.to_csv(error_csv)

    def global_stats(self, global_stats_csv):
        """
        Function that summarizes all the Action vs. XDL data.

        Global Statistics
        -------------------
        Accuracy = Matching Actions / Total Actions
            Total Actions will contain extras from "Action" and "XDL"
        Error % = Number of Errors / Total Number of Steps
        """
        # find Action and XDL columns
        action_idx = []
        xdl_idx = []
        for j in range(self.data.shape[1]):  # iterate through columns
            value = self.data.iloc[0, j]
            if value == "Action":
                action_idx.append(j)
            if value == "XDL":
                xdl_idx.append(j)

        idx = 0
        ord_match = 0
        total_actions = 0
        error = 0
        while idx < len(action_idx):
            a_idx = action_idx[idx]
            x_idx = xdl_idx[idx]
            for i in range(self.data.shape[0]):  # iterate through rows
                # create list of actions
                a_value = self.data.iloc[i, a_idx]
                x_value = self.data.iloc[i, x_idx]
                if isinstance(a_value, float) or isinstance(x_value, float):
                    continue
                # count number of errors
                if x_value[0:5] == "Error":
                    error += 1
                a_list = a_value.split(",")
                x_list = x_value.split(",")
                print(a_list, x_list)

                # add max number of actions to total_actions (includes extra actions)
                total_actions += max(len(a_list), len(x_list))

                # ordered matches
                if len(a_list) < len(x_list):
                    for order in range(len(a_list)):
                        if a_list[order] == x_list[order]:
                            ord_match += 1
                else:
                    for order in range(len(x_list)):
                        if a_list[order] == x_list[order]:
                            ord_match += 1

            idx += 1
        print(
            "Matches:",
            ord_match,
            "Total # of Actions:",
            total_actions,
            "Accuracy(%):",
            (ord_match * 100 / total_actions),
            "Errors:",
            error,
            "Accuracy_w/o_errors(%):",
            (ord_match * 100 / (total_actions - error)),
        )


stats = xdl_stats(TIMESTEP_CSV)
# stats.create_error_csv(ERROR_CSV)
stats.global_stats(GLOBAL_STATS_CSV)
