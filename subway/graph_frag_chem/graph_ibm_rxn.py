from rxn4chemistry import RXN4ChemistryWrapper

API_KEY = "apk-ccbb3f6f74af6962119bc9f7461e784a263b14fef8080c320610c016b69057164d7317ae96eff2deda922d13b59ccc330e924ab679d8c7bdada2610ab690867df7277c559da08fc3ee10d269016ee747"
PROJECT_ID = "60d3575e99348f00013922c7"

# setup project
rxn4chemistry_wrapper = RXN4ChemistryWrapper(api_key=API_KEY)
# rxn4chemistry_wrapper.create_project("subway_graph")

rxn4chemistry_wrapper.set_project("60d3575e99348f00013922c7")

response = rxn4chemistry_wrapper.predict_reaction(
    "Fc1cc(-c2c(F)cccc2-c2ccc3c(ccn3-c3cnccn3)c2)c(F)cc1-c1c(F)cccc1-c1ccc2c(ccn2-c2cnccn2)c1.c1ccc2c(c1)[nH]c1ccccc12"
)

results = rxn4chemistry_wrapper.get_predict_reaction_results(response["prediction_id"])
print(results["response"]["payload"]["attempts"][0]["smiles"])


class IBM_RXN:
    """ Class that contains functions to utilize the rxn4chemistry by IBM
    """

    def __init__(self, api_key, project_id):
        self.api_key = api_key
        self.project_id = project_id
        self.wrapper = RXN4ChemistryWrapper(api_key=api_key)
        self.wrapper.set_project(project_id)
    
    def run(self):
        # batch run?
        # per rxn run?

