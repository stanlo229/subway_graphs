from rxn4chemistry import RXN4ChemistryWrapper
from rdkit import Chem
from rdkit.Chem import rdChemReactions
import time

API_KEY = "apk-ccbb3f6f74af6962119bc9f7461e784a263b14fef8080c320610c016b69057164d7317ae96eff2deda922d13b59ccc330e924ab679d8c7bdada2610ab690867df7277c559da08fc3ee10d269016ee747"
PROJECT_ID = "60e8605799348f0001561731"


# setup project
# rxn4chemistry_wrapper = RXN4ChemistryWrapper(api_key=API_KEY)
# rxn4chemistry_wrapper.create_project("subway_graph_batch")
# print(f"Identifier for the project {rxn4chemistry_wrapper.project_id}")

# rxn4chemistry_wrapper.set_project("60e8b0f399348f000156e858")

# reactants_list = ["BrBr.c1ccc2cc3ccccc3cc2c1", "Cl.c1ccc2cc3ccccc3cc2c1"] * 4

# response = rxn4chemistry_wrapper.predict_reaction_batch(precursors_list=reactants_list)
# print(response)
# time.sleep(5)
# results = rxn4chemistry_wrapper.get_predict_reaction_batch_results(response["task_id"])
# index = 0
# reaction_list = []
# while index < len(reactants_list):
#     reaction_list.append(results["predictions"][index]["smiles"])
#     index += 1

# print(reaction_list)

# results = rxn4chemistry_wrapper.get_predict_reaction_results(response["prediction_id"])
# smile = results["response"]["payload"]["attempts"][0]["smiles"]
# smile_modify = "C[N+]12CC(=O)O[B-]1(c1sccc1Br)OC(=O)C2.OB(O)c1ccc2ncccc2c1>>C[N+]12CC(=O)OB(c3sccc3-c3ccc4ncccc4c3)OC(=O)C1"
# index = 0
# product_next = 0
# while index < len(smile_modify):
#     curr_char = smile_modify[index]
#     if curr_char == ">":
#         product_next += 1
#     if product_next == 2:
#         if curr_char == "B":
#             next_char = smile_modify[index + 1]
#             if next_char != "r":
#                 smile_modify = smile_modify[:index] + "[" + smile_modify[index:]
#                 smile_modify = (
#                     smile_modify[: index + 2] + "-]1" + smile_modify[index + 2 :]
#                 )
#                 index += 4
#     index += 1
# print(smile_modify)

# rxn = rdChemReactions.ReactionFromSmarts(
#     "C[N+]12CC(=O)O[B-]1(c1sccc1Br)OC(=O)C2.OB(O)c1ccc2ncccc2c1>>C[N+]12CC(=O)OB(c3sccc3-c3ccc4ncccc4c3)OC(=O)C1"
# )
# print(rxn.GetProductTemplate(0))


class IBM_RXN:
    """ Class that contains functions to utilize the rxn4chemistry by IBM
    """

    def __init__(self, api_key, project_id):
        self.api_key = api_key
        self.project_id = project_id
        self.wrapper = RXN4ChemistryWrapper(api_key=api_key)
        self.wrapper.set_project(project_id)

    def check_bmida(self, reaction_smiles):
        """ Function that checks whether the product has a bmida substructure, and then returns product.
        If bmida is present, it modifies the string to have correct smiles syntax.
        Parameters
        ----------
        reaction_smiles: reaction SMILES

        Return
        ------
        product_smiles: smiles of product molecule
        """
        index = 0
        product_next = 0
        bmida = False
        # modify string to change B to [B-]2 if contains bmida substructure
        while index < len(reaction_smiles):
            curr_char = reaction_smiles[index]
            if curr_char == ">":
                product_next += 1
            # check for product
            if product_next == 2:
                # check if bmida structure exists in product ([N+])
                # definitely a better way (this is just a hotfix)
                if curr_char == "N":
                    next_char = reaction_smiles[index + 1]
                    if next_char == "+":
                        bmida = True
            index += 1
        # modify string
        index = 0
        product_next = 0
        while index < len(reaction_smiles):
            curr_char = reaction_smiles[index]
            if curr_char == ">":
                product_next += 1
            # check for product
            if product_next == 2:
                if bmida:
                    if curr_char == "B":
                        next_char = reaction_smiles[index + 1]
                        if next_char != "r":
                            reaction_smiles = (
                                reaction_smiles[:index] + "[" + reaction_smiles[index:]
                            )
                            reaction_smiles = (
                                reaction_smiles[: index + 2]
                                + "-]2"
                                + reaction_smiles[index + 2 :]
                            )
                            index += 4
            index += 1
        # get product smile string
        index = 0
        product_next = 0
        product_smile = ""
        while index < len(reaction_smiles):
            curr_char = reaction_smiles[index]
            # if 2 ">" are found, add char
            if product_next == 2:
                product_smile += curr_char
            # check if product
            if curr_char == ">":
                product_next += 1
            index += 1

        return product_smile

    def get_product(self, reactant_smiles):
        """Function that takes in same input as graph_frag_chem, and uses IBM RXN to generate the products
        Parameters
        ----------
        reactant_smiles: list of reactant smiles

        Return
        ------
        product_list: list of product rdkit molecules
        """
        rxn_reactants_smiles = ""
        index = 0
        while index < len(reactant_smiles):
            rxn_reactants_smiles += reactant_smiles[index]
            index += 1
            if index != len(reactant_smiles):
                rxn_reactants_smiles += "."
        response = self.wrapper.predict_reaction(rxn_reactants_smiles)
        # delay for 10 seconds because IBM RXN servers have a 5calls/min rule
        time.sleep(20)
        results = self.wrapper.get_predict_reaction_results(response["prediction_id"])
        smile = results["response"]["payload"]["attempts"][0]["smiles"]
        product_smile = self.check_bmida(smile)
        # handles errors when IBM RXN does not produce a viable product
        product_list = []
        product_mol = Chem.MolFromSmiles(product_smile)
        if product_mol == None:
            print("error: ", product_smile)
            product_list = product_smile
        else:
            product_list.append(product_smile)
        # delay for 10 seconds because IBM RXN servers have a 5calls/min rule
        time.sleep(5)

        return product_list

    def get_product_batch(self, reactants_list):
        """Function that takes in list of reactants, and uses IBM RXN to generate the products
        Parameters
        ----------
        reactant_list: list of reactant smile pairings

        Return
        ------
        product_list: list of product rdkit molecules
        """
        product_list = []
        ibm_reaction_list = []
        index = 0

        # call a batch of reactions
        response = self.wrapper.predict_reaction_batch(reactants_list)
        # delay for 20 seconds because IBM RXN servers have a 5calls/min rule and wait for reaction to finish
        time.sleep(20)
        results = self.wrapper.get_predict_reaction_batch_results(response["task_id"])
        # get product_smile from reaction_smile of each reaction and add to product_list
        while index < len(reactants_list):
            # add time buffer for predictions that take longer
            try:
                reaction_smile = results["predictions"][index]["smiles"]
            except:
                print("Reaction taking longer than usual")
                time.sleep(20)
                results = self.wrapper.get_predict_reaction_batch_results(
                    response["task_id"]
                )
                print(results)
                reaction_smile = results["predictions"][index]["smiles"]
            # NOTE: IBM rxn4chemistry predicts multiple products with different confidence levels.
            # So, it's possible to get more of them.
            # look for product smiles by SMARTS reaction
            try:
                rxn = rdChemReactions.ReactionFromSmarts(reaction_smile)
            except:
                product_smiles = None
            else:
                try:
                    product_mol = rxn.GetProductTemplate(0)  # only 1 product
                except:
                    product_smiles = None
                else:
                    product_smiles = Chem.MolToSmiles(product_mol)
                    product_mol_smiles = Chem.MolFromSmiles(product_smiles)
                    if product_mol_smiles == None:
                        product_smiles = None
            # NOTE: HAS TO BE SMILES OUTPUT because rdkit screws
            # up the smile when transferred over to other files
            product_list.append(product_smiles)
            ibm_reaction_list.append(reaction_smile)
            index += 1
        # delay for 12 seconds because IBM RXN servers have a 5calls/min rule
        time.sleep(5)
        return product_list, ibm_reaction_list


# ws1 = "OB(O)c1ccc2ncccc2c1"
# ws2 = "CC(C)(C)OC(=O)n1ccc2cc(B3OC(C)(C)C(C)(C)O3)ccc21"
ws3 = "CC(C)C1=CC(=C(C(=C1)C(C)C)C2=CC(=CC=C2)P(C3CCCCC3)C4CCCCC4)C(C)C~C1=CC=C([C-]=C1)C2=CC=CC=C2N~Cl[Pd+]"
ws4 = "[O-]P(=O)([O-])[O-]~[K+]~[K+]~[K+]"
# ws5 = "CN1CC(=O)OB(c2cccc(Br)c2F)OC(=O)C1"
ps1 = "O=C1CNCC(=O)O1"
ps2 = "Fc1c(F)c(Br)c2nsnc2c1Br"
# ps3 = "Clc1ccc(Cl)c(-c2ccc(-c3ccc(-c4ccc(-c5cc(Cl)ccc5Cl)cc4)nc3)cc2)c1"
# ps4 = "Brc1ccc(Br)nc1"

product_list, ibm_reaction_list = IBM_RXN(API_KEY, PROJECT_ID).get_product_batch(
    [ps1 + "." + ps2 + "." + ws3 + "." + ws4]
)

print(product_list, ibm_reaction_list)

