from rxn4chemistry import RXN4ChemistryWrapper
import time

API_KEY = "apk-ccbb3f6f74af6962119bc9f7461e784a263b14fef8080c320610c016b69057164d7317ae96eff2deda922d13b59ccc330e924ab679d8c7bdada2610ab690867df7277c559da08fc3ee10d269016ee747"
PROJECT_ID = "60e8605799348f0001561731"

rxn4chemistry_wrapper = RXN4ChemistryWrapper(api_key=API_KEY)
rxn4chemistry_wrapper.set_project(PROJECT_ID)
time.sleep(5)
results = rxn4chemistry_wrapper.paragraph_to_actions(
    r"In accordance with modified literature procedures, [S6] a 250-mL flask equipped with a large olive-shaped magnetic stirring bar was charged with anhydrous tin chloride (SnCl2, 19.0 g, 100 mmol, 2.5 equiv) and concentrated hydrochloric acid (100 mL). The resulting solution was stirred at room temperature to ensure complete dissolution of SnCl2. In parallel, a 500-mL three-necked flask equipped with an overhead mechanical stirrer, a dropping funnel and an internal thermometer was charged with the solid grinded 4-iodoaniline (8.76 g, 40.0 mmol, 1.0 equiv) followed by hydrochloric acid (6 M, 100 mL). The resulting suspension was vigorously stirred for 20 min at –5 °C. A solution of sodium nitrite (NaNO2, 3.31 g, 48.0 mmol, 1.2 equiv) in water (20 mL) was then added dropwise to the reaction mixture over 30 min (keeping internal reaction temperature below –5 °C). The obtained diazonium salt solution was further stirred at –5 °C for 1 h, and the SnCl2 solution was subsequently added dropwise with vigorous stirring over 30 min (keeping internal reaction temperature below –5 °C). The resulting mixture was further stirred for 1 h at room temperature. The resulting off-white suspension was filtered over a fritted funnel, and the filter cake was washed with isopropanol followed by Et2O. The obtained hydrazinium salt was further dried under high vacuum (10–2 mbar) for 15 h (the crude product was protected from light with an aluminium foil). The title compound S1 (8.05 g, 29.8 mmol, 75%) was obtained as a light brown powder and was used for the next step without further purification"
)

for action in results["actions"]:
    print(action)
