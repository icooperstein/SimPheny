
import pandas as pd
from pyhpo import Ontology
_=Ontology()
phenotype_data = pd.read_csv('input.csv')

''' Create dictionary of all patients' term lists'''

## "Curated_Terms" column should have list of HPO terms in the form 'Term1; Term2; Term3'
## "ID" column should have patient identifier
patient_list = list(set(phenotype_data['ID']))

patient_terms_dict = {}
for patient in patient_list:
    patient_terms_dict[patient] = phenotype_data.loc[phenotype_data['ID']==patient]['Curated_Terms'].item().split('; ')


'''Perform phenotype category assignments'''
top_tier = ['Abnormal cellular phenotype', 'Abnormality of blood and blood-forming tissues', 'Abnormality of head or neck', 'Abnormality of limbs', 'Abnormality of metabolism/homeostasis', 'Abnormality of prenatal development or birth', 'Abnormality of the breast', 'Abnormality of the cardiovascular system', 'Abnormality of the digestive system', 'Abnormality of the ear', 'Abnormality of the endocrine system', 'Abnormality of the eye', 'Abnormality of the genitourinary system', 'Abnormality of the immune system', 'Abnormality of the integument', 'Abnormality of the musculoskeletal system', 'Abnormality of the nervous system', 'Abnormality of the respiratory system', 'Abnormality of the thoracic cavity', 'Abnormality of the voice', 'Constitutional symptom', 'Growth abnormality', 'Neoplasm']
print(len(top_tier))
winners = []
max_vals = []
all_winners=[]
for patient in patient_list:
    parentnames = []
    patient_terms_data =[]
    for term in patient_terms_dict[patient]:
        term_ic = Ontology.get_hpo_object(term).information_content['omim']
        if term in top_tier:
            parentnames.append(term)
            patient_terms_data.append([term, term_ic])

        else:
            parents = Ontology.get_hpo_object(term).all_parents
            parent_names = [parent.name for parent in parents]
            overlap = list(set(parent_names) & set(top_tier))
            for parent_terms in overlap:
                patient_terms_data.append([parent_terms, term_ic])
        
    results_table = pd.DataFrame(patient_terms_data, columns=['Term', 'Sum']).groupby('Term').sum().sort_values(by='Sum', ascending=False)
    winner = list(results_table.head(1).index)[0]
    winners.append(winner)
    max_val = results_table.head(1)['Sum'].item()
    max_vals.append(max_val)
    all_winners.append([patient, winner, max_val])
    vals = results_table.values.tolist()
    if vals.count(max_val) > 1 :
        print('error')
    result_terms = results_table.index.tolist()


results = pd.DataFrame(all_winners, columns = ['ID', 'Disease_Category', 'max_IC'])

results.to_csv('output.csv', index=None)