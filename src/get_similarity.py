import traceback

import pandas as pd
from chembl_webresource_client.new_client import new_client

similarity_query = new_client.similarity
molecule = new_client.molecule

molecule_fields = ['molecule_chembl_id', 'pref_name', 'molecule_structures', 'similarity']

# region The list of compounds that are in 2+ phase clinical trials
clinical_trial_compounds = set()
for phase in (2, 3, 4):
    clinical_trial_compounds.update(
        [c[molecule_fields[0]] for c in molecule.filter(max_phase=phase).only(molecule_fields[:1])])
print('The number of 2+ clinical trial compounds: %s' % len(clinical_trial_compounds))
# endregion

# The list of drugs that are currently investigated for treating COVID-19
covid_drugs = pd.read_csv('clinical_trial_covid_drugs_with_smiles_drugbank.csv', index_col='name')
print('The number of COVID-19 drugs in clinical trials: %s' % len(covid_drugs))

# region Calculate the similarity between compounds and COVID-19 drugs
candidates = {}  # candidates[compound_chembl_id] = [compound_name, compound_smiles, [(drug_name, similarity)]
for drug_name, smile in covid_drugs.itertuples():
    try:
        similar_compounds = similarity_query.filter(smiles=smile, similarity=0).only(molecule_fields)

        n_similar_compounds = 0
        for similar_compound_dict in similar_compounds:
            chembl_id, comp_name, structure, similarity = [similar_compound_dict[f] for f in molecule_fields]

            if chembl_id in clinical_trial_compounds:
                candidates.setdefault(chembl_id, [comp_name, structure['canonical_smiles'], []])
                candidates[chembl_id][2].append((drug_name, float(similarity)))
                n_similar_compounds += 1

        print('%s/%s 2+ clinical trial compounds similar to covid drug %s'
              % (n_similar_compounds, len(similar_compounds), drug_name))

    except:
        traceback.print_exc()
# endregion

# Save the results
candidates = pd.DataFrame.from_dict(candidates, orient='index',
                                    columns=['pref_name', 'canonical_smiles', 'similar_covid_drugs'])
candidates['n_covid_drugs'] = candidates['similar_covid_drugs'].apply(len)
candidates['max_similarity'] = candidates['similar_covid_drugs'].apply(lambda l: max([t[1] for t in l]))
candidates.to_csv('phase2_drugs_with_similarity_09012021.tsv', sep='\t')
