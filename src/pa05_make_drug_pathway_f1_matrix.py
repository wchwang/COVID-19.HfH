# Created by woochanghwang at 01/02/2021

# Created by woochanghwang at 17/07/2020
'''
Make Drug - enrichmemt pathway score ( F1 score) matrix
all pathways: 206
- Remvoe Top level pathways(Metabolism of lipids, Apoptosis et al)
'''
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx

def make_matrix(hfh_result_df,whole_pa_list):

    pathway_addr_prefix = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Drug_pathways_meanTarget_800_filtered/enriched_reactome_pathways_{}.csv"

    pathway = whole_pa_list
    print(len(pathway), pathway[:3])

    sorted_pathway = sorted(pathway)
    print(sorted_pathway[:3])

    drug_names = hfh_result_df['Drug'].to_list()

    # print(drug_names)
    print(sorted_pathway)
    print(len(sorted_pathway))
    drug_pathway_matrix = []

    for drug in drug_names:
        drug_pathway_enriched_df = pd.read_csv(pathway_addr_prefix.format(drug))
        drug_pathway_values = [drug]
        # print(drug)
        # print(drug_pathway_enriched_df.head())
        for path in sorted_pathway:
            path_df = drug_pathway_enriched_df[drug_pathway_enriched_df['term_name'] == path][['precision', 'recall']]
            if len(path_df) > 0:
                recall = path_df.iloc[0]['recall']
                precision = path_df.iloc[0]['precision']
                f1_score = 2 * ((precision * recall) / (precision + recall))
                drug_pathway_values.append(f1_score)
            else:
                drug_pathway_values.append(0.0)
        drug_pathway_matrix.append(drug_pathway_values)

    drug_pathway_df = pd.DataFrame(drug_pathway_matrix, columns=['Drug_name'] + sorted_pathway)
    print(drug_pathway_df)

    # sns.clustermap(drug_pathway_df)
    #
    # plt.show()

    drug_pathway_df.to_csv(
        "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Drug_pathways_meanTarget_800_filtered/F1/hfh_drug_to_low_level_pathways_0.01_matrix.csv",
        index=False)

def get_whole_enriched_pathways(hfh_result_df):
    hfh_drugs = hfh_result_df['Drug'].tolist()
    pa_addr_pref = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Drug_pathways_meanTarget_800_filtered/enriched_reactome_pathways_{}.csv"
    whole_pa_list = []
    for drug in hfh_drugs:
        pa_df = pd.read_csv(pa_addr_pref.format(drug))
        pa_df = pa_df[pa_df['p_value'] < 0.01]  #filter
        pa_list = pa_df['term_name'].to_list()
        whole_pa_list+=pa_list

    whole_pa_list = list(set(whole_pa_list))
    print(len(whole_pa_list))

    return whole_pa_list

def make_reactome_hierachi():
    reactome_hierachi_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/Reactome/ReactomePathwaysRelation_HSA.txt"
    reactome_path_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/Reactome/ReactomePathways.txt"

    with open(reactome_hierachi_addr) as r_h_f:
        reactome_pair = [x.strip().split('\t') for x in r_h_f.readlines()]

    with open(reactome_path_addr) as r_p_f:
        reactome_path_list = [x.strip().split('\t') for x in r_p_f.readlines()]

    reactome_path_map_dict = dict()
    reactome_name_to_id_dict = dict()
    for path in reactome_path_list:
        if path[2] == 'Homo sapiens':
            reactome_path_map_dict[path[0]] = path[1]
            reactome_name_to_id_dict[path[1]] = path[0]
    reactome_diGraph = nx.DiGraph(reactome_pair)

    # print(reactome_diGraph.nodes)

    return reactome_diGraph, reactome_path_map_dict, reactome_name_to_id_dict


def find_reactome_parents(path_name):
    reactome_diGraph, reactome_path_map_dict, reactome_name_to_id_dict = make_reactome_hierachi()
    reactome_id = reactome_name_to_id_dict.get(path_name,'NA')
    if reactome_id != 'NA':
        ancestor_temr_id = list(nx.ancestors(reactome_diGraph,reactome_id))
        parents = [reactome_path_map_dict[id] for id in ancestor_temr_id]
        parents = ' | '.join(parents)
        return parents
    else:   return 'NA'

# def find_reactome_hierachi_top2(path_name):
#     reactome_diGraph, reactome_path_map_dict, reactome_name_to_id_dict = make_reactome_hierachi()
#     # if "R-HSA" not in path_name:
#     reactome_id = reactome_name_to_id_dict.get(path_name,'NA')
#     # else:
#     #     reactome_id = path_name
#     print(path_name, reactome_id)
#     if reactome_id != 'NA':
#         ancestor_temr_id = list(nx.ancestors(reactome_diGraph,reactome_id))
#         print(ancestor_temr_id)
#         if len(ancestor_temr_id) == 1 :
#             print("return:",path_name, reactome_path_map_dict[ancestor_temr_id[0]])
#             ancestor = reactome_path_map_dict[ancestor_temr_id[0]]
#             return ancestor
#
#         elif len(ancestor_temr_id) > 1:
#             # predecessors = nx.predecessor(reactome_diGraph,reactome_id)
#             # predecessor_id = reactome_name_to_id_dict.get(predecessors[0], 'NA')
#             predecessors = list(reactome_diGraph.predecessors(reactome_id))
#             print("predecessors:", predecessors)
#             # predecessor_id = reactome_name_to_id_dict.get(predecessors[0],'NA')
#             predecessor_name = reactome_path_map_dict[predecessors[0]]
#             find_reactome_hierachi_top2(predecessor_name)
#             # return 'NA'
#
#         # parents = ' | '.join(parents)
#     else:   return 'NA'

def find_reactome_hierachi_top2(path_name):
    reactome_diGraph, reactome_path_map_dict, reactome_name_to_id_dict = make_reactome_hierachi()
    reactome_id = reactome_name_to_id_dict.get(path_name,'NA')
    print(path_name, reactome_id)
    if reactome_id != 'NA':
        ancestor_temr_id = list(nx.ancestors(reactome_diGraph,reactome_id))
        print(ancestor_temr_id, len(ancestor_temr_id))
        if len(ancestor_temr_id) == 1 :
            # print("return:", reactome_path_map_dict[reactome_id])
            return_path = reactome_path_map_dict[reactome_id]
            return return_path

        elif len(ancestor_temr_id) > 1:
            # predecessors = nx.predecessor(reactome_diGraph,reactome_id)
            # predecessor_id = reactome_name_to_id_dict.get(predecessors[0], 'NA')
            predecessors = list(reactome_diGraph.predecessors(reactome_id))
            # print("pre:",predecessors)
            # predecessor_id = reactome_name_to_id_dict.get(predecessors[0],'NA')
            predecessor_name = reactome_path_map_dict[predecessors[0]]
            return find_reactome_hierachi_top2(predecessor_name)
            # return 'NA'

        # parents = ' | '.join(parents)
    else:   return 'NA'

def find_ancestor_enriched_pathways(whole_pathways):
    hfh_pathway_df= pd.DataFrame(whole_pathways,columns=['Pathways'])
    print(hfh_pathway_df)
    # hfh_pathway_df = hfh_pathway_df[:4]
    hfh_pathway_df['Parents'] = hfh_pathway_df['Pathways'].apply(find_reactome_parents)
    hfh_pathway_df['Top2'] = hfh_pathway_df['Pathways'].apply(find_reactome_hierachi_top2)
    print(hfh_pathway_df)
    hfh_pathway_df.to_excel(
        "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Drug_pathways_meanTarget_800_filtered/Reactome_pathways_hierachi_0.01_top2.xlsx")


def find_low_level_pathways(moa_to_pathway_df):
    # print(moa_to_pathway_df)
    top2_groupby = moa_to_pathway_df.groupby('Top2')['Pathways'].agg(list).reset_index(name='Pathways')
    top2_groupby_size = moa_to_pathway_df.groupby('Top2')['Pathways'].size().reset_index(name='counts')
    top2_pathways_size_df = top2_groupby.merge(top2_groupby_size)

    high_level_top2_df = top2_pathways_size_df[top2_pathways_size_df['counts']>1]
    high_level_top2 = high_level_top2_df['Top2'].to_list()
    high_level_top2 = [x.strip() for x in high_level_top2]
    # print(high_level_top2)

    # moa_to_pathway_df = moa_to_pathway_df[['Pathways','Parents','Top2','SOM','C10','SOM_MoA','Range']]
    # moa_to_pathway_df['Parents_counts'] = moa_to_pathway_df['Parents'].str.split('|').str.len()

    low_level_pathways_df = moa_to_pathway_df.dropna()
    # print(low_level_pathways_df)

    low_level_pathways_df = low_level_pathways_df[~low_level_pathways_df['Pathways'].isin(high_level_top2)]
    print(low_level_pathways_df)
    low_level_pathways_df.to_excel("/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Drug_pathways_meanTarget_800_filtered/low_level_pathways_0.01.xlsx")

    return low_level_pathways_df['Pathways'].to_list()



def main():
    threshold = 800
    hfh_result_df = pd.read_excel(
        "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance_druginfo_ATC_{}.xlsx".format(
            threshold), sheet_name="filter_drugs")
    # print(hfh_result_df)

    # whole_pa_list = get_whole_enriched_pathways(hfh_result_df)
    # find_ancestor_enriched_pathways(whole_pa_list)    # one time
    drug_to_pathway_df= pd.read_excel("/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Drug_pathways_meanTarget_800_filtered/Reactome_pathways_hierachi_0.01_top2.xlsx")
    low_level_pathways = find_low_level_pathways(drug_to_pathway_df)


    make_matrix(hfh_result_df,low_level_pathways)


if __name__ == '__main__':
    main()
