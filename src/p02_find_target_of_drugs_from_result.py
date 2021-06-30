# Created by woochanghwang at 14/12/2020
import pandas as pd
import toolbox.data_handler as dh
import numpy as np

def make_HfH_drug_to_target(threshold):
    similar_drugs_df = pd.read_csv("../result/Similar_drugs/phase2_drugs_first_pass_filtering_with_drugbank_ids_v3.csv")
    print(similar_drugs_df)
    # threshold = 700

    stitch_drugbank_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/STITCH/2020.12/9606.protein_chemical.links.v5.0.drugbank.v5.1.6.target_symbol.s{}.onlyTarget.tsv".format(
        threshold)
    # stitch_drugbank_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/STITCH/9606.protein_chemical.links.v5.0.drugbank.target_symbol.s{}.tsv".format(threshold)
    stitch_drugbank_df = pd.read_csv(stitch_drugbank_addr, sep='\t')
    # print(stitch_drugbank_df)
    stitch_drugban_df_groupby = stitch_drugbank_df.groupby('DrugBank ID').agg(list).reset_index()
    print(stitch_drugban_df_groupby)
    stitch_pubchem_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/STITCH/2020.12/9606.protein_chemical.links.v5.0.pubchem.symbol.s{}.tsv".format(
        threshold)
    stitch_pubchem_df = pd.read_csv(stitch_pubchem_addr, sep='\t')

    # print(stitch_pubchem_df)
    stitch_pubchem_df_groupby = stitch_pubchem_df.groupby('Pubchem').agg(list).reset_index()
    print(stitch_pubchem_df_groupby)

    similar_drugs_df = similar_drugs_df[['chemblID', 'pref_name','ID', 'drugbank_id', 'Pubchem']]
    similar_drugs_drugbankT_df = pd.merge(left=similar_drugs_df, right=stitch_drugban_df_groupby,
                                          left_on='drugbank_id',
                                          right_on='DrugBank ID',
                                          how='left')
    print(similar_drugs_drugbankT_df)

    similar_drugs_drugbankT_stitchT_df = pd.merge(left=similar_drugs_drugbankT_df, right=stitch_pubchem_df_groupby,
                                                  left_on='Pubchem',
                                                  right_on='Pubchem',
                                                  how='left')

    similar_drugs_drugbankT_stitchT_df.to_csv(
        "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/data/COVID.HfH_9606.STITCH.v5.0.drugbank.v5.1.6.onlytarget_symbol.s{}.csv".format(
            threshold), index=False)
    # similar_drugs_list = similar_drugs_df['drugbank_id'].tolist()

def integrated_target(row):


    if (row['Target_x'] == 'NA') and (row['Target_y'] == 'NA'):
        return 'NA'
    if (row['Target_x'] != 'NA') & (row['Target_y'] =='NA'):
        return str(row['Target_x'])[1:-1].replace('\'','').replace(' ','')
    if (row['Target_x'] == 'NA') & (row['Target_y'] !='NA'):
        return str(row['Target_y'])[1:-1].replace('\'','').replace(' ','')


def make_integraed_target_file(HfH_drug_target_file_addr,HfH_drug_integrated_target_file_addr):
    HfH_drug_target_df = pd.read_csv(HfH_drug_target_file_addr)
    HfH_drug_target_df = HfH_drug_target_df.fillna('NA')
    # conditions = [
    #     (HfH_drug_target_df.Target_x.isnull()) & (HfH_drug_target_df.Target_y.isnull()),
    #     (HfH_drug_target_df.Target_x.isnull()) & ~(HfH_drug_target_df.Target_y.isnull()),
    #     ~(HfH_drug_target_df.Target_x.isnull()) & (HfH_drug_target_df.Target_y.isnull())]
    # choices = [np.nan,HfH_drug_target_df.Target_x,HfH_drug_target_df.Target_y]

    HfH_drug_target_df['Targets'] = HfH_drug_target_df.apply(integrated_target,axis=1)
    # HfH_drug_target_df['Targets'] = np.select(conditions,choices)
    print(HfH_drug_target_df)
    HfH_drug_target_df.to_csv(HfH_drug_integrated_target_file_addr,index=False)

def get_one_neighbor(graph, proteins):
    one_neighbors = []

    for protein in proteins:
        one_neighbors += [n for n in graph.neighbors(protein)]

    return list(set(one_neighbors))

def get_subgraph(targets, SIP_Graph):

    one_neighbours = get_one_neighbor(SIP_Graph, targets)
    target_subgraph = SIP_Graph.subgraph(one_neighbours)
    return target_subgraph

def find_targets_in_SIP_N_subnetwork(HfH_drug_integrated_target_file_addr,HfH_drug_target_in_SIP_dict_addr):
    SIP_Graph = dh.load_obj("/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/network/SP_based/COVID_All_Structure_All_Shortest_Paths_graph_24hr")
    SIP_nodes = SIP_Graph.nodes()
    # print(SIP_nodes)

    HfH_drugs_targets_df = pd.read_csv(HfH_drug_integrated_target_file_addr)
    HfH_drugs_IDs = HfH_drugs_targets_df['ID'].tolist()
    HfH_drugs_targets_df = HfH_drugs_targets_df.set_index('ID')
    HfH_drugs_targets_df = HfH_drugs_targets_df.fillna('NA')

    HfH_drugs_target_dict = dict()

    # drug_targets = HfH_drugs_targets_df.loc['DB11842', 'Targets'].tolist()
    # print(drug_targets)
    # print(HfH_drugs_targets_df.head(10))

    for drug in HfH_drugs_IDs:
        drug_targets = HfH_drugs_targets_df.loc[drug,'Targets']

        if drug_targets is 'NA':
            HfH_drugs_target_dict[drug] = drug_targets
        else:
            drug_targets = str(drug_targets)
            # print(drug_targets)
            drug_targets = drug_targets.split(',')
            drug_targets_in_SIP = list(set(drug_targets)&set(SIP_nodes))
            if len(drug_targets_in_SIP)==0:
                HfH_drugs_target_dict[drug] = 'NA'
                print("not in graph:",drug)
            else:
                target_subgraph = get_subgraph(drug_targets_in_SIP, SIP_Graph)
                HfH_drugs_target_dict[drug] = target_subgraph

    print(HfH_drugs_target_dict)
    dh.save_obj(HfH_drugs_target_dict, HfH_drug_target_in_SIP_dict_addr)

def main():
    threshold = 900
    HfH_drug_target_file_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/data/COVID.HfH_9606.STITCH.v5.0.drugbank.v5.1.6.onlytarget_symbol.s{}.csv".format(threshold)
    HfH_drug_integrated_target_file_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/data/COVID.HfH_9606.STITCH.v5.0.drugbank.v5.1.6.onlytarget_symbol.s{}.integrated.csv".format(
        threshold)

    # make_HfH_drug_to_target(threshold)
    #
    make_integraed_target_file(HfH_drug_target_file_addr,HfH_drug_integrated_target_file_addr)


    #########################
    HfH_drug_target_in_SIP_dict_addr =  "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/data/COVID.HfH_drug_to_target_SIP_subnetwork.s{}".format(threshold)
    find_targets_in_SIP_N_subnetwork(HfH_drug_integrated_target_file_addr, HfH_drug_target_in_SIP_dict_addr)
if __name__ == '__main__':
    main()
