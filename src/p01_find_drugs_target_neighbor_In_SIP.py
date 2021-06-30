# Created by woochanghwang at 18/12/2020
import toolbox.data_handler as dh
import pandas as pd
import numpy as np


def get_neighbor(key_proteins, SIP_G):
    key_proteins = str(key_proteins)
    if key_proteins == 'nan':
        return "NA"
    else:
        key_proteins = key_proteins.split(',')
        key_proteins_neighbor = []
        for key in key_proteins:
            key_proteins_neighbor.append(key)
            key_proteins_neighbor += [n for n in SIP_G.neighbors(key)]

        return ','.join(key_proteins_neighbor)



def main():
    threshold = 900
    HfH_drug_pair_targets_distance_addr = "../result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance_{}.xlsx".format(threshold)
    HfH_drugs_target_df = pd.read_excel(HfH_drug_pair_targets_distance_addr)
    print(HfH_drugs_target_df)
    HfH_drugA_target_df = HfH_drugs_target_df[['Control_id','Control_drugname','Targets_x_SIP']]
    print(HfH_drugA_target_df)
    HfH_drugB_target_df = HfH_drugs_target_df[['Similar_id','Similar_drugname','Targets_y_SIP']]
    print(HfH_drugB_target_df)
    HfH_drugA_target_df = HfH_drugA_target_df.drop_duplicates()
    HfH_drugB_target_df = HfH_drugB_target_df.drop_duplicates()
    # print(HfH_drugA_target_df)
    # print(HfH_drugB_target_df)
    HfH_drugA_target_df = HfH_drugA_target_df.rename(columns={
        'Control_id': 'Drug_ID',
        'Control_drugname': 'Name',
        'Targets_x_SIP': 'Targets'
    })
    HfH_drugB_target_df = HfH_drugB_target_df.rename(columns={
        'Similar_id': 'Drug_ID',
        'Similar_drugname': 'Name',
        'Targets_y_SIP': 'Targets'
    })
    HfH_drugAB_target_df = HfH_drugA_target_df.append(HfH_drugB_target_df)
    HfH_drugAB_target_df = HfH_drugAB_target_df.drop_duplicates()
    print(HfH_drugAB_target_df)

    SIP_Graph = dh.load_obj(
        "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/network/SP_based/COVID_All_Structure_All_Shortest_Paths_graph_24hr")

    HfH_drugAB_target_df['Neighbor'] = HfH_drugAB_target_df.apply(lambda row:
                                                                  get_neighbor(row['Targets'],
                                                                               SIP_Graph),
                                                                  axis=1)
    HfH_drugAB_target_df = HfH_drugAB_target_df.dropna()
    HfH_drugAB_target_df.to_csv("../result/Similar_drugs/Control_and_Similar_drugs_targets_neighbor_{}.tsv".format(threshold),
                                sep='\t',
                                index=False)



if __name__ == '__main__':
    main()
