# Created by woochanghwang at 01/02/2021

import pandas as pd
import toolbox.network_utilities as nu
import toolbox.data_handler as dh

def calc_DDI_cloeset_distance(SIP_G, drugA_targets, drugB_targets):
    if drugA_targets == "NA" or drugB_targets == "NA":
        return "NA"
    else:
        nodes = SIP_G.nodes()
        drugA_targets = drugA_targets.split(',')
        # drugB_targets = drugB_targets.split(',')
        drugA_targets = list(set(drugA_targets)&set(nodes))
        drugB_targets = list(set(drugB_targets) & set(nodes))
        return nu.calculate_closest_distance(SIP_G,drugA_targets,drugB_targets)

def calc_DDI_cloeset_min_distance(SIP_G, drugA_targets, drugB_targets):
    if drugA_targets == "NA" or drugB_targets == "NA":
        return "NA"
    else:
        nodes = SIP_G.nodes()
        drugA_targets = drugA_targets.split(',')
        drugB_targets = drugB_targets.split(',')
        drugA_targets = list(set(drugA_targets)&set(nodes))
        drugB_targets = list(set(drugB_targets) & set(nodes))
        return nu.calculate_closest_min_distance(SIP_G,drugA_targets,drugB_targets)

def calc_DDI_shortest_avg_distance(SIP_G, drugA_targets, drugB_targets):
    if drugA_targets == "NA" or drugB_targets == "NA":
        return "NA"
    else:
        nodes = SIP_G.nodes()
        drugA_targets = drugA_targets.split(',')
        drugB_targets = drugB_targets.split(',')
        drugA_targets = list(set(drugA_targets)&set(nodes))
        drugB_targets = list(set(drugB_targets) & set(nodes))
        return nu.calculate_shortest_avg_distance(SIP_G,drugA_targets,drugB_targets)

def add_DDI_distance(DDI_target_df, DDI_distance_addr):
    SIP_Graph = dh.load_obj(
        "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/network/SP_based/COVID_All_Structure_All_Shortest_Paths_graph_24hr")
    # DDI_target_df = pd.read_excel(DDI_target_addr)
    DDI_target_df = DDI_target_df.fillna('NA')

    moderate_up_degs_df = pd.read_excel("/Users/woochanghwang/PycharmProjects/LifeArc/COVID19_vaccine/result/all_study_key_DEGs_only_Moderate_UP.xlsx")
    moderate_up_degs_df = moderate_up_degs_df[moderate_up_degs_df["All"]>1]
    moderate_up_degs = moderate_up_degs_df['Gene'].tolist()

    severe_up_degs_df = pd.read_excel("/Users/woochanghwang/PycharmProjects/LifeArc/COVID19_vaccine/result/all_study_key_DEGs_only_Severe_UP.xlsx")
    severe_up_degs_df = severe_up_degs_df[severe_up_degs_df['All']>1]
    severe_up_degs = severe_up_degs_df['Gene'].tolist()
    print(moderate_up_degs_df)

    DDI_target_df['CM_Closest_mean'] = DDI_target_df.apply(lambda row:
                                                    calc_DDI_cloeset_distance(SIP_Graph,
                                                                      row['Targets_x'],
                                                                      moderate_up_degs),
                                                    axis=1)
    DDI_target_df['CS_Closest_mean'] = DDI_target_df.apply(lambda row:
                                                           calc_DDI_cloeset_distance(SIP_Graph,
                                                                                     row['Targets_x'],
                                                                                     severe_up_degs),
                                                           axis=1)

    DDI_target_df['SM_Closest_mean'] = DDI_target_df.apply(lambda row:
                                                           calc_DDI_cloeset_distance(SIP_Graph,
                                                                                     row['Targets_y'],
                                                                                     moderate_up_degs),
                                                           axis=1)
    DDI_target_df['SS_Closest_mean'] = DDI_target_df.apply(lambda row:
                                                           calc_DDI_cloeset_distance(SIP_Graph,
                                                                                     row['Targets_y'],
                                                                                     severe_up_degs),
                                                           axis=1)

    # DDI_target_df['Closest_min'] = DDI_target_df.apply(lambda row:
    #                                                 calc_DDI_cloeset_min_distance(SIP_Graph,
    #                                                                           row['Targets_x'],
    #                                                                           row['Targets_y']),
    #                                                 axis=1)
    # DDI_target_df['Shortest_avg'] = DDI_target_df.apply(lambda row:
    #                                                 calc_DDI_shortest_avg_distance(SIP_Graph,
    #                                                                           row['Targets_x'],
    #                                                                           row['Targets_y']),
    #                                                 axis=1)

    DDI_target_df.to_excel(DDI_distance_addr,index=False)

def find_DD_shard_proteins( drugA_targets, drugB_targets):
    if drugA_targets == "NA" or drugB_targets == "NA":
        return "NA"
    else:
        # print(drugA_targets)
        drugA_targets = drugA_targets.split(',')
        # drugB_targets = drugB_targets.split(',')

        shard_proteins = list(set(drugA_targets)&set(drugB_targets))
        if len(shard_proteins) > 0: return ';'.join(shard_proteins)
        else:    return 'NA'

def share_proteins(DD_shared_df, DD_shared_addr):

    '''
    arunachalam	overmyer	jain	xu
    '''
    DD_shared_df = DD_shared_df.fillna('NA')

    covid_studies = ['arunachalam',	'overmyer',	'jain',	'xu', 'All']
    threshold = 0
    for study in covid_studies:
        if study == 'All': threshold = 1
        else: threshold = 0
        moderate_up_degs_df = pd.read_excel(
            "/Users/woochanghwang/PycharmProjects/LifeArc/COVID19_vaccine/result/all_study_key_DEGs_only_Moderate_UP.xlsx")
        moderate_up_degs_df = moderate_up_degs_df[moderate_up_degs_df[study] > threshold]
        moderate_up_degs = moderate_up_degs_df['Gene'].tolist()

        severe_up_degs_df = pd.read_excel(
            "/Users/woochanghwang/PycharmProjects/LifeArc/COVID19_vaccine/result/all_study_key_DEGs_only_Severe_UP.xlsx")
        severe_up_degs_df = severe_up_degs_df[severe_up_degs_df[study] > threshold]
        severe_up_degs = severe_up_degs_df['Gene'].tolist()
        print(moderate_up_degs_df)

        DD_shared_df["Moderate_shared_targetY_{}".format(study)] = DD_shared_df.apply(lambda row:
                                                               find_DD_shard_proteins(row['Extened_targets'],
                                                                                         moderate_up_degs),
                                                               axis=1)
        DD_shared_df["Severe_shared_targetY_{}".format(study)] = DD_shared_df.apply(lambda row:
                                                               find_DD_shard_proteins(row['Extened_targets'],
                                                                                         severe_up_degs),
                                                               axis=1)
    DD_shared_df.to_excel(DD_shared_addr, index=False)
def main():
    threshold = 800

    # HfH_drug_pair_targets_df = pd.read_excel(
    #     "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance_druginfo_ATC_PA_{}.xlsx".format(
    #         threshold))

    # HfH_drug_pair_targets_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance_druginfo_ATC_{}.xlsx".format(
    #     threshold)
    HfH_drug_pair_targets_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Drug_Pagerank/filtered_800/filtered_drugs_rwr_extened_targets_{}.xlsx".format(threshold)


    ###################
    ## 1. check distance
    ###################
    # add_DDI_distance(HfH_drug_pair_targets_df, HfH_drug_pair_targets_distance_addr)

    ####################
    ## 2. check share proteins
    ###################



    # HfH_drug_pair_targets_distance_shared_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance_druginfo_ATC_shared_{}.xlsx".format(
    #     threshold)
    HfH_drug_pair_targets_distance_shared_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Drug_Pagerank/filtered_800/filtered_drugs_rwr_extened_targets_shared_{}.xlsx".format(threshold)

    HfH_drug_pair_targets_distance_df = pd.read_excel(HfH_drug_pair_targets_addr)
    # HfH_drug_pair_targetsY_distance_shared_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance_druginfo_ATC_PA_DEGs_shared_{}.xlsx".format(
    #     threshold)
    share_proteins(HfH_drug_pair_targets_distance_df, HfH_drug_pair_targets_distance_shared_addr)


def main_pre():
    threshold = 800

    HfH_drug_pair_targets_df = pd.read_excel(
        "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance_druginfo_ATC_PA_{}.xlsx".format(
            threshold))

    HfH_drug_pair_targets_distance_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance_druginfo_ATC_PA_DEGs_{}.xlsx".format(
        threshold)


    ###################
    ## 1. check distance
    ###################
    # add_DDI_distance(HfH_drug_pair_targets_df, HfH_drug_pair_targets_distance_addr)

    ####################
    ## 2. check share proteins
    ###################



    HfH_drug_pair_targets_distance_shared_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance_druginfo_ATC_PA_DEGs_shared_{}.xlsx".format(
        threshold)

    HfH_drug_pair_targets_distance_df = pd.read_excel(HfH_drug_pair_targets_distance_shared_addr)
    HfH_drug_pair_targetsY_distance_shared_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance_druginfo_ATC_PA_DEGs_shared_{}.xlsx".format(
        threshold)
    share_proteins(HfH_drug_pair_targets_distance_df, HfH_drug_pair_targetsY_distance_shared_addr)

if __name__ == '__main__':
    main()
