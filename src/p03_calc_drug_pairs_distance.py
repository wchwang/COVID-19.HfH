# Created by woochanghwang at 18/12/2020
import pandas as pd
import toolbox.network_utilities as nu
import toolbox.data_handler as dh

def check_targets_in_Graph(SIP_G, targets):
    if targets == "NA":
        return "NA"
    else:
        nodes = SIP_G.nodes()
        targets= str(targets).split(',')
        targets = list(set(targets)&set(nodes))
        return ','.join(targets)

def calc_num_targets_in_Graph(targets):
    targets = str(targets)
    if targets == "NA":
        return "NA"
    else:
        targets= targets.split(',')
        return len(targets)

def calc_shared_targets(drugA_targets, drugB_targets):
    if drugA_targets == "NA" or drugB_targets == "NA":
        return "NA"
    else:

        drugA_targets = drugA_targets.split(',')
        drugB_targets = drugB_targets.split(',')

        return len(set(drugA_targets)&set(drugB_targets))
def calc_DDI_cloeset_distance(SIP_G, drugA_targets, drugB_targets):
    if drugA_targets == "NA" or drugB_targets == "NA":
        return "NA"
    else:
        nodes = SIP_G.nodes()
        drugA_targets = drugA_targets.split(',')
        drugB_targets = drugB_targets.split(',')
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
    DDI_target_df['Closest_mean'] = DDI_target_df.apply(lambda row:
                                                    calc_DDI_cloeset_distance(SIP_Graph,
                                                                      row['Targets_x'],
                                                                      row['Targets_y']),
                                                    axis=1)

    DDI_target_df['Closest_min'] = DDI_target_df.apply(lambda row:
                                                    calc_DDI_cloeset_min_distance(SIP_Graph,
                                                                              row['Targets_x'],
                                                                              row['Targets_y']),
                                                    axis=1)
    DDI_target_df['Shortest_avg'] = DDI_target_df.apply(lambda row:
                                                    calc_DDI_shortest_avg_distance(SIP_Graph,
                                                                              row['Targets_x'],
                                                                              row['Targets_y']),
                                                    axis=1)

    DDI_target_df.to_excel(DDI_distance_addr,index=False)

def main():
    threshold = 700
    HfH_drug_integrated_target_file_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/data/COVID.HfH_9606.STITCH.v5.0.drugbank.v5.1.6.onlytarget_symbol.s{}.integrated.csv".format(
        threshold)

    HfH_drugs_targets_df = pd.read_csv(HfH_drug_integrated_target_file_addr)
    print(HfH_drugs_targets_df)
    drug_pair_addr= "../result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_v2.xlsx"
    HfH_drug_pair_df = pd.read_excel(drug_pair_addr)
    print(HfH_drug_pair_df)

    HfH_drug_pair_targets_df= pd.merge(left=HfH_drug_pair_df,
                                       right=HfH_drugs_targets_df,
                                       left_on='Control_id',
                                       right_on='ID',
                                       how='left')
    HfH_drug_pair_targets_df=HfH_drug_pair_targets_df.drop(columns=['chemblID', 'pref_name', 'ID', 'drugbank_id', 'Pubchem', 'DrugBank ID', 'Target_x', 'Target_y'])
    print(list(HfH_drug_pair_targets_df))

    HfH_drug_pair_targets_df = pd.merge(left=HfH_drug_pair_targets_df,
                                        right=HfH_drugs_targets_df,
                                        left_on='Similar_id',
                                        right_on='ID',
                                        how='left')
    HfH_drug_pair_targets_df = HfH_drug_pair_targets_df.drop(
        columns=['chemblID', 'pref_name', 'ID', 'drugbank_id', 'Pubchem', 'DrugBank ID', 'Target_x', 'Target_y'])
    print(list(HfH_drug_pair_targets_df))
    SIP_Graph = dh.load_obj(
        "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/network/SP_based/COVID_All_Structure_All_Shortest_Paths_graph_24hr")

    HfH_drug_pair_targets_df['Targets_x_SIP'] = HfH_drug_pair_targets_df.apply(lambda row:
                                                                             check_targets_in_Graph(SIP_Graph,
                                                                                                    row['Targets_x']),
                                                                             axis=1)
    HfH_drug_pair_targets_df['Targets_y_SIP'] = HfH_drug_pair_targets_df.apply(lambda row:
                                                                               check_targets_in_Graph(SIP_Graph,
                                                                                                      row['Targets_y']),
                                                                               axis=1)

    HfH_drug_pair_targets_df['N_Targets_x_SIP'] = HfH_drug_pair_targets_df.apply(lambda row:
                                                                               calc_num_targets_in_Graph(row['Targets_x_SIP']),
                                                                               axis=1)
    HfH_drug_pair_targets_df['N_Targets_y_SIP'] = HfH_drug_pair_targets_df.apply(lambda row:
                                                                               calc_num_targets_in_Graph(row['Targets_y_SIP']),
                                                                               axis=1)
    HfH_drug_pair_targets_df['Shared_targets'] = HfH_drug_pair_targets_df.apply(lambda row:
                                                                               calc_shared_targets(row['Targets_x_SIP'],
                                                                                                   row['Targets_y_SIP']),
                                                                               axis=1)
    HfH_drug_pair_targets_distance_addr = "../result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance_{}.xlsx".format(threshold)
    add_DDI_distance(HfH_drug_pair_targets_df, HfH_drug_pair_targets_distance_addr)

    #################################################
    # To add various shortest
    #################################################

    # hfh_result_df = pd.read_excel("/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance_druginfo_route.xlsx")
    # print(hfh_result_df)
    # hfh_result_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_various_distance.xlsx"
    # add_DDI_distance(hfh_result_df, hfh_result_addr )



if __name__ == '__main__':
    main()
