# Created by woochanghwang at 16/12/2020
import pandas as pd
import networkx as nx
import toolbox.network_utilities as nu
import toolbox.data_handler as dh

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

def calc_drug_similarity_by_network(control_id, similar_id, HfH_drug_target_dict,similar_method):
    control_G = HfH_drug_target_dict[control_id]
    case_G = HfH_drug_target_dict[similar_id]
    if (control_G == 'NA') or (case_G == 'NA'):
        return 'NA'
    else:
        return nu.calc_network_similarity(control_G,case_G,method=similar_method)


def main():
    threshold = 900
    HfH_drug_target_in_SIP_dict_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/data/COVID.HfH_drug_to_target_SIP_subnetwork.s{}".format(
        threshold)
    HfH_drug_target_in_SIP_dict = dh.load_obj(HfH_drug_target_in_SIP_dict_addr)

    # control_case_pair_df = pd.read_csv("/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID.csv")
    # print(control_case_pair_df)

    # control_case_pair_df = pd.read_excel("/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_{}.xlsx".format(threshold))
    control_case_pair_df = pd.read_excel(
        "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance_{}.xlsx".format(
            threshold))

    control_case_pair_df['Similarity_VEO'] = control_case_pair_df.apply(lambda row:
                                                                    calc_drug_similarity_by_network(row['Control_id'],
                                                                                                    row['Similar_id'],
                                                                                                    HfH_drug_target_in_SIP_dict,
                                                                                                    'VEO'),
                                                                    axis=1)
    control_case_pair_df['Similarity_EO'] = control_case_pair_df.apply(lambda row:
                                                                        calc_drug_similarity_by_network(
                                                                            row['Control_id'],
                                                                            row['Similar_id'],
                                                                            HfH_drug_target_in_SIP_dict,
                                                                            'EO'),
                                                                        axis=1)
    # control_case_pair_df['Similarity_GED'] = control_case_pair_df.apply(lambda row:
    #                                                                     calc_drug_similarity_by_network(
    #                                                                         row['Control_id'],
    #                                                                         row['Similar_id'],
    #                                                                         HfH_drug_target_in_SIP_dict,
    #                                                                         'GED'),
    #                                                                     axis=1)
    # nu.calc_network_similarity()
    control_case_pair_df['Similarity_EOVL'] = control_case_pair_df.apply(lambda row:
                                                                       calc_drug_similarity_by_network(
                                                                           row['Control_id'],
                                                                           row['Similar_id'],
                                                                           HfH_drug_target_in_SIP_dict,
                                                                           'EOVL'),
                                                                       axis=1)
    print(control_case_pair_df)

    control_case_pair_df.to_excel("/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_{}.xlsx".format(threshold),
                                  index=False)
if __name__ == '__main__':
    main()
