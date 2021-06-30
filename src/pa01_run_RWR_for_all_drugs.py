# Created by woochanghwang at 25/01/2021
import pandas as pd
import networkx as nx
import toolbox.data_handler as dh

def get_targets_in_network(SIP_G, drugA_targets, drugB_targets):
    print(drugA_targets, drugB_targets)
    if drugA_targets == "NA" or drugB_targets == "NA":
        return "NA"
    else:
        nodes = SIP_G.nodes()
        drugA_targets = drugA_targets.split(',')
        drugB_targets = drugB_targets.split(',')
        drugA_targets = list(set(drugA_targets)&set(nodes))
        drugB_targets = list(set(drugB_targets) & set(nodes))
        return drugA_targets,drugB_targets

def main():
    threshold = 800
    hfh_result_df = pd.read_excel(
        "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance_druginfo_ATC_{}.xlsx".format(800),sheet_name="filtered")
    print(hfh_result_df)

    SIP_Graph = dh.load_obj(
        "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/network/SP_based/COVID_All_Structure_All_Shortest_Paths_graph_24hr")

    ########
    # Only for Dexamethasone, triamcinolone
    #######


    # hfh_drug_df = hfh_result_df[['Control_drugname','Similar_drugname']]
    # print(hfh_drug_df)
    hfh_result_df = hfh_result_df.fillna('NA')
    for index,row in hfh_result_df.iterrows():
        # if (row['Control_drugname'] != "Dexamethasone"):  continue
        print(row)
        drugA_targets,drugB_targets = get_targets_in_network(SIP_Graph,row['Targets_x_SIP'],row['Targets_y_SIP'])
        # print(drugA_targets)
        start_genes_for_PR = dict()
        for target in drugA_targets:
            start_genes_for_PR[target] = 1
        PR_score = nx.pagerank(SIP_Graph, personalization=start_genes_for_PR)
        PR_score_df = pd.DataFrame.from_dict(PR_score, orient='index', columns=['PageRank'])
        RWR_result_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Drug_Pagerank/Control_{}/COVID_HfH_{}_Pagerank.xlsx".format(threshold,row['Control_drugname'])
        PR_score_df.to_excel(RWR_result_addr)

        start_genes_for_PR = dict()
        for target in drugB_targets:
            start_genes_for_PR[target] = 1
        PR_score = nx.pagerank(SIP_Graph, personalization=start_genes_for_PR)
        PR_score_df = pd.DataFrame.from_dict(PR_score, orient='index', columns=['PageRank'])
        RWR_result_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Drug_Pagerank/Similar_{}/COVID_HfH_{}_Pagerank.xlsx".format(
            threshold,row['Similar_drugname'])
        PR_score_df.to_excel(RWR_result_addr)

def main_pre():
    hfh_result_df = pd.read_excel(
        "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance_druginfo_route.xlsx")
    print(hfh_result_df)

    SIP_Graph = dh.load_obj(
        "/Users/woochanghwang/PycharmProjects/LifeArc/COVID-19/result/network/SP_based/COVID_All_Structure_All_Shortest_Paths_graph_24hr")


    # hfh_drug_df = hfh_result_df[['Control_drugname','Similar_drugname']]
    # print(hfh_drug_df)
    hfh_result_df = hfh_result_df.fillna('NA')
    for index,row in hfh_result_df.iterrows():
        # print(row)
        drugA_targets,drugB_targets = get_targets_in_network(SIP_Graph,row['Targets_x'],row['Targets_y'])
        # print(drugA_targets)
        start_genes_for_PR = dict()
        for target in drugA_targets:
            start_genes_for_PR[target] = 1
        PR_score = nx.pagerank(SIP_Graph, personalization=start_genes_for_PR)
        PR_score_df = pd.DataFrame.from_dict(PR_score, orient='index', columns=['PageRank'])
        RWR_result_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Drug_Pagerank/Control/COVID_HfH_{}_Pagerank.xlsx".format(
            '_'.join(row['Control_drugname']))
        PR_score_df.to_excel(RWR_result_addr)

        start_genes_for_PR = dict()
        for target in drugB_targets:
            start_genes_for_PR[target] = 1
        PR_score = nx.pagerank(SIP_Graph, personalization=start_genes_for_PR)
        PR_score_df = pd.DataFrame.from_dict(PR_score, orient='index', columns=['PageRank'])
        RWR_result_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Drug_Pagerank/Similar/COVID_HfH_{}_Pagerank.xlsx".format(
            '_'.join(row['Similar_drugname']))
        PR_score_df.to_excel(RWR_result_addr)

        # network_property_df = pd.DataFrame(
        #     columns=['Eigen', 'Degree', 'RWR'])
        # for node in network_nodes:
        #     network_property_df.loc[node] = [eigenvector_centrality_dict[node], degree_centrality_dict[node],
        #                                      PR_score[node]]
        #
        # network_property_df.to_csv(result_save_dir)
        # print(network_property_df.head())
    # DDI_target_df = pd.read_excel(DDI_target_addr)

    # DDI_target_df['Distance'] = DDI_target_df.apply(lambda row:
    #                                                 calc_DDI_cloeset_distance(SIP_Graph,
    #                                                                   row['Targets_x'],
    #                                                                   row['Targets_y']),
    #                                                 axis=1)
    #
    # DDI_target_df['Cloest_min'] = DDI_target_df.apply(lambda row:
    #                                                   calc_DDI_cloeset_min_distance(SIP_Graph,
    #                                                                                 row['Targets_x'],
    #                                                                                 row['Targets_y']),
    #                                                   axis=1)
    # DDI_target_df['Shortest_avg'] = DDI_target_df.apply(lambda row:
    #                                                     calc_DDI_shortest_avg_distance(SIP_Graph,
    #                                                                                    row['Targets_x'],
    #                                                                                    row['Targets_y']),
    #                                                     axis=1)
    #
    # DDI_target_df.to_excel(DDI_distance_addr, index=False)


if __name__ == '__main__':
    main()
