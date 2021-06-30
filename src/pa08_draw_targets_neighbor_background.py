# Created by woochanghwang at 20/02/2021

import pandas as pd
import networkx as nx
import toolbox.data_handler as dh
import random
import seaborn as sns
import matplotlib.pyplot as plt

def check_target_exist(iterables, key_to_find):
    key_list = list(filter(lambda x: key_to_find in x, iterables))
    # print(key_list)
    # print(len(key_list))
    # return key_list
    if len(key_list) > 0:
        return iterables
    else:
        return ['NA']

def get_targets_a_drug(drug):
    threshold = 800
    HfH_drug_pair_targets_df = pd.read_excel(
        "/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance_druginfo_ATC_PA_{}.xlsx".format(
            threshold))
    # print(HfH_drug_pair_targets_df)
    a_drug_df = HfH_drug_pair_targets_df[HfH_drug_pair_targets_df['Similar_drugname']==drug]
    a_drug_targets = a_drug_df['Targets_y_SIP'].tolist()
    a_drug_targets = a_drug_targets[0].split(',')
    # print(a_drug_df['Targets_y'])

    return a_drug_targets

def get_COVID_Jain_DEGs():
    jain_DEGs_moderate = pd.read_excel("/Users/woochanghwang/PycharmProjects/MTIProject/COVID19_vaccine/data/Jain/Jain_CSBJ_supple2_v2.xlsx",sheet_name="Control_moderate")
    # print(jain_DEGs_moderate)
    jain_DEGs_moderate = jain_DEGs_moderate[['Gene','log2FoldChange']]
    jain_DEGs_moderate['COVID'] = "Moderate"
    jain_DEGs = jain_DEGs_moderate
    # print(jain_DEGs_moderate)
    # jain_DEGs_severe = pd.read_excel(
    #     "/Users/woochanghwang/PycharmProjects/LifeArc/COVID19_vaccine/data/Jain/Jain_CSBJ_supple2.xlsx",
    #     sheet_name="Control_Severe")
    # # print(jain_DEGs_severe)
    # jain_DEGs_severe = jain_DEGs_severe[['Gene','log2FoldChange']]
    # jain_DEGs_severe['COVID'] = "Severe"
    # jain_DEGs = jain_DEGs_moderate.append(jain_DEGs_severe)
    #
    # jain_DEGs_mild = pd.read_excel(
    #     "/Users/woochanghwang/PycharmProjects/LifeArc/COVID19_vaccine/data/Jain/Jain_CSBJ_supple2.xlsx",
    #     sheet_name="Control_mild")
    # # print(jain_DEGs_severe)
    # jain_DEGs_mild = jain_DEGs_mild[['Gene', 'log2FoldChange']]
    # jain_DEGs_mild['COVID'] = "Mild"
    # jain_DEGs = jain_DEGs.append(jain_DEGs_mild)
    #
    jain_DEGs = jain_DEGs.reset_index()
    jain_DEGs = jain_DEGs.rename(columns={
        'log2FoldChange':'log2FC'
    })
    jain_DEGs = jain_DEGs[['Gene','log2FC','COVID']]
    print(jain_DEGs)

    return jain_DEGs

def get_COVID_Arun_DEGs():
    arun_DEGs_df= pd.read_excel("/Users/woochanghwang/PycharmProjects/MTIProject/COVID19_vaccine/data/Arunachalam/Arunachalam_DEGs_results.xlsx")
    arun_DEGs_moderate = arun_DEGs_df[['Symbol','log2FoldChange_MH']]
    arun_DEGs_moderate = arun_DEGs_moderate.rename(columns={ 'Symbol':'Gene',
                                                             'log2FoldChange_MH' : 'log2FC'})
    arun_DEGs_moderate['COVID'] ='Moderate'
    return arun_DEGs_moderate


    # arun_DEGs_severe = arun_DEGs_df[['Symbol', 'log2FoldChange_SH']]
    # arun_DEGs_severe = arun_DEGs_severe.rename(columns={'Symbol': 'Gene',
    #                                                         'log2FoldChange_SH': 'log2FC'})
    # arun_DEGs_severe['COVID'] = 'Severe'
    #
    # arun_DEGs = arun_DEGs_moderate.append(arun_DEGs_severe)

    # return arun_DEGs

def get_COVID_Overmyer_DEGs():

    DEGs_moderate= pd.read_excel("/Users/woochanghwang/PycharmProjects/MTIProject/COVID19_vaccine/data/Overmyer/Overmyer_DEGs_M_S.xlsx",
                                      sheet_name="Moderate")
    DEGs_moderate = DEGs_moderate[['Gene','log2FC']]
    # DEGs_moderate = DEGs_moderate.rename(columns={ 'Symbol':'Gene',
    #                                                          'log2FoldChange_MH' : 'log2FC'})
    DEGs_moderate['COVID'] ='Moderate'

    return DEGs_moderate

def get_COVID_xu_DEGs():
    DEGs_moderate = pd.read_excel("/Users/woochanghwang/PycharmProjects/MTIProject/COVID19_vaccine/data/Xu/Xu_DEGs_scRNA.xlsx")
    DEGs_moderate = DEGs_moderate[DEGs_moderate['compare']=='MH']
    DEGs_moderate = DEGs_moderate[['gene','avg_logFC']]
    DEGs_moderate = DEGs_moderate.rename(columns={ 'gene':'Gene',
                                                    'avg_logFC' : 'log2FC'})
    DEGs_moderate['COVID'] = 'Moderate'

    return DEGs_moderate

def get_paths_with_targets_N_control(drug):
    all_shortest_results = dh.load_obj(
        "/Users/woochanghwang/PycharmProjects/MTIProject/COVID-19/result/network/SP_based/COVID_All_Structure_All_Shortest_Paths_24hr")

    all_shortest_paths = []
    for band, shortest_path in all_shortest_results.items():
        list_paths = list(shortest_path)
        list_paths = [list(path) for path in list_paths]
        all_shortest_paths += list_paths
    #
    all_shortest_paths = list(map(tuple, all_shortest_paths))
    # drug_names = ['gallopamil', 'triamcinolone']
    # drug_targets = get_targets_a_drug(drug_names[1])
    drug_targets = get_targets_a_drug(drug)

    # paths_including_targets = []
    #
    # for target in drug_targets:
    #     print("target:", target)
    #     for path in all_shortest_paths:
    #         path_containing_target = check_target_exist(path, target)
    #         # print(path_containing_target)
    #         if path_containing_target[0] != 'NA':
    #             paths_including_targets.append(path_containing_target)
    #
    # print(len(paths_including_targets))
    #
    # paths_not_in_targets = list(set(all_shortest_paths) - set(paths_including_targets))
    #
    # print(paths_not_in_targets[:3])
    # random_paths = random.sample(paths_not_in_targets, len(paths_including_targets))
    # print(random_paths[:2])
    #
    # paths_including_targets_edges = []
    # for path in paths_including_targets:
    #     paths_including_targets_edges += path
    #
    # paths_random_edges = []
    # for path in random_paths:
    #     paths_random_edges += path
    #
    # graph_including_targets = nx.Graph(paths_including_targets_edges)
    # random_graph = nx.Graph(paths_random_edges)
    #
    # nodes_in_graph_with_targets = graph_including_targets.nodes()
    # nodes_in_random_graph = random_graph.nodes()
    #
    # nodes_from_SIP = dict()
    # nodes_from_SIP['withTarget'] = nodes_in_graph_with_targets
    # nodes_from_SIP['Control'] =  nodes_in_random_graph
    #
    # return nodes_from_SIP

    ##############################
    ## per target
    #############################
    nodes_from_SIP = dict()

    for target in drug_targets:
        print("target:", target)
        paths_including_atarget = []
        for path in all_shortest_paths:
            path_containing_target = check_target_exist(path, target)
            # print(path_containing_target)
            if path_containing_target[0] != 'NA':
                paths_including_atarget.append(path_containing_target)

        print(len(paths_including_atarget))

        paths_not_in_targets = list(set(all_shortest_paths) - set(paths_including_atarget))

        print(paths_not_in_targets[:3])
        random_paths = random.sample(paths_not_in_targets, len(paths_including_atarget))
        print(random_paths[:2])

        paths_including_targets_edges = []
        for path in paths_including_atarget:
            paths_including_targets_edges += path

        paths_random_edges = []
        for path in random_paths:
            paths_random_edges += path

        graph_including_targets = nx.Graph(paths_including_targets_edges)
        random_graph = nx.Graph(paths_random_edges)

        nodes_in_graph_with_targets = graph_including_targets.nodes()
        nodes_in_random_graph = random_graph.nodes()

        # ############
        # ## exclude same nodes
        # ############
        #
        # nodes_in_random_graph = list(set(nodes_in_random_graph) - set( nodes_in_graph_with_targets))
        #
        # ##########
        #
        # nodes_from_SIP[target] = nodes_in_graph_with_targets
        # nodes_from_SIP['{}_Control'.format(target)] =  nodes_in_random_graph
        #
        # print("Targeet & random:", len(set(nodes_in_random_graph)&set(nodes_in_graph_with_targets)),
        #       len(nodes_in_random_graph),
        #       len(nodes_in_graph_with_targets))

        ################
        # one-neighbor
        ################

        targets_neighbor = list(graph_including_targets.neighbors(target))
        targets_neighbor.append(target)
        print(targets_neighbor)

        # flag = True
        # while flag:
        #     random_targets = random.sample(nodes_in_random_graph, len(targets_neighbor))
        #     if len(set(targets_neighbor)&set(random_targets)) == 0:
        #         flag = False


        nodes_from_SIP[target] = targets_neighbor
        # nodes_from_SIP['{}_Control'.format(target)] = random_targets

        print("target, random:", len(targets_neighbor))

    return nodes_from_SIP

def main_box_hist_sep():
    drug_names = ['gallopamil', 'triamcinolone']
    drug = drug_names[1]
    nodes_from_SIP = get_paths_with_targets_N_control(drug)

    # print(nodes_from_SIP)

    # covid_degs_df = get_COVID_Jain_DEGs()
    covid_degs_df = get_COVID_Arun_DEGs()
    # covid_degs_df = get_COVID_Overmyer_DEGs()
    # covid_degs_df = get_COVID_xu_DEGs()

    covid_degs_path_df = pd.DataFrame(columns=['Gene', 'log2FC', 'COVID'])

    #############################################################
    # for path,genes in nodes_from_SIP.items():
    #     print(path)A
    #     covid_path_df = covid_degs_df[covid_degs_df['Gene'].isin(genes)]
    #     covid_path_df['Target'] = path
    #     # print(covid_path_df)
    #
    #     covid_degs_path_df = covid_degs_path_df.append(covid_path_df, ignore_index=True)
    #     # print(covid_degs_path_df)
    #
    # print(covid_degs_path_df)
    # path_groups_sorted = covid_degs_path_df.groupby(['Target']).median().sort_values(ascending=False,
    #                                                                          by='log2FC').index.tolist()
    # covid_degs_path_df = covid_degs_path_df
    # # covid_degs_path_df.to_excel("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/pathway_DEGs/COVID_Jain_DEGs_Pathways_for_twodrugs.xlsx",index=False)
    # # plt.figure(figsize=(70, 15))
    # ax = sns.boxplot(x="Target", y="log2FC", hue="COVID",
    #                  data=covid_degs_path_df,order=path_groups_sorted)
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    # plt.tight_layout()
    #
    # plt.savefig("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/paths_with_target_DEGs/{}_drugs_target_neighbor_FC_boxplot_MS_Arun_median_noControl.pdf".format(drug))
    # plt.show()
    ##############################################################

    ############################################################
    # Tartet neithbot vs All genes
    #############################################################
    for path,genes in nodes_from_SIP.items():
        # print(path)
        covid_path_df = covid_degs_df[covid_degs_df['Gene'].isin(genes)]
        covid_path_df['Target'] = path
        # print(covid_path_df)

        covid_degs_path_df = covid_degs_path_df.append(covid_path_df, ignore_index=True)
        # print(covid_degs_path_df)

    #########
    # Add background
    ##########
    covid_backgroud_df = covid_degs_df
    covid_backgroud_df['Target'] = "Backgroud"
    ######
    ## check outlier
    ######
    covid_backgroud_df = covid_backgroud_df.dropna()
    # covid_backgroud_df = covid_backgroud_df[~covid_backgroud_df['Gene'].str.contains(".")]
    covid_backgroud_df = covid_backgroud_df[(covid_backgroud_df['log2FC']>-5)&(covid_backgroud_df['log2FC']<5)]
    print(covid_backgroud_df)

    covid_degs_path_df = covid_degs_path_df.append(covid_backgroud_df, ignore_index=True) # add backgroud into target df
    ##########
    print(covid_degs_path_df)
    path_groups_sorted = covid_degs_path_df.groupby(['Target']).median().sort_values(ascending=False,
                                                                             by='log2FC').index.tolist()
    covid_degs_path_df = covid_degs_path_df

    #########################
    ## Draw Boxplot
    ##########
    # covid_degs_path_df.to_excel("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/pathway_DEGs/COVID_Jain_DEGs_Pathways_for_twodrugs.xlsx",index=False)
    # plt.figure(figsize=(70, 15))
    # ax = sns.boxplot(y="Target", x="log2FC", hue="COVID",
    #                  data=covid_degs_path_df,order=path_groups_sorted)
    ax = sns.boxplot(y="Target", x="log2FC", data=covid_degs_path_df)

    # ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    plt.setp(ax.artists, edgecolor='k', facecolor='w')
    plt.setp(ax.lines, color='k')
    plt.tight_layout()
    plt.savefig(
        "/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/paths_with_target_DEGs/{}_drugs_FC_boxplot_MS_Arun_background_BW.pdf".format(
            drug))
    plt.show()
    #################################
    ## Histogram
    ###############################

    # ax = sns.displot(data=covid_backgroud_df, x="log2FC", hue="Target",kind="kde", fill=True)
    # ax = sns.displot(data=covid_backgroud_df, x="log2FC", hue="Target")
    ax = sns.histplot(data=covid_backgroud_df, x="log2FC", bins=50)
    plt.tight_layout()
    plt.savefig("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/paths_with_target_DEGs/drugs_FC_displot_MS_Arun_background_bin50.pdf".format(drug))
    plt.show()
    ##############################################################

    #############################################################
    ## Check target neighbor & RNA-Seq
    #############################################################

    # background_genes = covid_backgroud_df['Gene'].tolist()
    # targets = covid_degs_path_df['Target'].tolist()
    # targets = list(set(targets)-set(["Background"]))
    #
    # for target in targets:
    #     print("Target:",target , len(target))
    #     target_genes = covid_degs_path_df[covid_degs_path_df['Target']==target]['Gene'].tolist()
    #     target_in_rnaseq = list(set(target_genes)&set(background_genes))
    #     print("Target: {}".format(len(target_genes)))
    #     print("Target & Backgoud:{}".format(len(target_in_rnaseq)))
    #     # print("genelist:",target_in_rnaseq)



    # #################################################################
    # # box + histogram
    # #################################################################
    # f, (ax_box, ax_hist) = plt.subplots(2, sharex=True, gridspec_kw={"height_ratios": (.15, .85)})
    #
    # # Add a graph in each part
    # target_list = get_targets_a_drug(drug)
    #
    # # sns.boxplot(y="Target", x="log2FC", hue="COVID", data=covid_degs_path_df, ax=ax_box)
    # # sns.histplot(data=covid_backgroud_df, x="log2FC", hue="Target", ax=ax_hist)
    # sns.boxplot(y="Target", x="log2FC", data=covid_degs_path_df, ax=ax_box)
    # sns.histplot(data=covid_backgroud_df, x="log2FC", ax=ax_hist)
    #
    # # Remove x axis name for the boxplot
    # ax_box.set(xlabel='')
    #
    # plt.show()

def main():
    '''
    Draw boxplots on the histogram
    :return:
    '''


    # covid_degs_df = get_COVID_Jain_DEGs()
    covid_degs_df = get_COVID_Arun_DEGs()
    # covid_degs_df = get_COVID_Overmyer_DEGs()
    # covid_degs_df = get_COVID_xu_DEGs()

    drug_names = ['gallopamil', 'triamcinolone']
    drug = drug_names[0]
    nodes_from_SIP = get_paths_with_targets_N_control(drug)

    covid_degs_path_gall_df = pd.DataFrame(columns=['Gene', 'log2FC', 'COVID'])


    ############################################################
    # Tartet neithbot vs All genes  , Gallopamil
    #############################################################
    for path,genes in nodes_from_SIP.items():
        # print(path)
        covid_path_df = covid_degs_df[covid_degs_df['Gene'].isin(genes)]
        covid_path_df['Target'] = path
        # print(covid_path_df)

        covid_degs_path_gall_df = covid_degs_path_gall_df.append(covid_path_df, ignore_index=True)
        # print(covid_degs_path_df)

    ####################################
    # drug: Triamcinolone
    ####################################
    drug = drug_names[1]
    nodes_from_SIP = get_paths_with_targets_N_control(drug)

    covid_degs_path_tria_df = pd.DataFrame(columns=['Gene', 'log2FC', 'COVID'])

    ############################################################
    # Tartet neithbot vs All genes  , Gallopamil
    #############################################################
    for path, genes in nodes_from_SIP.items():
        # print(path)
        covid_path_df = covid_degs_df[covid_degs_df['Gene'].isin(genes)]
        covid_path_df['Target'] = path
        # print(covid_path_df)

        covid_degs_path_tria_df = covid_degs_path_tria_df.append(covid_path_df, ignore_index=True)
        # print(covid_degs_path_df)

    #########
    # Add background
    ##########
    covid_backgroud_df = covid_degs_df
    covid_backgroud_df['Target'] = "Backgroud"
    ## check outlier
    covid_backgroud_df = covid_backgroud_df.dropna()
    # covid_backgroud_df = covid_backgroud_df[~covid_backgroud_df['Gene'].str.contains(".")]
    covid_backgroud_df = covid_backgroud_df[(covid_backgroud_df['log2FC']>-5)&(covid_backgroud_df['log2FC']<5)]
    print(covid_backgroud_df)



    #################################################################
    # box + histogram
    #################################################################
    f, (ax_box_gall,ax_box_tria, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (.15, .15, .70)})

    # Add a graph in each part
    # plt.figure(figsize=((30, 30)))

    # sns.boxplot(y="Target", x="log2FC", hue="COVID", data=covid_degs_path_df, ax=ax_box)
    # sns.histplot(data=covid_backgroud_df, x="log2FC", hue="Target", ax=ax_hist)
    sns.boxplot(y="Target", x="log2FC", data=covid_degs_path_gall_df, ax=ax_box_gall, palette="Set2")
    sns.boxplot(y="Target", x="log2FC", data=covid_degs_path_tria_df, ax=ax_box_tria)
    sns.histplot(data=covid_backgroud_df, x="log2FC",bins=50, ax=ax_hist)

    # Remove x axis name for the boxplot
    ax_box_gall.set(xlabel='')
    ax_box_tria.set(xlabel='')
    # make boxplot color black and white
    plt.setp(ax_box_gall.artists, edgecolor='k', facecolor='w')
    plt.setp(ax_box_gall.lines, color='k')
    plt.setp(ax_box_tria.artists, edgecolor='k', facecolor='w')
    plt.setp(ax_box_tria.lines, color='k')
    plt.tight_layout()
    # plt.savefig("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/paths_with_target_DEGs/two_drugs_targetFC_boxplot_on_backgroud_histplot_MS_Arun.pdf")
    # f.savefig("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/paths_with_target_DEGs/two_drugs_targetFC_boxplot_on_backgroud_histplot_MS_Arun_BW_bin50.pdf")
    plt.show()

    covid_degs_path_tria_df.to_excel("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/one_neighbor/triamcinolone_target_neighbor_covid_Arun_moderate.xlsx")
    covid_degs_path_gall_df.to_excel("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/one_neighbor/gallopamil_target_neighbor_covid_Arun_moderate.xlsx")


if __name__ == '__main__':
    main()
    # main_box_hist_sep()