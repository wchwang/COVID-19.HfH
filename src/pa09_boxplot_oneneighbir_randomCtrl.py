# Created by woochanghwang at 02/03/2021

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


def get_target_RWR_result_drugs(drugs):
    drug_target_df = pd.read_excel(
        "/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/Drug_Pagerank/filtered_800/filtered_drugs_rwr_extened_targets_800.xlsx")

    drug_target_df = drug_target_df[drug_target_df['Drug'].isin(drugs)]
    # print(drug_target_df)
    return drug_target_df

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



def make_boxplot_input(target_neighbor_omics_dict):
    covid_logfc_all_targets_df = pd.DataFrame(columns=['log2FC'])

    ############################################################
    for target, values in target_neighbor_omics_dict.items():
        print(target)
        if 'Control' in target:
            # print(values[:2])
            combine_values = [value for set in values for value in set]
            # combine_values = [sum(value)/len(value) for value in values]
            covid_logfc_df = pd.DataFrame(combine_values, columns=['log2FC'])
        else:
            covid_logfc_df = pd.DataFrame(values, columns=['log2FC'])

        covid_logfc_df['Target'] = target
        # print(covid_path_df)

        covid_logfc_all_targets_df = covid_logfc_all_targets_df.append(covid_logfc_df, ignore_index=True)
        # print(covid_logfc_all_targets_df)

    print(covid_logfc_all_targets_df)

    return covid_logfc_all_targets_df

def draw_boxplot(drug,covid_logfc_all_targets_df, save_addr):

    path_groups_sorted = covid_logfc_all_targets_df.groupby(['Target']).median().sort_values(ascending=False,
                                                                             by='log2FC').index.tolist()


    # covid_logfc_all_targets_df = covid_logfc_all_targets_df
    # covid_logfc_all_targets_df.to_excel("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/pathway_DEGs/COVID_Jain_DEGs_Pathways_for_twodrugs.xlsx",index=False)
    # plt.figure(figsize=(70, 15))
    # ax = sns.boxplot(x="Target", y="log2FC",
    #                  data=covid_logfc_all_targets_df,order=path_groups_sorted)

    # plt.figure(figsize=(8, 5))
    plt.figure(figsize=(5, 8))

    sns.set_context( context="paper", font_scale=1.8)

    flierprops = dict(marker='.',markerfacecolor='1', markersize=2,
                      linestyle='none')

    # ax = sns.boxplot(y="Target", x="log2FC", width=.6,
    #                  data=covid_logfc_all_targets_df,flierprops=flierprops)
    ax = sns.boxplot(x="Target", y="log2FC", width=.6,
                     data=covid_logfc_all_targets_df,flierprops=flierprops)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

    ########
    mybox = ax.artists[2]
    print(mybox)
    # Change the appearance of that box

    for i, box in enumerate(ax.artists):
        if i%2 == 0:
            case_col = 'black'
            box.set_edgecolor(case_col)
            box.set_facecolor('white')
            for j in range(6 * i, 6 * (i + 1)):
                ax.lines[j].set_color(case_col)
                ax.lines[j].set_mfc(case_col)
                ax.lines[j].set_mec(case_col)


        else:
            ctl_col = 'grey'
            box.set_edgecolor(ctl_col)
            box.set_facecolor('white')
            for j in range(6 * i, 6 * (i + 1)):
                ax.lines[j].set_color(ctl_col)
                ax.lines[j].set_mfc(ctl_col)
                ax.lines[j].set_mec(ctl_col)

    #########################
    plt.tight_layout()

    # plt.savefig("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/one_neighbor/{}_drugs_target_neighbor_FC_boxplot_MS_Arun_control.pdf".format(drug))
    plt.savefig( save_addr)
    plt.show()


    #############################################################

def make_oneneighbor_randomCtrl_matrix(drug):

    all_shortest_results = dh.load_obj(
        "/Users/woochanghwang/PycharmProjects/MTIProject/COVID-19/result/network/SP_based/COVID_All_Structure_All_Shortest_Paths_24hr")

    all_shortest_paths = []
    for band, shortest_path in all_shortest_results.items():
        list_paths = list(shortest_path)
        list_paths = [list(path) for path in list_paths]
        all_shortest_paths += list_paths
    #
    all_shortest_paths = list(map(tuple, all_shortest_paths))


    drug_targets = get_targets_a_drug(drug)
    print(drug_targets)

    covid_degs_df = get_COVID_Arun_DEGs()
    covid_degs_df = covid_degs_df.dropna()

    covid_omics_logfc=covid_degs_df['log2FC'].tolist()
    covid_degs_path_df = pd.DataFrame(columns=['Gene', 'log2FC', 'COVID'])

    # #############################
    # # Check DEP
    # ##############################
    # dep_genes = pd.read_csv("/Users/woochanghwang/PycharmProjects/MTIProject/COVID-19/data/DEP/christian_ms_24h.p0.05.fc0.5.txt",names=['Gene'])['Gene'].tolist()

    ##############
    ## Key pathway
    ##############
    if drug == "gallopamil":
        reactom_addr = "/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/data/reactom_GPCR_downstream.txt"
    elif drug == "triamcinolone":
        reactom_addr = "/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/data/reactom_metabolism_of_lipids.txt"
    reactom_enriched_path = pd.read_csv(reactom_addr)['Gene'].tolist()

    ##############################
    ## per target
    #############################
    target_neighbor_omics_dict = dict()
    random_set_dict = dict()
    for target in drug_targets:
        print("target:", target)
        paths_including_atarget = []
        for path in all_shortest_paths:
            path_containing_target = check_target_exist(path, target)
            # print(path_containing_target)
            if path_containing_target[0] != 'NA':
                paths_including_atarget.append(path_containing_target)


        paths_including_targets_edges = []
        for path in paths_including_atarget:
            paths_including_targets_edges += path


        graph_including_targets = nx.Graph(paths_including_targets_edges)
        # graph_including_targets = nx.DiGraph(paths_including_targets_edges)


        ################
        ## one-neighbor
        ################
        targets_neighbor = list(graph_including_targets.neighbors(target))
        targets_neighbor.append(target)
        # print(targets_neighbor)
        #
        # covid_target_neighbor = covid_degs_df[covid_degs_df['Gene'].isin(targets_neighbor)]
        # covid_target_neighbor["Target"] = target
        # covid_degs_path_df = covid_degs_path_df.append(covid_target_neighbor, ignore_index=True)
        #
        # target_neighbor_logfc = covid_target_neighbor['log2FC'].tolist()
        #
        # target_neighbor_omics_dict[target] = target_neighbor_logfc
        #
        # print("target:", len(targets_neighbor))
        #
        # target_len = len(targets_neighbor)
        # random_num = 100
        # random_set = []
        # for i in range(random_num):
        #     random_targets = random.sample(covid_omics_logfc,target_len)
        #     random_set.append(random_targets)
        # target_neighbor_omics_dict['{}_Control'.format(target)] = random_set
        # random_set_dict[target] = random_set
        #
        # print(random_set)

        #################################
        ## istead of target neiggbor => common genes with DEP
        ## DEP --> enriched pathway(Gallopamil, GPCR downstream)
        #################################

        # nodes = graph_including_targets.nodes()
        # keyPath_in_paths = list(set(reactom_GPCR_downstream)&(set(nodes)))
        keyPath_in_paths = list(set(reactom_enriched_path) & (set(targets_neighbor)))
        keyPath_in_paths.append(target)

        covid_target_neighbor = covid_degs_df[covid_degs_df['Gene'].isin(keyPath_in_paths)]
        covid_target_neighbor["Target"] = target
        covid_degs_path_df = covid_degs_path_df.append(covid_target_neighbor, ignore_index=True)

        target_neighbor_logfc = covid_target_neighbor['log2FC'].tolist()

        target_neighbor_omics_dict[target] = target_neighbor_logfc

        print("target:", len(keyPath_in_paths))

        target_len = len(keyPath_in_paths)
        random_num = 100
        random_set = []
        for i in range(random_num):
            random_targets = random.sample(covid_omics_logfc,target_len)
            random_set.append(random_targets)
        target_neighbor_omics_dict['{}_Control'.format(target)] = random_set
        random_set_dict[target] = random_set

        print(random_set)
    print(covid_degs_path_df.head())
    covid_degs_path_df.to_excel("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/one_neighbor/{}_target_oneneighbor_Reactome_covid_Arun_moderate.xlsx".format(drug))

    return target_neighbor_omics_dict,random_set_dict

def make_target_RWR_result_randomCtrl_matrix(drugs):

    drug_targets = get_target_RWR_result_drugs(drugs)
    print(drug_targets)

    covid_degs_df = get_COVID_Arun_DEGs()
    covid_degs_df = covid_degs_df.dropna()

    ###########################
    # between -5 ~ 5
    ##########################
    covid_degs_df = covid_degs_df[(covid_degs_df['log2FC']<5)&(covid_degs_df['log2FC']>-5)]
    ###########################

    covid_omics_logfc=covid_degs_df['log2FC'].tolist()
    covid_degs_path_df = pd.DataFrame(columns=['Gene', 'log2FC', 'COVID'])

    target_extended_omics_dict = dict()
    random_set_dict = dict()



    for drug in drugs:
        extended_tragets = drug_targets[drug_targets['Drug']==drug]["Extened_targets"].tolist()

        extended_tragets = extended_tragets[0].split(',')
        print(extended_tragets)
        covid_target_neighbor = covid_degs_df[covid_degs_df['Gene'].isin(extended_tragets)]
        covid_target_neighbor["Target"] = drug
        covid_degs_path_df = covid_degs_path_df.append(covid_target_neighbor, ignore_index=True)

        target_neighbor_logfc = covid_target_neighbor['log2FC'].tolist()

        target_extended_omics_dict[drug] = target_neighbor_logfc

        print("target:", len(extended_tragets))

        target_len = len(extended_tragets)
        random_num = 100
        random_set = []
        for i in range(random_num):
            random_targets = random.sample(covid_omics_logfc,target_len)
            random_set.append(random_targets)
        target_extended_omics_dict['{}_Control'.format(drug)] = random_set
        random_set_dict[drug] = random_set

        print(random_set)
    print(covid_degs_path_df.head())


    # covid_degs_path_df.to_excel("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/RWR_extended/{}_target_oneneighbor_Reactome_covid_Arun_moderate.xlsx".format(drug))

    print(target_extended_omics_dict)
    return target_extended_omics_dict,random_set_dict

def make_kinase_randomCtrl_matrix(genelist):

    # drug_targets = get_target_RWR_result_drugs(drugs)
    # print(drug_targets)
    drug = "kinase"
    covid_degs_df = get_COVID_Arun_DEGs()
    covid_degs_df = covid_degs_df.dropna()

    covid_omics_logfc=covid_degs_df['log2FC'].tolist()
    covid_degs_path_df = pd.DataFrame(columns=['Gene', 'log2FC', 'COVID'])

    target_extended_omics_dict = dict()
    random_set_dict = dict()

    extended_tragets = genelist

    print(extended_tragets)
    covid_target_neighbor = covid_degs_df[covid_degs_df['Gene'].isin(extended_tragets)]
    covid_target_neighbor["Target"] = drug
    covid_degs_path_df = covid_degs_path_df.append(covid_target_neighbor, ignore_index=True)

    target_neighbor_logfc = covid_target_neighbor['log2FC'].tolist()

    target_extended_omics_dict[drug] = target_neighbor_logfc

    print("target:", len(extended_tragets))

    target_len = len(extended_tragets)
    random_num = 100
    random_set = []
    for i in range(random_num):
        random_targets = random.sample(covid_omics_logfc, target_len)
        random_set.append(random_targets)
    target_extended_omics_dict['{}_Control'.format(drug)] = random_set
    random_set_dict[drug] = random_set

    print(random_set)


    print(covid_degs_path_df.head())

    covid_degs_path_df.to_excel("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/kinase/{}_0to24_target_covid_Arun_moderate.xlsx".format(drug))

    print(target_extended_omics_dict)
    return target_extended_omics_dict,random_set_dict


def make_oneneighbor_SOM_path_randomCtrl_matrix(drug,path_name,path_genes):

    all_shortest_results = dh.load_obj(
        "/Users/woochanghwang/PycharmProjects/MTIProject/COVID-19/result/network/SP_based/COVID_All_Structure_All_Shortest_Paths_24hr")

    all_shortest_paths = []
    for band, shortest_path in all_shortest_results.items():
        list_paths = list(shortest_path)
        list_paths = [list(path) for path in list_paths]
        all_shortest_paths += list_paths
    #
    all_shortest_paths = list(map(tuple, all_shortest_paths))


    drug_targets = get_targets_a_drug(drug)
    print(drug_targets)

    covid_degs_df = get_COVID_Arun_DEGs()
    covid_degs_df = covid_degs_df.dropna()

    ###########################
    # between -5 ~ 5
    ##########################
    covid_degs_df = covid_degs_df[(covid_degs_df['log2FC']<5)&(covid_degs_df['log2FC']>-5)]
    ###########################

    covid_omics_logfc=covid_degs_df['log2FC'].tolist()
    covid_degs_path_df = pd.DataFrame(columns=['Gene', 'log2FC', 'COVID'])

    # #############################
    # # Check DEP
    # ##############################
    # dep_genes = pd.read_csv("/Users/woochanghwang/PycharmProjects/MTIProject/COVID-19/data/DEP/christian_ms_24h.p0.05.fc0.5.txt",names=['Gene'])['Gene'].tolist()

    ##############
    ## Key pathway
    ##############
    # if drug == "gallopamil":
    #     reactom_addr = "/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/data/reactom_GPCR_downstream.txt"
    # elif drug == "triamcinolone":
    #     reactom_addr = "/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/data/reactom_metabolism_of_lipids.txt"
    # reactom_enriched_path = pd.read_csv(reactom_addr)['Gene'].tolist()

    reactom_enriched_path = path_genes
    ##############################
    ## per target
    #############################
    target_neighbor_omics_dict = dict()
    random_set_dict = dict()
    for target in drug_targets:
        print("target:", target)
        paths_including_atarget = []
        for path in all_shortest_paths:
            path_containing_target = check_target_exist(path, target)
            # print(path_containing_target)
            if path_containing_target[0] != 'NA':
                paths_including_atarget.append(path_containing_target)


        paths_including_targets_edges = []
        for path in paths_including_atarget:
            paths_including_targets_edges += path


        graph_including_targets = nx.Graph(paths_including_targets_edges)
        # graph_including_targets = nx.DiGraph(paths_including_targets_edges)


        ################
        ## one-neighbor
        ################
        targets_neighbor = list(graph_including_targets.neighbors(target))
        targets_neighbor.append(target)
        # print(targets_neighbor)
        #
        # covid_target_neighbor = covid_degs_df[covid_degs_df['Gene'].isin(targets_neighbor)]
        # covid_target_neighbor["Target"] = target
        # covid_degs_path_df = covid_degs_path_df.append(covid_target_neighbor, ignore_index=True)
        #
        # target_neighbor_logfc = covid_target_neighbor['log2FC'].tolist()
        #
        # target_neighbor_omics_dict[target] = target_neighbor_logfc
        #
        # print("target:", len(targets_neighbor))
        #
        # target_len = len(targets_neighbor)
        # random_num = 100
        # random_set = []
        # for i in range(random_num):
        #     random_targets = random.sample(covid_omics_logfc,target_len)
        #     random_set.append(random_targets)
        # target_neighbor_omics_dict['{}_Control'.format(target)] = random_set
        # random_set_dict[target] = random_set
        #
        # print(random_set)

        #################################
        ## istead of target neiggbor => common genes with DEP
        ## DEP --> enriched pathway(Gallopamil, GPCR downstream)
        #################################

        # nodes = graph_including_targets.nodes()
        # keyPath_in_paths = list(set(reactom_GPCR_downstream)&(set(nodes)))
        keyPath_in_paths = list(set(reactom_enriched_path) & (set(targets_neighbor)))
        keyPath_in_paths.append(target)

        covid_target_neighbor = covid_degs_df[covid_degs_df['Gene'].isin(keyPath_in_paths)]
        covid_target_neighbor["Target"] = target
        covid_degs_path_df = covid_degs_path_df.append(covid_target_neighbor, ignore_index=True)

        target_neighbor_logfc = covid_target_neighbor['log2FC'].tolist()

        target_neighbor_omics_dict[target] = target_neighbor_logfc

        print("target:", len(keyPath_in_paths))

        target_len = len(keyPath_in_paths)
        random_num = 100
        random_set = []
        for i in range(random_num):
            random_targets = random.sample(covid_omics_logfc,target_len)
            random_set.append(random_targets)
        target_neighbor_omics_dict['{}_Control'.format(target)] = random_set
        random_set_dict[target] = random_set

        print(random_set)
    print(covid_degs_path_df.head())
    covid_degs_path_df.to_excel("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/one_neighbor/{}_target_oneneighbor_SOM_{}_covid_Arun_moderate.xlsx".format(drug,path_name))

    return target_neighbor_omics_dict,random_set_dict

def calc_significancy_boxplot(drug, covid_logfc_all_targets_df, save_addr):
    from scipy import stats
    drug_targets = get_targets_a_drug(drug)
    ####################
    # For SOM
    ##############
    # drug_targets = ['gallopamil', 'triamcinolone']
    ########

    significancy_result_dict = dict()
    for target in drug_targets:
        print(target)
        target_logfc_list = covid_logfc_all_targets_df[covid_logfc_all_targets_df['Target'] == target]['log2FC'].tolist()
        target_control_logfc_list = covid_logfc_all_targets_df[covid_logfc_all_targets_df['Target'] == target+"_Control"]['log2FC'].tolist()

        # t,p = stats.ttest_ind(target_logfc_list, target_control_logfc_list)
        # print("Ttest Target, t, p:", target, t, p)
        # ks,ks_p = stats.ks_2samp(target_logfc_list,target_control_logfc_list)
        # Mann - Whitney - Wilcoxon test
        t,p = stats.mannwhitneyu(target_logfc_list,target_control_logfc_list)
        print("mann whitneyu Target, t, p:", target, t,p)
        significancy_result_dict[target]= [t,p]

    stat_result_df = pd.DataFrame.from_dict(significancy_result_dict, orient="index", columns=['mannwhitneyu','p-value'])

    print(stat_result_df)
    # stat_result_df.to_excel("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/one_neighbor/{}_drugs_target_neighbor_control_mannwhitney.xlsx".format(drug))
    stat_result_df.to_excel( save_addr )

def get_paths_with_targets_graphs(drug):
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

    ##############################
    ## per target
    #############################
    nodes_from_SIP = dict()
    graph_for_target_dict = dict()

    for target in drug_targets:
        print("target:", target)
        paths_including_atarget = []
        for path in all_shortest_paths:
            path_containing_target = check_target_exist(path, target)
            # print(path_containing_target)
            if path_containing_target[0] != 'NA':
                paths_including_atarget.append(path_containing_target)

        print(len(paths_including_atarget))

        paths_including_targets_edges = []
        for path in paths_including_atarget:
            paths_including_targets_edges += path


        # graph_including_targets = nx.Graph(paths_including_targets_edges)
        graph_including_targets = nx.DiGraph(paths_including_targets_edges)

        nodes_in_graph_with_targets = graph_including_targets.nodes()

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
        graph_for_target_dict[target] = graph_including_targets

        print("target, random:", len(targets_neighbor))

    return nodes_from_SIP, graph_for_target_dict

def one_neighbor_analysis():
    covid_degs_df = get_COVID_Arun_DEGs()

    drug_names = ['gallopamil', 'triamcinolone']
    drug = drug_names[0]
    nodes_from_SIP,graph_including_targets = get_paths_with_targets_graphs(drug)

    covid_degs_path_gall_df = pd.DataFrame(columns=['Gene', 'log2FC', 'COVID'])

    ############################################################
    # Tartet neithbot vs All genes  , Gallopamil
    #############################################################
    for path, genes in nodes_from_SIP.items():
        # print(path)
        covid_path_df = covid_degs_df[covid_degs_df['Gene'].isin(genes)]
        covid_path_df['Target'] = path
        # print(covid_path_df)

        covid_degs_path_gall_df = covid_degs_path_gall_df.append(covid_path_df, ignore_index=True)
        # print(covid_degs_path_df)

    gall_between_df = pd.read_excel("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/one_neighbor/gallapamil_drugs_target_neighbor_between_centrality.xlsx",
                                    names=['Gene','BC'])
    gall_nodes = graph_including_targets.nodes()
    print(len(gall_nodes))


    # ####################################
    # # drug: Triamcinolone
    # ####################################
    # drug = drug_names[1]
    # nodes_from_SIP, graph_including_targets = get_paths_with_targets_N_control(drug)
    #
    # covid_degs_path_tria_df = pd.DataFrame(columns=['Gene', 'log2FC', 'COVID'])
    #
    # ############################################################
    # # Tartet neithbot vs All genes  , Gallopamil
    # #############################################################
    # for path, genes in nodes_from_SIP.items():
    #     # print(path)
    #     covid_path_df = covid_degs_df[covid_degs_df['Gene'].isin(genes)]
    #     covid_path_df['Target'] = path
    #     # print(covid_path_df)
    #
    #     covid_degs_path_tria_df = covid_degs_path_tria_df.append(covid_path_df, ignore_index=True)
    #     # print(covid_degs_path_df)

    # trim_between = nx.betweenness_centrality(graph_including_targets)

    # print(gall_between)
    # print(trim_between)

    # gall_between_df = pd.DataFrame.from_dict(gall_between,orient="index")
    # trim_between_df = pd.DataFrame.from_dict(trim_between,orient="index")

    # gall_between_df.to_excel("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/one_neighbor/gallapamil_drugs_target_neighbor_between_centrality.xlsx")
    # trim_between_df.to_excel("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/one_neighbor/triamcinolone_drugs_target_neighbor_between_centrality.xlsx")


def draw_one_neighbor_boxplot():
    drug_names = ['gallopamil', 'triamcinolone']
    drug = drug_names[1]
    # # nodes_from_SIP = get_paths_with_targets_N_control(drug)
    # # covid_degs_df = get_COVID_Arun_DEGs()
    #
    # one_neighbor_analysis()
    ##################
    ## v1 , one-neighbor & a enriched path
    ######################
    # target_neighbor_omics_dict,random_set_dict= make_oneneighbor_randomCtrl_matrix(drug)
    #
    # covid_logfc_all_targets_df = make_boxplot_input(target_neighbor_omics_dict)
    # # covid_logfc_all_targets_df.to_excel("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/one_neighbor/{}_drugs_target_neighbor_control_FC_boxplot_input.xlsx".format(drug))
    # covid_logfc_all_targets_df.to_excel(
    #     "/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/one_neighbor/{}_drugs_target_Reactome_control_FC_boxplot_input.xlsx".format(
    #         drug))
    # draw_boxplot(drug,covid_logfc_all_targets_df)
    #
    # calc_significancy_boxplot(drug,covid_logfc_all_targets_df)

    ##################
    ## v2, one-neighbor & SOM path
    ##################
    reactom_name_dict, reactom_id_dict = get_reactom_paths()

    # two_drugs_paths = pd.read_csv("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/SOM/tria_gal_pathways.txt", sep='\t',names=['path'])['path'].tolist()

    # control_paths = \
    # pd.read_csv("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/SOM/dex_pathways.txt", sep='\t',
    #             names=['path'])['path'].tolist()

    two_drugs_paths = ["HDL remodeling"]
    for path in two_drugs_paths:
    # for path in control_paths:
        path = path.strip()
        print(reactom_name_dict[path])
        path_genes = reactom_name_dict[path]
        target_neighbor_omics_dict, random_set_dict = make_oneneighbor_SOM_path_randomCtrl_matrix(drug, path,
                                                                                                  path_genes)

        covid_logfc_all_targets_df = make_boxplot_input(target_neighbor_omics_dict)

        path = path.replace(",", ".")
        # covid_logfc_all_targets_addr = "/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/one_neighbor/SOM_pathways/v3/{}_drugs_target_Reactome_{}_control_FC_boxplot_input.xlsx".format(
        #     drug, path)
        # boxplot_addr = "/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/one_neighbor/SOM_pathways/v3/{}_drugs_target_oneneighbor_reactome_{}_FC_boxplot_MS_Arun_control.pdf".format(
        #     drug, path)
        # significancy_boxplot_addr = "/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/one_neighbor/SOM_pathways/v3/{}_drugs_target_oneneighbor_Reactome_{}_control_mannwhitneyu.xlsx".format(
        #     drug, path)

        covid_logfc_all_targets_addr = "/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/one_neighbor/SOM_pathways/v3/{}_drugs_target_Reactome_{}_control_FC_boxplot_input_vertical.xlsx".format(
            drug, path)
        boxplot_addr = "/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/one_neighbor/SOM_pathways/v3/{}_drugs_target_oneneighbor_reactome_{}_FC_boxplot_MS_Arun_control_vertical_grey.pdf".format(
            drug, path)
        significancy_boxplot_addr = "/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/one_neighbor/SOM_pathways/v3/{}_drugs_target_oneneighbor_Reactome_{}_control_mannwhitneyu_vertical.xlsx".format(
            drug, path)

        # covid_logfc_all_targets_df.to_excel("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/one_neighbor/{}_drugs_target_neighbor_control_FC_boxplot_input.xlsx".format(drug))
        # covid_logfc_all_targets_df.to_excel(covid_logfc_all_targets_addr)
        draw_boxplot(drug, covid_logfc_all_targets_df, boxplot_addr)

        # calc_significancy_boxplot(drug, covid_logfc_all_targets_df, significancy_boxplot_addr)

def draw_SOM_input_boxplot():
    drug_names = ['gallopamil', 'triamcinolone']

    target_neighbor_omics_dict, random_set_dict = make_target_RWR_result_randomCtrl_matrix(drug_names)


    covid_logfc_all_targets_df = make_boxplot_input(target_neighbor_omics_dict)


    covid_logfc_all_targets_addr = "/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/RWR_extended/twodrugs_extended_targets_control_FC_boxplot_input.xlsx"
    boxplot_addr = "/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/RWR_extended/twodrugs_extended_targets_FC_boxplot_MS_Arun_control.pdf"
    significancy_boxplot_addr = "/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/RWR_extended/twodrugs_extended_targets_control_mannwhitney.xlsx"

    drug = "twodrugs"
    # covid_logfc_all_targets_df.to_excel(covid_logfc_all_targets_addr)
    draw_boxplot(drug, covid_logfc_all_targets_df, boxplot_addr)

    calc_significancy_boxplot(drug, covid_logfc_all_targets_df, significancy_boxplot_addr)


def draw_kinase_boxplot():
    # drug_names = ['gallopamil', 'triamcinolone']

    krogan_kinase_list = pd.read_csv(
        "/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/kinase/from_Meabh/0hto24hchangedkinases.txt",
        names=['Gene'])['Gene'].tolist()

    target_neighbor_omics_dict, random_set_dict = make_kinase_randomCtrl_matrix(krogan_kinase_list)

    covid_logfc_all_targets_df = make_boxplot_input(target_neighbor_omics_dict)

    covid_logfc_all_targets_addr = "/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/kinase/kinase_0to24_control_FC_boxplot_input.xlsx"
    boxplot_addr = "/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/kinase/kinase_0to24_FC_boxplot_MS_Arun_control.pdf"
    significancy_boxplot_addr = "/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/kinase/kinase_0to24_control_ttest.xlsx"

    drug = "kinase"
    covid_logfc_all_targets_df.to_excel(covid_logfc_all_targets_addr)
    draw_boxplot(drug, covid_logfc_all_targets_df, boxplot_addr)

    # calc_significancy_boxplot(drug, covid_logfc_all_targets_df, significancy_boxplot_addr)


def make_expression_patterns_geneSet():

    krogan_kinase_list = pd.read_csv("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/kinase/from_Meabh/Kinases_viral_infection_highly_changed_list.txt",
                                     names=['Gene'])['Gene'].tolist()
    covid_degs_df = get_COVID_Arun_DEGs()
    covid_kinase_df = covid_degs_df[covid_degs_df['Gene'].isin(krogan_kinase_list)]

    covid_kinase_df.to_excel("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/kinase/covid_kinase_krogan_moderate.xlsx",index=False)



def main():
    draw_one_neighbor_boxplot()

    # draw_SOM_input_boxplot()
    #############
    ## Check kinase(Meabh) expression pattern
    ################
    # make_expression_patterns_geneSet()
    # draw_kinase_boxplot()




def get_reactom_paths():
    reactome_addr= "/Users/woochanghwang/PycharmProjects/MTIProject/General/data/Reactome/Current/ReactomePathways.gmt"
    reactom_name_dict = dict()
    reactom_id_dict = dict()

    with open(reactome_addr) as reactom_f:
        reactom_data = [x.strip().split('\t') for x in reactom_f.readlines()]

    for path in reactom_data:
        reactom_name_dict[path[0]] = path[2:]
        reactom_id_dict[path[1]] = path[2:]
    # print(reactom_data[:2])
    # print(reactom_name_dict)

    ####################
    # Test for HfH
    ####################
    two_drugs_paths = pd.read_csv("/Users/woochanghwang/PycharmProjects/MTIProject/Hackerton/result/SOM/tria_gal_pathways.txt", sep='\t',names=['path'])['path'].tolist()
    for path in two_drugs_paths:
        path = path.strip()
        print(reactom_name_dict[path])


    return reactom_name_dict, reactom_id_dict

if __name__ == '__main__':
    main()
    # get_reactom_paths()