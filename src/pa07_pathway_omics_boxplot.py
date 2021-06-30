# Created by woochanghwang at 16/02/2021
'''
Draw boxplot for drug enriched pathways
'''

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def read_reactom(reactom_addr):
    with open(reactom_addr) as reactom_f:
        reactom_r = [x.strip().split('\t') for x in reactom_f.readlines()]

    # print(reactom_r[:3])
    reactom_name_dict = dict()
    for path in reactom_r:
        reactom_name_dict[path[1]] = path[2:]

    return reactom_name_dict

def get_COVID_Jain_DEGs():
    jain_DEGs_moderate = pd.read_excel("/Users/woochanghwang/PycharmProjects/LifeArc/COVID19_vaccine/data/Jain/Jain_CSBJ_supple2_v2.xlsx",sheet_name="Control_moderate")
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
    arun_DEGs_moderate= pd.read_excel("/Users/woochanghwang/PycharmProjects/LifeArc/COVID19_vaccine/data/Arunachalam/Arunachalam_DEGs_M_S.xlsx",
                                      sheet_name="Moderate")
    arun_DEGs_moderate = arun_DEGs_moderate[['Symbol','log2FoldChange_MH']]
    arun_DEGs_moderate = arun_DEGs_moderate.rename(columns={ 'Symbol':'Gene',
                                                             'log2FoldChange_MH' : 'log2FC'})
    arun_DEGs_moderate['COVID'] ='Moderate'

    return arun_DEGs_moderate

def get_COVID_Overmyer_DEGs():

    DEGs_moderate= pd.read_excel("/Users/woochanghwang/PycharmProjects/LifeArc/COVID19_vaccine/data/Overmyer/Overmyer_DEGs_M_S.xlsx",
                                      sheet_name="Moderate")
    DEGs_moderate = DEGs_moderate[['Gene','log2FC']]
    # DEGs_moderate = DEGs_moderate.rename(columns={ 'Symbol':'Gene',
    #                                                          'log2FoldChange_MH' : 'log2FC'})
    DEGs_moderate['COVID'] ='Moderate'

    return DEGs_moderate

def get_COVID_xu_DEGs():
    DEGs_moderate = pd.read_excel("/Users/woochanghwang/PycharmProjects/LifeArc/COVID19_vaccine/data/Xu/Xu_DEGs_scRNA.xlsx")
    DEGs_moderate = DEGs_moderate[DEGs_moderate['compare']=='MH']
    DEGs_moderate = DEGs_moderate[['gene','avg_logFC']]
    DEGs_moderate = DEGs_moderate.rename(columns={ 'gene':'Gene',
                                                    'avg_logFC' : 'log2FC'})
    DEGs_moderate['COVID'] = 'Moderate'

    return DEGs_moderate

def get_two_drug_pathways():
    drug_names = ['gallopamil', 'triamcinolone']
    drug_addr_prefix = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Drug_pathways_meanTarget_800_filtered/enriched_reactome_pathways_{}.csv"

    drug_enriched_pathways = []
    for drug in drug_names:
        drug_addr = drug_addr_prefix.format(drug)
        drug_enriched_df = pd.read_csv(drug_addr)
        paths = drug_enriched_df['term_name'].tolist()
        drug_enriched_pathways += paths
        # print(drug_enriched_df)
        # print(list(drug_enriched_df))

    drug_enriched_pathways = list(set(drug_enriched_pathways))

    return drug_enriched_pathways

def get_all_enriched_pathways():
    f1_matrix_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Drug_pathways_meanTarget_800_filtered/F1/hfh_drug_to_low_level_pathways_0.01_matrix.csv"
    f1_matrix_df = pd.read_csv(f1_matrix_addr,index_col=0)
    f1_matrix_df = f1_matrix_df.T
    f1_matrix_df['Sum'] = f1_matrix_df.sum( axis=1)
    print(f1_matrix_df)

    pathways = f1_matrix_df.index.tolist()
    print(pathways)
    return pathways

def main():

    # drug_enriched_pathways = get_all_enriched_pathways()
    drug_enriched_pathways = get_two_drug_pathways()
    # print(len(drug_enriched_pathways))

    # print(drug_enriched_pathways[:10])

    reactom_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/Reactome/ReactomePathways_for_gProfiler.gmt"
    reactom_dict = read_reactom(reactom_addr)

    # covid_degs_df = get_COVID_Jain_DEGs()
    # covid_degs_df = get_COVID_Arun_DEGs()
    # covid_degs_df = get_COVID_Overmyer_DEGs()
    covid_degs_df = get_COVID_xu_DEGs()

    covid_degs_path_df = pd.DataFrame(columns=['Gene', 'log2FC', 'COVID'])

    for path in drug_enriched_pathways:
        print(path)
        path_genes = reactom_dict[path]
        print(len(reactom_dict[path]))
        covid_path_df = covid_degs_df[covid_degs_df['Gene'].isin(path_genes)]
        covid_path_df['Pathway'] = path
        # print(covid_path_df)

        covid_degs_path_df = covid_degs_path_df.append(covid_path_df, ignore_index=True)
        # print(covid_degs_path_df)

    print(covid_degs_path_df)
    path_groups_sorted = covid_degs_path_df.groupby(['Pathway']).median().sort_values(ascending=False,
                                                                             by='log2FC').index.tolist()
    covid_degs_path_df = covid_degs_path_df
    # covid_degs_path_df.to_excel("/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/pathway_DEGs/COVID_Jain_DEGs_Pathways_for_twodrugs.xlsx",index=False)
    plt.figure(figsize=(70, 15))
    ax = sns.boxplot(x="Pathway", y="log2FC", hue="COVID",
                     data=covid_degs_path_df,order=path_groups_sorted)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    plt.tight_layout()

    plt.savefig("/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/pathway_DEGs/two_drugs_all_path_deg_boxplot_M_Xu_median.pdf")
    plt.show()



def main_only_SOM():
    drug_enriched_pathways_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/SOM/tria_gal_pathways.txt"
    reactom_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/Reactome/ReactomePathways_for_gProfiler.gmt"

    reactom_dict = read_reactom(reactom_addr)

    # drug_enriched_pathways = pd.read_csv(drug_enriched_pathways_addr,sep='\t',names=['Path'])
    # print(drug_enriched_pathways)
    drug_enriched_pathways = pd.read_csv(drug_enriched_pathways_addr, sep='\t', names=['Path'])['Path'].tolist()
    drug_enriched_pathways = [x.strip() for x in drug_enriched_pathways]
    print(drug_enriched_pathways)

    covid_degs_df = get_COVID_Jain_DEGs()

    covid_degs_path_df = pd.DataFrame(columns=['Gene', 'log2FC', 'COVID'])

    for path in drug_enriched_pathways:
        print(path)
        path_genes = reactom_dict[path]
        print(len(reactom_dict[path]))
        covid_path_df = covid_degs_df[covid_degs_df['Gene'].isin(path_genes)]
        covid_path_df['Pathway'] = path
        # print(covid_path_df)

        covid_degs_path_df = covid_degs_path_df.append(covid_path_df, ignore_index=True)
        # print(covid_degs_path_df)

    print(covid_degs_path_df)
    covid_degs_path_df = covid_degs_path_df
    # covid_degs_path_df.to_excel("/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/pathway_DEGs/COVID_Jain_DEGs_Pathways_for_twodrugs.xlsx",index=False)
    plt.figure(figsize=(10, 15))
    ax = sns.boxplot(x="Pathway", y="log2FC", hue="COVID",
                     data=covid_degs_path_df)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    plt.tight_layout()

    plt.savefig("/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/pathway_DEGs/path_deg_boxplot_MMS.pdf")
    plt.show()


if __name__ == '__main__':
    main()
