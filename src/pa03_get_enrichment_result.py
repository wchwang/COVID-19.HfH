# Created by woochanghwang at 29/01/2021
import pandas as pd

def add_number_of_enrichment(filtered_drugs):


    enrichment_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Drug_pathways_800_filtered/enriched_reactome_pathways_{}.csv"

    drug_num_enrichment_dict = dict()

    for drug in filtered_drugs:
        drug_enrichment_df = pd.read_csv(enrichment_addr.format(drug))
        number_of_enrichment = len(drug_enrichment_df['term_name'].tolist())
        # print(list(drug_enrichment_df))
        drug_num_enrichment_dict[drug] = number_of_enrichment

    print(drug_num_enrichment_dict)
    drug_num_enrichment_df = pd.DataFrame.from_dict(drug_num_enrichment_dict,orient="index",
                                                    columns={'Num_Reactom'})
    drug_num_enrichment_df = drug_num_enrichment_df.reset_index()
    drug_num_enrichment_df = drug_num_enrichment_df.rename(
        columns={'index' : 'Drug'}
    )

    drug_num_enrichment_df = drug_num_enrichment_df.sort_values(by='Num_Reactom')
    print(list(drug_num_enrichment_df))

    drug_target_df = pd.read_excel("/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Drug_Pagerank/filtered_800/filtered_drugs_rwr_extened_targets_800.xlsx")
    # print(drug_target_df)
    print(list(drug_target_df))

    drug_target_pathways_df = pd.merge( left=drug_target_df,
                                        right = drug_num_enrichment_df,
                                        how="left",
                                        left_on="Drug",
                                        right_on="Drug")

    print(drug_target_pathways_df)

    drug_target_pathways_df.to_excel("/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Drug_Pagerank/filtered_800/filtered_drugs_rwr_extened_targets_num_pathway_800.xlsx")


def make_all_enriched_pathways(filtered_drugs):
    min_paths = 131
    all_enrichment_pathways = []
    enrichment_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Drug_pathways_meanTarget_800_filtered/enriched_reactome_pathways_{}.csv"

    for drug in filtered_drugs:
        drug_enrichment_df = pd.read_csv(enrichment_addr.format(drug))
        # print(list(drug_enrichment_df))
        drug_enrichment_df = drug_enrichment_df.sort_values(by='p_value')
        # enrichment_pathways = drug_enrichment_df['term_name'].tolist()[:min_paths]
        drug_enrichment_df = drug_enrichment_df[drug_enrichment_df['p_value']<0.01]
        enrichment_pathways = drug_enrichment_df['term_name'].tolist()
        print(drug, len(enrichment_pathways))
        all_enrichment_pathways += enrichment_pathways

    all_enrichment_pathways = list(set(all_enrichment_pathways))
    print(len(all_enrichment_pathways))

    return all_enrichment_pathways



def main():
    threshold = 800
    hfh_result_df = pd.read_excel(
        "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance_druginfo_ATC_{}.xlsx".format(
            800), sheet_name="filtered")
    # print(hfh_result_df)

    control_drugs = list(set(hfh_result_df['Control_drugname'].tolist()))
    similar_drugs = list(set(hfh_result_df['Similar_drugname'].tolist()))
    filtered_drugs = list(set(control_drugs+similar_drugs))

    ################################################
    ## 1. find min number of enrichment pathwayss
    ################################################
    # add_number_of_enrichment(filtered_drugs)    # Min number of enrichmeent is 131

    ################################################
    ## 2. make union enrichment pathways()
    ## FDR < 0.05
    ################################################
    all_enriched_pathways = make_all_enriched_pathways(filtered_drugs)


if __name__ == '__main__':
    main()
