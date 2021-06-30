# Created by woochanghwang at 06/01/2021
import pandas as pd

def main_labellers():
    similar_drugs_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance.xlsx"
    drug_info_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/Drugbank/v5.1.6/drugbank_routes_labellers.tsv"

    drug_info_df = pd.read_csv(drug_info_addr,sep='\t')
    similar_drugs_df = pd.read_excel(similar_drugs_addr)

    # print(similar_drugs_df)
    print(drug_info_df)

    similar_drugs_info_df = pd.merge(left=similar_drugs_df,
                                     right=drug_info_df,
                                     how="left",
                                     left_on="Similar_id",
                                     right_on="drugbank_id")
    print(similar_drugs_info_df)

    similar_drugs_info_df.to_excel("/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance_druginfo_route.xlsx")

def main():
    threshold = 900
    similar_drugs_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance_{}.xlsx".format(threshold)
    drug_info_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/Drugbank/v5.1.6/drugbank.tsv"

    drug_info_df = pd.read_csv(drug_info_addr, sep='\t')
    similar_drugs_df = pd.read_excel(similar_drugs_addr)

    # print(similar_drugs_df)
    drug_info_df = drug_info_df[['drugbank_id','groups','atc_codes']]
    print(list(drug_info_df))

    similar_drugs_info_df = pd.merge(left=similar_drugs_df,
                                     right=drug_info_df,
                                     how="left",
                                     left_on="Similar_id",
                                     right_on="drugbank_id")
    print(similar_drugs_info_df)

    similar_drugs_info_df.to_excel(
        "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance_druginfo_ATC_{}.xlsx".format(threshold))

if __name__ == '__main__':
    main()
