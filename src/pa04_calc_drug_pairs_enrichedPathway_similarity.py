# Created by woochanghwang at 20/01/2021
import pandas as pd

def jaccard_index (L1,L2):
    similarity = len(set(L1).intersection(set(L2))) / len(set(L1).union(set(L2)))
    return similarity

def overlap(L1,L2):
    similarity = len(set(L1).intersection(set(L2))) / min(len(set(L1)), len(set(L2)))
    return similarity

def calc_pa_similarity(control,similar, way):
    pa_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Drug_pathways_meanTarget_800_filtered/enriched_reactome_pathways_{}.csv"
    control_df = pd.read_csv(pa_addr.format(control))
    similar_df = pd.read_csv(pa_addr.format(similar))
    print(list(control_df))
    control_pa  = control_df['term_name'].tolist()
    similar_pa = similar_df['term_name'].tolist()
    if way ==  "jaccard":
        return jaccard_index(control_pa,similar_pa)
    elif way == "overlap":
        return overlap(control_pa,similar_pa)





    return 1
def main():
    threshold = 800
    hfh_result_df = pd.read_excel(
        "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance_druginfo_ATC_{}.xlsx".format(
            800), sheet_name="filtered")
    print(hfh_result_df)

    # hfh_drug_df = hfh_result_df[['Control_drugname','Similar_drugname']]
    # print(hfh_drug_df)
    # # hfh_result_df = hfh_result_df.fillna('NA')
    # for index, row in hfh_drug_df.iterrows():
    #     contol = row['Control_drugname']
    #     similar = row['Similar_drugname']

    hfh_result_df['PA_similar_Jaccard'] = hfh_result_df.apply(lambda row:
                                                      calc_pa_similarity(
                                                          row['Control_drugname'],
                                                          row['Similar_drugname'],
                                                          way="jaccard"
                                                      ),
                                                      axis=1)

    hfh_result_df['PA_similar_overlap'] = hfh_result_df.apply(lambda row:
                                                              calc_pa_similarity(
                                                                  row['Control_drugname'],
                                                                  row['Similar_drugname'],
                                                                  way="overlap"
                                                              ),
                                                              axis=1)

    hfh_result_df.to_excel("/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance_druginfo_ATC_PA_{}.xlsx".format(
            800),index=False)




if __name__ == '__main__':
    main()
