# Created by woochanghwang at 29/01/2021
import pandas as pd

def main():
    threshold = 800
    hfh_result_df = pd.read_excel(
        "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance_druginfo_ATC_{}.xlsx".format(
            800), sheet_name="filtered")
    print(hfh_result_df)

    control_drugs = list(set(hfh_result_df['Control_drugname'].tolist()))
    similar_drugs = list(set(hfh_result_df['Similar_drugname'].tolist()))
    filtered_drugs = list(set(control_drugs+similar_drugs))

    # print(len(set(control_drugs)&set(similar_drugs)))   # There is no common drugs between control & similar.
    print("control:", len(control_drugs))
    print("similar:", len(similar_drugs))
    print("all", len(filtered_drugs))
    filtered_drug_rwr_dict = dict()

    hfh_drug_targets = pd.read_excel(
        "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Similar_drugs/Control_and_identified_name_and_drugbankID_network_similarity_distance_druginfo_ATC_{}.xlsx".format(
            800), sheet_name="filter_drugs")

    print(hfh_drug_targets)
    mean_N_target = 14
    RWR_result_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Drug_Pagerank/{}_{}/COVID_HfH_{}_Pagerank.xlsx"
    #######
    # Control
    #######
    type = "Control"
    for drug in control_drugs:
        n_drug_target = int(hfh_drug_targets[hfh_drug_targets['Drug'] == drug]['N_Targets_SIP'])

        rwr_result_df = pd.read_excel(RWR_result_addr.format(type,threshold,drug))

        rwr_result_df = rwr_result_df.rename(columns={
            "Unnamed: 0" : "Genes"
        })
        rwr_result_df= rwr_result_df.sort_values(by='PageRank',ascending=False)
        # print(rwr_result_df.head(3))
        n_extended_target = max(mean_N_target,n_drug_target)
        # print("mean, ndrug, ext",mean_N_target,n_drug_target,n_extended_target)
        extended_targets = rwr_result_df['Genes'].tolist()[:n_extended_target]    # Max number of targets : 125 --> mean:14
        print(len(extended_targets), drug)
        filtered_drug_rwr_dict[drug] = ','.join(extended_targets)

    type = 'Similar'
    for drug in similar_drugs:
        n_drug_target= int(hfh_drug_targets[hfh_drug_targets['Drug']==drug]['N_Targets_SIP'])
        # print(n_drug_target)
        rwr_result_df = pd.read_excel(RWR_result_addr.format(type,threshold,drug))

        rwr_result_df = rwr_result_df.rename(columns={
            "Unnamed: 0" : "Genes"
        })
        rwr_result_df= rwr_result_df.sort_values(by='PageRank',ascending=False)
        # print(rwr_result_df.head(3))
        n_extended_target = max(mean_N_target, n_drug_target)
        # print("mean, ndrug, ext", mean_N_target, n_drug_target, n_extended_target)
        extended_targets = rwr_result_df['Genes'].tolist()[:n_extended_target]
        print(len(extended_targets), drug)
        # print(extended_targets[:5])

        filtered_drug_rwr_dict[drug] = ','.join(extended_targets)


    filtered_drugs_df = pd.DataFrame.from_dict(filtered_drug_rwr_dict, orient="index",
                                               columns=['Extened_targets'])
    filtered_drugs_df = filtered_drugs_df.reset_index()
    filtered_drugs_df  = filtered_drugs_df.rename(columns={
        "index" : "Drug"
    })
    print(filtered_drugs_df)

    filtered_drugs_df.to_csv("/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Drug_Pagerank/filtered_800/filtered_drugs_rwr_extened_targets_mean_{}.csv".format(threshold),index=False)
    filtered_drugs_df.to_excel(
        "/Users/woochanghwang/PycharmProjects/LifeArc/Hackerton/result/Drug_Pagerank/filtered_800/filtered_drugs_rwr_extened_targets_mean_{}.xlsx".format(
            threshold), index=False)

if __name__ == '__main__':
    main()
