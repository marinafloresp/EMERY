import os
import numpy as np
import streamlit as st
import pandas as pd
from common import saveTable,alt_ver_barplot,alt_scatterplot,alt_boxplot, ecdf,alt_line_chart

genes_path = os.path.join('..', 'Entrez_to_Ensg99.mapping_table.tsv')

"""
Class: Data: creates an object with the results data contained in the uploaded pickle file
Input: needs the pickle file data and the name of the file
"""
class Data():
    def __init__(self,pickle_data, name):
        self.file_name = name
        self.expr_df_filtered = pickle_data[list(pickle_data)[0]]
        self.expr_df_selected = pickle_data[list(pickle_data)[1]]
        self.drug_response = pickle_data[list(pickle_data)[2]]
        self.class_labels = pickle_data[list(pickle_data)[3]]
        self.cluster_labels = pickle_data[list(pickle_data)[4]]
        self.patient_df = pickle_data[list(pickle_data)[5]]
        self.importance_df = pickle_data[list(pickle_data)[8]]
        self.stacked_posterior = pickle_data[list(pickle_data)[6]]
        self.clusters_names = sorted(list(self.patient_df['cluster_number'].unique()))
        self.baskets_names = sorted(list(self.patient_df['tissues'].unique()))

    #Maps genes' ENSEMBL IDs to their names
    @staticmethod
    def IdGenes(df):
        features_EN = df.columns.tolist()
        genes_df = pd.read_csv(genes_path, sep='\t')
        matching_genes = genes_df[genes_df['ensg_v99'].isin(features_EN)].reset_index(drop=True)
        matching_genesEN = matching_genes['ensg_v99'].tolist()
        genes = []
        for col in features_EN:
            found_match = False
            for i in range(len(matching_genesEN)):
                if col == matching_genesEN[i]:
                    genes.append(matching_genes['gene_name_v99'][i])
                    found_match = True
                    break
            if not found_match:
                genes.append(col)
        return genes

    #Changes responsive feature to be categorical (Responsive vs Non-responsive) and maps ENSEMBL IDs for Genes' names
    def setData(self):
        self.patient_df["responsive"] = self.patient_df["responsive"].replace([0, 1], ['Non-responsive', 'Responsive'])
        selected_genes = Data.IdGenes(self.expr_df_selected)
        filtered_genes = Data.IdGenes(self.expr_df_filtered)
        self.expr_df_selected.columns = selected_genes
        self.importance_df.index = selected_genes
        self.expr_df_filtered.columns = filtered_genes

    #Displays information about the file in a table format
    def fileInfo(self):
        num_samples = len(self.expr_df_filtered.axes[0])
        num_feat_filt = len(self.expr_df_filtered.axes[1])
        num_feat_select = len(self.expr_df_selected.axes[1])
        clusters = len(self.clusters_names)
        tissues = len(self.baskets_names)
        drug = self.file_name.split('_')[2]
        info = {'File name: ': {'information': self.file_name}, 'Drug: ': {'information': drug},'Num. cell lines' : {'information':num_samples},
                'Num. transcripts (after filtering): ': {'information': num_feat_filt},
                'Num. transcripts (after feature selection): ': {'information': num_feat_select},
                'Num.clusters: ': {'information': clusters}, 'Num. tissues/baskets: ': {'information': tissues}}
        df = pd.DataFrame(data=info).T
        style = df.style.hide_index()
        style.hide_columns()
        st.dataframe(df, use_container_width=True)

    #Shows information about the number of samples per group in a bar plot. If RD (Response Data) is True, samples are grouped by response within groups.
    def countPlot(self, RD, title,feature,x_lab):
        self.patient_df[feature] = pd.Categorical(self.patient_df[feature])
        size = len(np.unique(self.patient_df[feature]))
        if RD:
            df_grouped = self.patient_df.groupby(["responsive", feature]).size().reset_index(name='Count')
            alt_ver_barplot(df_grouped, feature, 'Count', 2, x_lab, "Number of samples", "responsive", title, "NumSamples",
                            ["responsive", feature, 'Count'])
            st.caption("The x-axis shows the levels of the grouping chosen (clusters or baskets/tissues). The y-axis shows the number of samples."
                       " Within levels, the number of samples that are responsive (blue) or non-responsive (red) to treatment are shown.")
        else:
            df_grouped = self.patient_df.groupby([feature]).size().reset_index(name='Count')
            alt_ver_barplot(df_grouped, feature, 'Count', size, x_lab, "Number of samples", feature, title, "NumSamples",
                            [feature,'Count'])
            st.caption(
                "The x-axis shows the levels of the grouping chosen (clusters or baskets/tissues). The y-axis shows the number of samples.")

    #Displays results in a table format, option to show results grouped by response (RD) within groups
    def showRawData(self, feature, x_variable, RD):
        if RD:
            features = ["responsive", feature]
            columns = ["responsive",x_variable, "Number of samples"]
        else:
            features = feature
            columns = [x_variable, "Number of samples"]
        raw_count = self.patient_df.groupby(features).size().to_frame(name = 'count').reset_index()
        raw_count.columns = columns
        if RD:
            raw_count['responsive'] = raw_count['responsive']== "Responsive"
        saveTable(raw_count, "NumOfS")
        st.dataframe(raw_count, use_container_width=True)

    #Shows samples' AAC response in a scatterplot or table if RawD is True
    def AAC_scatterplot(self, RawD):
        if RawD:
            saveTable(self.patient_df, "rawAAC")
            self.patient_df["responsive"] = self.patient_df["responsive"] == "Responsive"
            st.dataframe(self.patient_df, use_container_width=True)
        else:
            reseted_df = self.patient_df.reset_index()
            reseted_df["index"] = range(reseted_df.shape[0])
            alt_scatterplot(reseted_df, 'index', 'responses', "Sample index", "AAC response", "AAC response", "AAC", ["samples", "responses"])
            st.caption("The x-axis shows the sample index. The y-axis shows the real AAC response of the sample to the drug.")

    #Display AAC response statistics in a table format
    def raw_data_AAC(self,feature,x_variable):
        raw_count = self.patient_df[[feature, "responses"]].groupby(feature).mean()
        default_value = 0
        raw_count['SD'] = self.patient_df[[feature, "responses"]].groupby(feature).std().fillna(default_value)
        raw_count['Median'] = self.patient_df[[feature, "responses"]].groupby(feature).median()
        raw_count['Min'] = self.patient_df[[feature, "responses"]].groupby(feature).min()
        raw_count['Max'] = self.patient_df[[feature, "responses"]].groupby(feature).max()
        raw_count = raw_count.reset_index()
        raw_count.columns = [x_variable, "Mean", "SD", "Median", "Min", "Max"]
        return raw_count

    #Shows samples' AAC response in a boxplot (if RawD is False) or in a table (if RawD is True). Boxplot can also be displayed grouped by response within groups
    def AAC_response(self,feature, RD, x_lab, RawD):
        size = len(np.unique(self.patient_df[feature]))
        if RawD:
            df = Data.raw_data_AAC(self,feature, x_lab)
            saveTable(df, "rawAAC")
            st.dataframe(df, use_container_width=True)
        else:
            if RD:
                alt_boxplot(self.patient_df, feature, "responses", 2, x_lab, "AAC response", "responsive", "AAC response",
                            "AAC")
                st.caption(
                    "The x-axis shows the levels of the grouping chosen (clusters or baskets/tissues). The y-axis shows the real AAC response to the drug."
                    " Within levels, the AAC response by the responsive (blue) or non-responsive (red) samples to treatment are shown.")
            else:
                alt_boxplot(self.patient_df, feature, "responses", size, x_lab, "AAC response", feature, "AAC response", "AAC")
                st.caption("The x-axis shows the levels of the grouping chosen (clusters or baskets/tissues). The y-axis shows the real AAC response to the drug.")

    #Displays mean inferred response probabilities for clusters o baskets in a bar plot or table (if RawD_prob is True)
    def barInferredProb(self, feature,RawD_prob):
        stacked = self.stacked_posterior
        if feature == 'baskets':
            len_colors = len(self.baskets_names)
            inferred_basket = stacked.basket_p.values
            baskets = []
            for item in self.baskets_names:
                for _ in range(len(inferred_basket[0])):
                    baskets.append(item)
            inferred_basket = inferred_basket.flatten()
            df = pd.DataFrame({'mean probability': inferred_basket, 'tissue': baskets})
            title = 'Tissue/Basket'
            feature = 'tissue'
            subheader = "Inferred basket response"
        elif feature == "clusters":
            len_colors = len(self.clusters_names)
            inferred_cluster = stacked.cluster_p.values
            clusters = []
            for item in self.clusters_names:
                for _ in range(len(inferred_cluster[0])):
                    clusters.append(item)
            inferred_cluster = inferred_cluster.flatten()
            df = pd.DataFrame({'mean probability': inferred_cluster, 'cluster': clusters})
            title = 'Clusters'
            feature = 'cluster'
            subheader = "Inferred cluster response"
        if RawD_prob:
            saveTable(df, "raw-prob")
            st.dataframe(df, use_container_width=True)
        else:

            alt_boxplot(df, feature, 'mean probability',len_colors, title, "Inferred response probability", feature, subheader, "Inferred")
            st.caption("The x-axis shows the levels of the grouping chosen (clusters or baskets/tissues). "
                       "The y-axis shows the inferred response probability to the treatment."
            )

    #Displays ECDF results for a chosen cluster or basket in a plot or table (if RawD is true)
    def ecdf_indiv(self, feature, choice, index, RawD,cred_inter):
        intervals = np.array([100-cred_inter-5, 50, cred_inter+5])
        if feature == 'baskets':
            basket_data = self.stacked_posterior.basket_p[index]
            pct, pct_val = ecdf(basket_data,intervals)
            title = 'Basket ' + choice
        elif feature == "clusters":
            cluster_data = self.stacked_posterior.cluster_p[index]
            pct, pct_val = ecdf(cluster_data,intervals)
            title = 'Cluster '+ str(choice)
        st.write("""##### {}th, 50th and {}th Percentile values are""".format(intervals[0],intervals[2]) + ": {0:.2f}, {1:.2f} and {2:.2f}".format(*pct_val['x']))
        if RawD:
            saveTable(pct, "raw-ecdf")
            st.dataframe(pct, use_container_width=True)
        else:
            alt_line_chart(pct,pct_val, 'Probability', 'Percent', 'Probability', 'Percent', "Cumulative distribution for "+title,"ecdf")
            st.caption(
                "The x-axis represents the probabilities points. The y-axis represents the proportion or fraction of data points that are less than or equal to a given value.")
