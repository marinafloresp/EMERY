import os
import numpy as np
import streamlit as st
import pandas as pd
import altair as alt
from common import saveTable,alt_ver_barplot,alt_scatterplot,alt_boxplot, ecdf,alt_line_chart, colours, savePlot

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
        self.disease_types = sorted(list(self.patient_df['tissues'].unique()))

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

    #Changes resistant feature to be categorical (Resistant vs Non-resistant) and maps ENSEMBL IDs for Genes' names
    def setData(self):
        self.patient_df["resistance"] = self.patient_df["responsive"].replace([0, 1], ['Resistant', 'Non-resistant'])
        self.patient_df["disease"] = self.patient_df["tissues"]
        self.patient_df = self.patient_df.drop(['responsive', 'tissues'], axis=1)
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
        tissues = len(self.disease_types)
        drug = self.file_name.split('_')[2]
        info = {'File name: ': {'information': self.file_name}, 'Drug: ': {'information': drug},'Num. cell lines' : {'information':num_samples},
                'Num. transcripts (after filtering): ': {'information': num_feat_filt},
                'Num. transcripts (after feature selection): ': {'information': num_feat_select},
                'Num.clusters: ': {'information': clusters}, 'Num. disease types: ': {'information': tissues}}
        df = pd.DataFrame(data=info).T
        style = df.style.hide_index()
        style.hide_columns()
        st.dataframe(df, use_container_width=True)

    #Shows information about the number of samples per group in a bar plot. If RD (Response Data) is True, samples are grouped by response within groups.
    def countPlot(self, RD, title,feature,x_lab):
        self.patient_df[feature] = pd.Categorical(self.patient_df[feature])
        size = len(np.unique(self.patient_df[feature]))
        if RD:
            df_grouped = self.patient_df.groupby(["resistance", feature]).size().reset_index(name='Count')
            alt_ver_barplot(df_grouped, feature, 'Count', 2, x_lab, "Number of samples", "resistant", title, "NumSamples",
                            ["resistance", feature, 'Count'])
            st.caption("The x-axis shows the levels of the grouping chosen (clusters or disease types). The y-axis shows the number of samples."
                       " Within levels, the number of samples that are Non-resistant (green) or Resistant (red) to treatment are shown.")
        else:
            df_grouped = self.patient_df.groupby([feature]).size().reset_index(name='Count')
            alt_ver_barplot(df_grouped, feature, 'Count', size, x_lab, "Number of samples", feature, title, "NumSamples",
                            [feature,'Count'])
            st.caption(
                "The x-axis shows the levels of the grouping chosen (clusters or disease types). The y-axis shows the number of samples.")

    #Displays results in a table format, option to show results grouped by response (RD) within groups
    def showRawData(self, feature, x_variable, RD):
        if RD:
            features = ["resistance", feature]
            columns = ["resistance",x_variable, "Number of samples"]
        else:
            features = feature
            columns = [x_variable, "Number of samples"]
        raw_count = self.patient_df.groupby(features).size().to_frame(name = 'count').reset_index()
        raw_count.columns = columns
        if RD:
            raw_count['resistance'] = raw_count['resistance']== "Non-resistant"
        saveTable(raw_count, "NumOfS")
        st.dataframe(raw_count, use_container_width=True)

    #Shows samples' AAC response in a scatterplot or table if RawD is True
    def AAC_scatterplot(self, RawD):
        if RawD:
            saveTable(self.patient_df, "rawAAC")
            self.patient_df["resistance"] = self.patient_df["resistance"] == "Non-resistant"
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
                alt_boxplot(self.patient_df, feature, "responses", 2, x_lab, "AAC response", "resistance", "AAC response",
                            "AAC")
                st.caption(
                    "The x-axis shows the levels of the grouping chosen (clusters or disease types). The y-axis shows the real AAC response to the drug."
                    " Within levels, the AAC response by the non-resistant (blue) or resistant (red) samples to treatment are shown.")
            else:
                alt_boxplot(self.patient_df, feature, "responses", size, x_lab, "AAC response", feature, "AAC response", "AAC")
                st.caption("The x-axis shows the levels of the grouping chosen (clusters or disease types). The y-axis shows the real AAC response to the drug.")

    #Displays mean inferred response probabilities for clusters o baskets in a bar plot or table (if RawD_prob is True)
    def barInferredProb(self, feature,RawD_prob,cred_inter):
        stacked = self.stacked_posterior
        if feature == 'baskets':
            len_colors = len(self.disease_types)
            inferred_basket = stacked.basket_p.values
            baskets = []
            for item in self.disease_types:
                for _ in range(len(inferred_basket[0])):
                    baskets.append(item)
            inferred_basket = inferred_basket.flatten()
            df = pd.DataFrame({'mean probability': inferred_basket, 'disease': baskets})
            feature = 'disease'
        elif feature == "clusters":
            len_colors = len(self.clusters_names)
            inferred_cluster = stacked.cluster_p.values
            clusters = []
            for item in self.clusters_names:
                for _ in range(len(inferred_cluster[0])):
                    clusters.append(item)
            inferred_cluster = inferred_cluster.flatten()
            df = pd.DataFrame({'mean probability': inferred_cluster, 'cluster': clusters})
            feature = 'cluster'
        palette = colours(len_colors)
        intervals = [(50 - (cred_inter/2)) / 100, 0.5, (50+(cred_inter/2)) / 100]
        interval_data = (
            df.groupby(feature)['mean probability']
            .quantile(intervals)
            .reset_index()
            .rename(columns={'level_1': 'Percentile th', 'mean probability': 'Range'})
        )
        interval_data['Percentile th'] = interval_data['Percentile th']*100
        intervals =[int(x*100) for x in intervals]
        st.write("##### Percentiles of chosen credible interval are: {}th (lower bound/min of range), {}th (median) and {}th (upper bound/max of range".format(*intervals))
        if RawD_prob:
            saveTable(interval_data, "raw-prob")
            st.dataframe(interval_data, use_container_width=True)
        else:
            base = alt.Chart(interval_data, title="Inferred response probabilities").mark_boxplot(ticks=True, size=50).encode(
                x=feature +':N',
                y='Range:Q',
                color=feature +':N'
            ).properties(height=650, width=600
            ).configure_range(category=alt.RangeScheme(palette))
            savePlot(base, "probabilities")
            st.altair_chart(base, theme="streamlit", use_container_width=True)
            st.caption("The x-axis shows the different groups in the variable chosen (clusters or disease types)."
                       " The y-axis shows the range defined by the credible interval in which the overall probability of response falls. ")
