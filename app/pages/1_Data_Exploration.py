import streamlit as st
from analysis import Analysis,heatMap
from common import add_logo, hideRows, savePlot,sideBar, openGeneCard,searchTranscripts
from interpretation import Prototypes, DEA
from streamlit_option_menu import option_menu

hide_rows = hideRows()
add_logo()
sideBar()
st.header("Data exploration")
st.write("---")
menu = option_menu(None, ["Samples information", "pyBasket results","Statistics"],
    icons=["bi bi-clipboard-data", "bi bi-basket","bi bi-graph-up"],
    menu_icon="cast", default_index=0, orientation="horizontal")

if "data" in st.session_state:
    data = st.session_state["data"]
    #Create a heatmap object
    heatmap = heatMap(st.session_state["saved data"],st.session_state["File Name"])
    analysis = st.session_state["Analysis"]
    #Subpage showing information about the samples: Number and AAC response
    if menu == "Samples information":
        st.subheader("Number of samples")
        st.write(
            "The number of samples in each cluster or disease type is shown below. Using 'Group by resistance' option,"
            "samples will be coloured based on their resistance against the drug: Resistant (red) and Non-resistant (green).")
        st.write(" ")
        col11, col12 = st.columns((1, 1))
        with col11:
            #Option to select level of grouping: clusters or baskets
            option_tab1 = st.selectbox("Select how to group samples", ('Clusters', 'Disease types'),
                                       key="option-tab1")
        with col12:
            st.write(" ")
            #Option to group samples by response within groups
            RD = st.checkbox("Group by resistance", key="responsive")
            #Option to show data in a table format
            RawD = st.checkbox("Show raw data", key="raw-data")
        if option_tab1 == "Clusters":
            if RawD:
                data.showRawData("cluster_number","Cluster number", RD)
            else:
                data.countPlot(RD,"Number of samples","cluster_number","Cluster number")
        elif option_tab1 == 'Disease types':
            if RawD:
                data.showRawData("disease", "Disease type", RD)
            else:
                data.countPlot(RD, "Number of samples per disease type", "disease", "Disease types")
        st.write("---")
        st.subheader("AAC response")
        st.write(
            "The Area Above the Curve (AAC) is a measure used to analyse drug response and quantify the effect of"
            " a drug over a period of time. In the context of the Genomics of Drug Senstivity in Cancer (GDSC) project, "
            "the AAC drug response is the measure of the cell"
            " or cell line overall survival in response to the drug: the larger the AAC metric, the more resistance to the"
            " drug is shown. A larger AAC response implies a larger resistance to the drug."
            " "
            "The AAC response or resistance values per cluster o disease type are shown below, as well as within groups for Resistant (red) and Non-"
            "resistant (green) samples.")
        st.write(" ")
        col21, col22 = st.columns((1, 1))
        with col21:
            # Option to select level of colouring: none, clusters or baskets
            option_tab2 = st.selectbox("Select how to group samples", ('None', 'Clusters', 'Disease types'),
                                       key="option-tab2")
        with col22:
            st.write(" ")
            st.write(" ")
            #Option to show data in a table format
            RawD_AAC = st.checkbox("Show raw data", key="raw-data-AAC")
        if option_tab2 == 'None':
            data.AAC_scatterplot(RawD_AAC)
        else:
            with col22:
                # Option to group samples by response within groups
                RD_AAC = st.checkbox("Group by resistance", key="responsive-AAC")
            if option_tab2 == "Clusters":
                data.AAC_response("cluster_number", RD_AAC, "Cluster number", RawD_AAC)
            elif option_tab2 == 'Disease types':
                data.AAC_response("disease", RD_AAC, "Disease type", RawD_AAC)
    # Subpage showing information about pyBasket pipeline results
    elif menu == "pyBasket results":
        st.subheader("Inferred response probability distribution")
        st.write("In the second stage of the pyBasket pipeline, a hierarchical Bayesian model is used to estimate the"
                 " overall response probability distribution of each cluster or disease type."
                 " By choosing below a credible interval, the range containing a particular percentage of probable values is shown in the plot. "
                 "This range is delimited by the percentiles below (lower bound) and above (upper bound) the credible interval and the median (50th). Percentiles represent the value below "
                 "which a given percentage of probabilities fall, e.g. a 5th percentile is the value below which 5% of the data points fall.")
        col1, col2 = st.columns((2,2))
        with col1:
            # Option to select level of grouping: clusters or baskets
            option_page2 = st.selectbox("Select group", ('Clusters', 'Disease types'),
                                        key="option-page2")
        with col2:
            cred_inter = st.number_input('Credible interval', value=90)
            st.caption("90% credible interval shown by default")
        # Option to show data in a table format
        RawD_prob = st.checkbox("Show raw data", key="raw-data-prob")
        if option_page2 == 'Disease types' and cred_inter <101 and cred_inter>0:
            data.barInferredProb("baskets",RawD_prob,cred_inter)
        elif option_page2 == "Clusters" and cred_inter <101 and cred_inter>0:
            data.barInferredProb("clusters",RawD_prob,cred_inter)
        else:
            st.warning("Please select a credible interval between 0 and 100.")
    # Subpage for analysis: PCA, prototypes and DEA
    elif menu == "Statistics":
        tab21, tab22, tab23 = st.tabs(["Dimensionality reduction", "Prototypes", "Differential expression"])
        #Dimensionality reduction subpage (PCA)
        with tab21:
            st.subheader("Dimensionality reduction")
            st.write("The goal of dimensionality reduction is to project high-dimensionality data to a lower "
                     "dimensional subspace while preserving the essence of the data and the maximum amount of information.")
            st.write("Principal Component Analysis (PCA) is a dimensionality reduction method that enables the visualisation"
                     " of high-dimensional data as well as capturing as much as possible of its variation. "
                     "The results for PCA using 5 Principal Components on the data can be explored by selecting a variable to colour cell lines or samples by: "
                     "clusters, disease type or resistance.")
            st.write(" ")
            # Option to select level of grouping: clusters, baskets or response
            option_PCA = st.selectbox("Select how to group samples", ('Clusters', 'Disease types', 'Resistance'), key="PCA")
            #Option to show data in table format
            RawD = st.checkbox("Show raw data", key="raw-data-PCA")
            if option_PCA == "Clusters":
                analysis.PCA_analysis("cluster_number", RawD)
            elif option_PCA == "Resistance":
                analysis.PCA_analysis("resistance", RawD)
            elif option_PCA == 'Disease types':
                analysis.PCA_analysis("disease", RawD)
        #Prototypes samples subpage
        with tab22:
            st.subheader("Prototypes")
            st.write("The prototypical sample of each cluster, disease type or pattern of response has been calculated using"
                     " K-Medoids. K-Medoids is a partitioning technique that finds the sample (medoid or prototype)"
                    " that is the closest to the rest of samples in the same group. The dimension of the expression level of transcripts"
                     " has been reduced and plotted using Principal Component Analysis (PCA) as described in the 'Dimensionality reduction' tab.")
            st.write(" ")
            #Option to select level of grouping: clusters, baskets or response
            option_prototype = st.selectbox("Select how to group samples", ('Clusters', 'Disease types'), key="Prototypes")
            prototype = Prototypes(data)
            prototype.findPrototypes(option_prototype)
        #Differential expression analysis subpage
        with tab23:
            st.subheader("Differential expression")
            st.write(" ")
            st.write("The goal of Differential Expression Analysis (DEA) is to discover whether the expression "
                     "level of a feature (gene or transcript) is quantitatively different between experimental "
                     "groups or conditions. Below, any two clusters, types of disease or resistance groups can be compared. ")
            st.write("The T-test for the means of each feature in two independent groups or conditions is calculated."
                     "The null hypothesis is that the feature has identical average values across conditions.")
            st.write(" ")
            col51, col52 = st.columns((2, 2))
            with col51:
                # Option to select level of grouping: clusters, baskets or response
                option = st.selectbox("Select how to group samples", ('Clusters', 'Disease types', 'Resistance'), key="DEA")
                if option == 'Clusters':
                    subgroups = data.clusters_names
                    feature = "cluster_number"
                elif option == 'Disease types':
                    subgroups = data.disease_types
                    feature = "disease"
                else:
                    subgroups = ["Resistant","Non-resistant"]
                    feature = "resistance"
                #Option to select groups to perform DEA
                groups = st.multiselect(
                    'Please select up to 2 Groups to compare', subgroups, max_selections=2)
            if len(groups)<2:
                st.write("")
            else:
                with col51:
                    #User input for corrected p-value threshold
                    pthresh = st.number_input('P-value threshold for significance (0.05 by default)', value = 0.05)
                    st.caption("Applied on corrected p-values.")
                    #User input for log2 fold change threshold
                    logthresh = st.number_input('log2 FC threshold for significance (1 by default)', value=1.0)
                dea = DEA(data)
                st.write("#### DEA on {} compared: {} vs {}".format(option,groups[0], groups[1]))
                dea.diffAnalysis_simple(groups[0],groups[1],feature,pthresh,logthresh)
                results = dea.showResults(feature)
                with col52:
                    st.write(" ")
                    st.write(" ")
                    dea.infoTest(groups[0],groups[1],option,pthresh,logthresh)
                st.subheader("Single Transcript DEA")
                st.write("The difference in the expression level of a selected transcript between groups being compared can be explored below. "
                         "Further information about the transcript can be found by using the GeneCards button. ")
                st.write("")
                col53, col54 = st.columns((2,4))
                with col53:
                    # Option to select transcript results to be displayed
                    transcript = searchTranscripts(results["Feature"])
                    st.write(" ")
                    dea.infoTranscript(transcript)
                with col54:
                    dea.boxplot(groups[0], groups[1], feature, transcript)
















