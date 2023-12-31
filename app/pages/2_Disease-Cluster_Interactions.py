import streamlit as st
from analysis import Analysis, heatMap
from interpretation import Prototypes, DEA
from common import add_logo,savePlot,sideBar, saveTable, openGeneCard,savePlot_plt,searchTranscripts,ecdf_interaction
from streamlit_option_menu import option_menu

add_logo()
sideBar()
st.header("Disease-Cluster interactions")
st.write("---")
menu = option_menu(None, ["All interactions", "Selected interaction"],
    icons=["bi bi-asterisk", "bi bi-braces-asterisk"],
    menu_icon="cast", default_index=0, orientation="horizontal")

if "data" in st.session_state:
    data = st.session_state["data"]
    analysis_data = st.session_state["Analysis"]
    #Create a Heatmap object
    heatmap = heatMap(st.session_state["saved data"],st.session_state["File Name"])
    #Find samples in disease-cluster interaction
    disease, cluster = st.session_state["disease"], st.session_state["cluster"]
    subgroup, size = analysis_data.findInteraction(st.session_state["cluster"], st.session_state["disease"])
    #Subpage with overview of all interactions
    if menu == "All interactions":
        st.write("""This page shows results and information about the interactions between clusters and types of diseases. Below, the information to 
                 display in the heatmap can be selected, such as the number of samples, the number of Non-resistant samples or the inferred response in each disease-cluster interaction.""")
        #Option to select the information to display in the heatmap
        heatmap_info = st.radio("Select information to display",
                            ['Number of samples', 'Number of Non-resistant samples', "Inferred response"],
                            key="HM_info", horizontal=True)
        st.write("")
        #Option to mark interactions that have a minimum number of samples
        min_num = st.slider('Mark interactions with a minimum number of samples', 0, 70)
        st.info("\:star: : disease-cluster interactions with at least {} samples.\n"
                "\:large_red_square: : selected disease+cluster interaction.\n".format(min_num))
        st.write("")
        num_samples = heatmap.heatmapNum()
        if heatmap_info == 'Number of samples':
            st.write("#### Number of samples per disease-cluster interaction")
            st.write("Explore the number of samples in each type of disease and cluster combination.")
            #Option to show data in a table
            RawD = st.checkbox("Show raw data", key="raw-data-S")
            if RawD:
                saveTable(num_samples, "NumSamplesHM")
                st.dataframe(num_samples, use_container_width=True)
            else:
                HM_NoS = heatmap.heatmap_interaction(num_samples,num_samples,"Number of samples per interaction"
                                                         ,min_num, int(st.session_state["cluster"]), st.session_state["disease"])
                savePlot_plt(HM_NoS, str(cluster) + "_" + disease)
                st.pyplot(HM_NoS)
                st.caption("The x-axis shows the levels of clusters. The y-axis shows the types of disease. The colour scale represents "
                           "the number of samples.")
        elif heatmap_info == 'Number of Non-resistant samples':
            st.write("#### Number of samples per disease-cluster interaction showing non-resistance to the drug")
            st.write(
                "Explore the number of samples in each disease and cluster combination that are Non-resistant (or responsive) to the drug.")
            #Option to show data in a table
            RawD = st.checkbox("Show raw data", key="raw-data-R")
            response_df = heatmap.heatmapResponse()
            if RawD:
                saveTable(response_df, "responseHM")
                st.dataframe(response_df, use_container_width=True)
            else:
                HM_response = heatmap.heatmap_interaction(response_df,num_samples,"Non-resistant samples per interaction", min_num,
                                                          int(cluster), disease)
                savePlot_plt(HM_response, str(cluster) + "_" + disease)
                st.pyplot(HM_response)
                st.caption(
                    "The x-axis shows the levels of clusters. The y-axis shows the types of diseases. The colour scale represents "
                    "the number of samples that show non-resistance (responsive) to the drug.")
        elif heatmap_info == "Inferred response":
            st.write("#### Inferred response probability per disease-cluster interaction.")
            st.write(
                "Explore the inferred response probability of each disease-cluster interaction calculated by the HBM model in the second stage of "
                "the pyBasket pipeline based on observed responses.")
            # Option to show data in a table
            RawD = st.checkbox("Show raw data", key="raw-data-R")
            inferred_df = heatmap.HM_inferredProb()
            if RawD:
                saveTable(inferred_df,"inferredHM")
                st.dataframe(inferred_df, use_container_width=True)
            else:
                HM_inferred = heatmap.heatmap_interaction(inferred_df,num_samples,"Inferred disease*cluster interaction", min_num,
                                                          int(cluster), disease)
                savePlot_plt(HM_inferred, "inferred_heatmap")
                st.pyplot(HM_inferred)
                st.caption(
                    "The x-axis shows the levels of clusters. The y-axis shows the types of diseases. The colour scale represents "
                    "the inferred probability of response to the treatment by the disease-cluster interaction.")
    # Subpage to show information and perform analysis on the samples included in the selected disease-cluster interaction
    elif menu == "Selected interaction":
        if len(subgroup) > 1:
            st.text("")
            st.write("In this page, results specific to the disease-cluster interaction selected in the left sidebar can be explored.")
            st.text("")
            st.info("###### Samples in **cluster {}** & **{} disease**: {}".format(cluster, disease, size))
            st.text("")
            tab1, tab2, tab3, tab4 = st.tabs(["Overview", "PCA", "Prototypes", "Differential Expression"])
            #Subpage with overview of the selected disease-cluster interaction if enough samples are found
            with tab1:
                st.subheader("Response to drug")
                st.markdown("Number of samples in **cluster {}** & **disease {}**: {} ".format(str(cluster), disease, size))
                col21, col22 = st.columns((2, 2))
                with col21:
                    analysis_data.samplesCount(subgroup)
                with col22:
                    analysis_data.responseSamples(subgroup)
                    st.caption("Samples ordered from lestas to most resistance (lower to higher AAC response). Responsive (non-resistant) samples are marked with a ticked box.")
                st.subheader("Inferred response probabilities")
                st.write(
                    "In the second stage of the pyBasket pipeline, a hierarchical Bayesian model is used to estimate the overall probability "
                    "of a disease-cluster interaction to be responsive to the treatment based on the observed AAC drug response metrics. In order to "
                    "account for the number of samples in an interaction, the Empirical Cumulative Distribution Function (ECDF) has "
                    "been applied to the probabilities output by the pyBasket pipeline. The ECDF represents the cumulative probability "
                    "distribution of a dataset based on observed data points and is used to estimate or approximate the underlying Probability "
                    "Density Function (PDF)."
                    " The PDF is applied to these values in order to visualise the relative likelihood of different values occurring within the distribution"
                    " of response probabilities and inform about the distribution's shape, central tendencies, and dispersion.")
                st.write(" The credible interval can be chosen below, which is the range containing a particular percentage of probable values (shadowed in gray), e.g. for the 90% "
                    "credible interval, the range containing the 90% of the values is shown. The bottom of the range or lower bound (red vertical line), the median (black vertical "
                    "line) and the top of the range or upper bound (blue vertical line) are shown.")
                #User input option to specify credible interval for ECDF
                cred_inter = st.number_input('Credible interval', value=90)
                st.caption("90% credible interval shown by default")
                # Option to show data in a table
                RawD_ecdf = st.checkbox("Show raw data", key="raw-data-ecdf")
                if cred_inter <101 and cred_inter>0:
                    ecdf_interaction(data,disease, cluster, RawD_ecdf, cred_inter)
                else:
                    st.warning("Please select a credible interval between 0 and 100.")
                st.caption("The x-axis shows the distribution of the joint probability values. The y-axis shows the density of such distribution."
                           " The red-line curve is the PDF that describes the distribution of probabilities across the possible values of the variable."
                           " The black vertical line is the median. The red vertical line is the lower bound of the credible interval. The blue vertical"
                           " line is the upper bound of the credible interval.")
            # Subpage with PCA on samples in the selected disease-cluster interaction if enough samples are found
            with tab2:
                st.subheader("Advanced PCA")
                st.write("")
                st.write(
                    "Principal Component Analysis (PCA) is a dimensionality reduction method that enables the visualisation"
                    " of high-dimensional data. The results for PCA on the data can be explored for the samples in the selected"
                    " disease-cluster interaction, that are grouped by resistance: Resistant (red) vs Non-resistant (green).")
                st.write(" ")
                # st.write("##### PCA of samples in **cluster {}** & **disease {}**".format(cluster, disease))
                RawD = st.checkbox("Show raw data", key="raw-data")
                analysis_data.adv_PCA(subgroup, RawD)
            # Subpage with Prototypical samples in the selected disease-cluster interaction if enough samples are found
            with tab3:
                st.subheader("Prototypes of subgroup")
                st.write("")
                st.write("The prototypical samples of the selected disease-cluster interaction for the Resistant and the"
                         " Non-resistant groups have been calculated using KMedoids. KMedoids finds the sample that is"
                         "the closest to the rest of samples in the group and therefore representative of the group.")
                sub_prototypes = Prototypes(data)
                sub_prototypes.findPrototypes_sub(subgroup)
            # Subpage with DEA with samples in the selected disease-cluster interaction if enough samples are found
            with tab4:
                st.subheader("Differential expression analysis (DEA)")
                st.write("")
                st.write("The goal of Differential Expression Analysis (DEA) is to discover whether the expression "
                         "level of a feature (gene or transcript) is quantitatively different between experimental "
                         "groups or conditions."
                         " T-test for the means of each feature in two independent groups or conditions is calculated."
                         "The null hypothesis is that the feature has identical average values across conditions.")
                st.write(" ")
                try:
                    dea = DEA(data)
                except:
                    st.warning(
                        "Not enough samples in the disease-cluster interaction. Please select another combination.")
                col41, col42 = st.columns((2, 2))
                with col41:
                    # Option to select the type of DEA to perform
                    option = st.selectbox("Select analysis", (
                        "Samples in interaction vs rest of samples", "Within interaction: resistant vs non-resistant",
                        "Interaction vs Interaction"), key="DEA")
                    # User input to apply correct p-value threshold for significance
                    pthresh = st.number_input('P-value threshold for significance (0.05 by default)', value=0.05)
                    st.caption("Applied on corrected p-values.")
                    # User input to apply log2FC threshold for significance
                    logthresh = st.number_input('log2 FC threshold for significance (1 by default)', value=1.0)
                if option == "Samples in interaction vs rest of samples":
                    st.write(" ")
                    st.subheader("Samples in interaction vs rest of samples")
                    st.write(
                        "DEA performed with the samples in the selected disease-cluster interaction vs rest of samples.")
                    if size > 2:
                        dea.diffAnalysis_inter(subgroup, pthresh, logthresh)
                        results = dea.showResults("interaction")
                    else:
                        st.warning("Not enough samples (minimum 3 samples). Please try a different combination.")
                    with col42:
                        st.write(" ")
                        st.write(" ")
                        try:
                            dea.infoTest((cluster, disease), 'All', 'Interaction', pthresh, logthresh)
                        except:
                            st.write(" ")
                    try:
                        st.subheader("Individual transcripts DEA")
                        st.write(
                            "The difference in the expression level of a transcript/feature and the individual results from "
                            "DEA can be explored below. Further information about the transcript can be found in the button that links to its GeneCard profile.")
                        transcript = searchTranscripts(results["Feature"])
                        col53, col54 = st.columns((2, 3))
                        with col53:
                            st.write(" ")
                            dea.infoTranscript(transcript)
                        with col54:
                            dea.boxplot_inter(subgroup, transcript)
                    except:
                        st.warning("Not enough samples. Please try a different combination.")
                elif option == "Within interaction: resistant vs non-resistant":
                    st.write(" ")
                    st.subheader("Resistant vs non-resistant samples within disease-cluster interaction")
                    st.write(
                        "DEA has been performed within the selected interaction by comparing Resistant vs Non-resistant samples.")
                    if size > 3:
                        dea.diffAnalysis_response(subgroup, pthresh, logthresh)
                        try:
                            results = dea.showResults("interaction")
                            st.subheader("Individual transcripts DEA")
                            st.write(
                                "The difference in the expression level of a transcript/feature and the individual results from "
                                "DEA can be explored below. Further information about the transcript can be found in the button that links to its GeneCard profile.")
                            col53, col54 = st.columns((2, 4))
                            with col53:
                                feature = searchTranscripts(results["Feature"])
                                st.write(" ")
                                dea.infoTranscript(feature)
                            with col54:
                                dea.boxplot_resp(subgroup, feature)
                        except:
                            st.warning("Not enough samples. Please try a different combination.")
                    else:
                        st.warning("Not enough samples (minimum 4 samples). Please try a different combination.")
                    with col42:
                        st.write(" ")
                        st.write(" ")
                        try:
                            dea.infoTest('Resistant', 'Non-resistant', (cluster, disease), pthresh, logthresh)
                        except:
                            st.write(" ")
                elif option == "Interaction vs Interaction":
                    st.write(" ")
                    if subgroup.size > 3:
                        with col41:
                            st.write("##### Select a second interaction")
                            col43, col44 = st.columns((2, 2))
                            with col43:
                                cluster_inter2 = st.selectbox("Select a cluster", data.clusters_names, key="cluster2")
                            with col44:
                                disease_inter2 = st.selectbox("Select a disease", data.disease_types, key="disease2")
                            subgroup2, size2 = analysis_data.findInteraction(cluster_inter2,
                                                                             disease_inter2)
                            st.info("###### Samples in **cluster {}** & **{} disease**: {}".format(cluster_inter2,
                                                                                                   disease_inter2,
                                                                                                   size2))
                        if size2 > 3 and (cluster_inter2 != cluster or disease_inter2 != disease):
                            st.subheader(
                                "Samples in ({}-cluster {}) vs samples in ({}-cluster {})".format(disease, int(cluster),
                                                                                                  disease_inter2,
                                                                                                  int(cluster_inter2)))
                            st.write(
                                "DEA has been performed with samples in the selected interaction compared to samples in a second interaction.")
                            dea.diffAnalysis_adv(subgroup, subgroup2, pthresh, logthresh)
                            results = dea.showResults("interaction")
                            st.subheader("Individual transcripts DEA")
                            st.write(
                                "The difference in the expression level of a transcript/feature and the individual results from "
                                "DEA can be explored below. Further information about the transcript can be found in the button that links to its GeneCard profile.")
                            col53, col54 = st.columns((2, 4))
                            with col53:
                                feature = searchTranscripts(results["Feature"])
                                st.write(" ")
                                dea.infoTranscript(feature)
                            with col54:
                                dea.boxplot_adv(subgroup, subgroup2, disease, disease_inter2, cluster, cluster_inter2,
                                                feature)
                            with col42:
                                st.write(" ")
                                st.write(" ")
                                st.write(" ")
                                dea.infoTest((disease, int(cluster)), (disease_inter2, int(cluster_inter2)),
                                             "Comparison", pthresh, logthresh)
                        else:
                            st.warning(
                                "Not enough samples in second interaction selected (minimum 4 samples). Please try a different combination.")
                    else:
                        st.warning(
                            "Not enough samples in first interaction selected (minimum 4 samples). Please try a different combination.")
        else:
            st.warning("Not enough samples. Please try a different combination.")
