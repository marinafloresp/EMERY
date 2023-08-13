import streamlit as st
from importance import FI, Global_ALE, Shap_FI
from common import add_logo, saveTable, savePlot,sideBar,openGeneCard,searchTranscripts,savePlot_plt
from streamlit_option_menu import option_menu

add_logo()
sideBar()
st.header("Transcripts (Feature) importance")
st.write("---")

menu = option_menu(None, ["Overview", "Important Transcripts (Global)", "Important Transcripts (Local)"],
    icons=["bi bi-bar-chart", "bi bi-globe", "bi bi-pin-map"],
    menu_icon="cast", default_index=0, orientation="horizontal")

if "data" in st.session_state:
    data = st.session_state["data"]
    disease, cluster = st.session_state["disease"], st.session_state["cluster"]
    analysis_data = st.session_state["Analysis"]
    #Create main Feature Interaction object
    feat_inter = FI(data)
    #Create SHAP object
    shapFI = Shap_FI(data)
    explainer, values = shapFI.SHAP()
    if menu == "Overview":
        st.write(" ")
        st.subheader("Overview of Important Transcripts")
        st.write('The transcripts that are driving the prediction of the AAC response of samples overall can be explored using Feature Importance Analysis. '
                 'Feature Importance is a technique used to describe how relevant a feature is and its effect on the model '
                 'being used to predict an outcome and can be calculated using several Model-Agnostic (MA) methods.')
        col11, col12 = st.columns((3,2))
        with col11:
            #Option to select the method to use
            overview_model = st.selectbox("Select a method", ["RF feature importance", "SHAP", "Permutation Based"], key="overview-model")
        with col12:
            #Option to select the number of top features to display
            num_feat = st.number_input('Select top number of features to display (<30 is encouraged)', value=10)
        if overview_model == "RF feature importance":
            st.write("#### Random Forest")
            st.write("Ramdom Forest is a common ML model that combines the output of multiple decision trees to reach a single result. It has been used "
                     "as part of the pydisease pipeline to select the 500 most important features. Features importance has been ranked based on the decrease in the impurity measure. "
                     "Below is shown the top {} features that are most important, i.e. their inclusion in a Random Forest will lead to a significant decrease in impurity. "
                     "Features are ranked based on their importance, higher values indicate greater importance. ".format(num_feat))
            #Option to show plot or raw data in a table
            RawD = st.checkbox("Show raw data", key="rd-RF")
            feat_inter.plotImportance(RawD,num_feat)
            st.caption("The x-axis represents the feature's importance (decrease in impurity measure).  The y-axis represents the features with the highest importance. ")
        elif overview_model == "SHAP":
            st.write("#### SHAP values")
            st.write("SHAP (SHapley Additive exPlanations) is a technique that explains the prediction of an observation by "
                         "computing the contribution of each feature to the prediction. It is based on Shapley values from game theory: it uses fair "
                     "allocation results from cooperative game to quantify the contribution that each feature makes to the model's prediction."
                         "Below, the top {} features that most influence the prediction of the AAC response to the drug are shown".format(num_feat))
            #Option to show plot or raw data in a table
            RawD = st.checkbox("Show raw data", key="rd-SHAP")
            explainer, values = shapFI.SHAP()
            fig = shapFI.SHAP_summary(values, num_feat, RawD)
            st.caption("The x-axis represents the magnitude of the feature's impact on the model's output."
                       " The y-axis shows the features with the highest impact ordered in decreasing order."
                       "Each point is a Shapley value of an instance per feature/transcript and the colour represents the feature's value.")
        elif overview_model == "Permutation Based":
            st.write("#### Permutation based feature importance")
            st.write('Permutation based feature importance is a MA method that measures'
                     ' the importance of a feature by calculating the increase in the model’s prediction '
                     'error when re-shuffling each predictor or feature. '
                     'Below is shown how much impact the top {} features have in the model’s prediction for the AAC response.'.format(num_feat))
            #Option to show plot or raw data in a table
            RawD = st.checkbox("Show raw data", key="rd-PBI")
            feat_inter.permutationImportance(num_feat,RawD)
            st.caption("The x-axis represents the feature's importance for the model's prediction. The y-axis represents the features ordered in decreasing order"
                       " of importance.")
    elif menu == "Important Transcripts (Global)":
        st.write(" ")
        st.subheader("Important Transcripts (Global)")
        st.write(" ")
        st.write('Here, the transcripts that have the highest impact on the AAC response prediction for a group of cell lines on average can be individually explored.'
                 'Global Model-Agnostic (MA) methods are used to describe the average behaviour of a Machine Learning model.   ')
        #Option to select the global method to use
        global_method = st.selectbox("Select a global method", ['ALE', 'PDP (SHAP)'], key ='global-method')
        if global_method == 'ALE':
            st.write("#### Accumulated Local Effects")
            st.write(" ")
            st.write("Accumulated Local Effects (ALE) describe how a transcript (feature) influences the prediction of AAC drug response. This method has been implemented so that the impact of a transcript's values"
                     " can be explored for all samples or "
                     "for a group of samples with specified conditions, such as only samples included in a selected disease*cluster interaction, a specific cluster or disease type."
                     "In addition, the impact of a transcript in two different groups of samples, e.g. two different clusters, can also be compared. ")
            st.write(
                "The resulting ALE plot shows how the prediction of the AAC drug response change as the transcript's expression values (represented by dots) move across different intervals (represented by lines in the x-axis).")
            global_ALE = Global_ALE(data)
            #Option to select the group of samples to display the feature's impact
            global_samples = st.radio("", ['All samples','Use samples in the selected interaction', 'Select only samples in cluster',
                                     'Select only samples in disease type'],
                                key="global_samples", horizontal=True)
            #List of transcripts and option to select
            transcripts = global_ALE.transcripts
            feature = searchTranscripts(transcripts)
            if global_samples == 'All samples':
                # Option to split the ALE results on Resistant vs Non-resistant samples.
                response = st.checkbox("Split by Resistant vs Non-resistant samples")
                if response:
                    global_ALE.global_ALE_resp(feature)
                else:
                    global_ALE.global_ALE(feature)
            elif global_samples == 'Use samples in the selected interaction':
                # Option to split the ALE results on Resistant vs Non-resistant samples.
                try:
                    global_ALE.global_ALE_single(feature, cluster,disease, "interaction")
                except:
                    st.warning("Not enough samples in the selected disease*cluster interaction. Please try a different combination.")
            else:
                if global_samples == 'Select only samples in cluster':
                    #Option to select the subgroup or subgroups of samples
                    groups = st.multiselect(
                        'Please select a cluster or up to 2 clusters to compare.', data.clusters_names, max_selections=2)
                    option = "clusters"
                elif global_samples == 'Select only samples in disease type':
                    # Option to select the subgroup or subgroups of samples
                    groups = st.multiselect(
                        'Please select a disease type or up to 2 disease types to compare.',data.disease_types, max_selections=2)
                    option = "disease"
                try:
                    if len(groups)<2:
                        global_ALE.global_ALE_single(feature, groups[0],None,option)
                    else:
                        global_ALE.global_ALE_mult(feature, groups[0], groups[1], option)
                except:
                    st.warning(
                        "Please select at least a group.")
        elif global_method == 'PDP (SHAP)':
            st.write("#### Partial Dependence Plot (PDP) by SHAP")
            st.write("The partial dependence plot (PDP) calculated by SHAP shows "
                     "the marginal effect that one or two features, that might interact with each other,"
                     " have on the predictions made by the ML model (the predicted AAC response to the drug).  ")
            features = shapFI.SHAP_results(values,len(values))
            #Option to select transcript
            transcript_PDP = st.selectbox("Select feature/transcript", features['Transcript'], key="transcript_PDP")
            st.caption("values ordered by decreasing importance by SHAP")
            st.write("#### SHAP dependence plot")
            shapFI.SHAP_dependence(values, transcript_PDP)
            st.caption("The x-axis is the actual feature value. "
                       "The y-axis represents the SHAP value for the feature: how much knowing that feature value changes the "
                       "output of the model for that prediction."
                       "The color is a second feature that interacts with the chosen feature.")
    elif menu == "Important Transcripts (Local)":
        st.subheader("Important Transcripts (Local)")
        st.write(" ")
        st.write("The transcripts that have the highest impact on the prediction of the AAC response for a specific cell lines can be explored."
                 " Local Model-Agnostic (MA) interpretation methods aim to explain individual predictions made by a Machine Learning model.  ")
        st.write("#### Filter samples")
        col31, col32, col33 = st.columns((2,2,2))
        with col31:
            #Option to filter samples to show for selection
            group_samples = st.radio("", ['Use samples in interaction', 'Select only samples in a cluster',
                                      'Select only samples in a disease type'],
                                 key="Local-samples")
        if group_samples == 'Select only samples in a cluster':
            with col33:
                #Option to select a cluster
                cluster= st.selectbox("Select a cluster", data.clusters_names, key="only_cluster")
            disease = "None"
        elif group_samples == 'Select only samples in a disease type':
            with col33:
                #Option to select a disease
                disease = st.selectbox("Select a disease type", data.disease_types, key="only_disease")
            cluster = "None"
        elif group_samples == 'Use samples in selection':
            disease, cluster = disease, cluster
        transcripts, size = feat_inter.displaySamples(cluster,disease)
        st.info("###### Samples in **cluster {}** & **{} disease**: {}".format(cluster, disease, size))
        if size >1:
            with col32:
                #Option to filter for samples that are Resistant, Non-resistant or all
                responses = st.radio("",
                                 ['All', 'Only Non-resistant samples', "Only Resistant samples"],
                                 key="samples-resp")
            with col33:
                #Option to select number of top features to display
                n_features = st.number_input('Number of features', value=5)
            transc = feat_inter.filterSamples(transcripts, responses)
            sample = st.selectbox("Select a sample", transc, key="sample")
            # Option to select the local MA method to use
            local_model = st.selectbox("Select a local interpretable method", ["LIME", "SHAP"], key="local_model")
            if local_model == "LIME":
                st.write("#### LIME for individual predictions")
                st.write(" ")
                st.write("Local interpretable model-agnostic explanations (LIME) is a technique to explain individual predictions "
                         "made by a ML model. For this, an explanation is obtained by approximating the main model with a more interpretable one."
                         " The interpretable model is trained on small perturbations of the original observation to provide a good local approximation."
                         "For a selected sample, LIME results will show the features that have the strongest (positive or negative) impact on the "
                         "sample prediction."
                         )
                # Option to show plot or raw data in a table
                RawD = st.checkbox("Show raw data", key="rd-LIME")
                limedf = feat_inter.limeInterpretation(sample,n_features,RawD)
                st.caption("The x-axis represents the feature's importance on the model's prediction. The y-axis represents the features with highest importance in decreasing order."
                    "Green values: positive effect on prediction. Red values: negative effect on prediction. ")
                st.write(" ")
                st.write("##### Further information about the features displayed: ")
                feature = searchTranscripts(limedf['Feature'].tolist())
            elif local_model == "SHAP":
                st.write("  ")
                st.write("#### SHAP")
                st.write(" ")
                st.write("SHAP (SHapley Additive exPlanations) is a technique that explains the prediction of an observation by "
                         "computing the contribution of each feature to the prediction. It is based on Shapley values from game theory,"
                         " as it uses fair allocation results from cooperative game to allocate credit for a model's output among its input features."
                         )
                st.write("  ")
                st.write("##### Model's prediction")
                pred, true_value = shapFI.SHAP_pred(sample)
                st.write("**Model's predicted AAC response is:** {} ".format(pred))
                st.write("**True AAC response is:** {} ".format(true_value))
                st.write(" ")
                st.write("##### Decision plot for {}".format(sample))
                st.write(
                    "The SHAP decision plot shows the relationship between the {} most influential features/transcripts and the model's prediction for the "
                    "individual cell line {}. "
                    "They are a linear representation of SHAP values.".format(n_features, sample))
                st.write("")
                #Option to show plot or raw data in a table
                RawD = st.checkbox("Show raw data", key="rd-SHAP-dec")
                transcripts = shapFI.SHAP_decision(sample, explainer, values, n_features, RawD)
                st.caption(
                    "The x-axis represents the model's AAC prediction. The y-axis represents the model's features. "
                    "The grey vertical line represents the base value. Coloured line is the prediction and how each feature impacts it."
                    " Bracket values are the real features values for the chosen sample.")
                st.write("##### Further information about the features displayed: ")
                feature = searchTranscripts(transcripts)

