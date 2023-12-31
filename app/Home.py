import tempfile
from common import add_logo,sideBar,readPickle
import webbrowser
from interpretation import *
from streamlit_option_menu import option_menu
from explorer import Data
from analysis import Analysis

#Opens Wikipedia page of drug
def openWikipedia(drug):
    webpage_link = "https://en.wikipedia.org/wiki/"+drug
    webbrowser.open(webpage_link)

# Opens PubMed page of drug
def openPubMed(drug):
    webpage_link = "https://pubmed.ncbi.nlm.nih.gov/?term={}" + drug
    webbrowser.open(webpage_link)

#Opens DrugBank page of the drug
def openDrugBank(num):
    webpage_link = "https://go.drugbank.com/drugs/" + num
    webbrowser.open(webpage_link)

st.set_page_config(
    page_title="EMERY",
    page_icon="bi bi-basket",
    layout="wide"
)

add_logo()
sideBar()
st.header("Home")
st.write("---")
menu = option_menu(None, ["Overview", "Data Upload", "Trial information"],
    icons=["bi bi-basket", 'cloud-upload', 'capsule'],
    menu_icon="cast", default_index=0, orientation="horizontal")
files = "https://gla-my.sharepoint.com/personal/ronan_daly_glasgow_ac_uk/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fronan%5Fdaly%5Fglasgow%5Fac%5Fuk%2FDocuments%2FJoe%20Files%2FpyBasket&ga=1"

#Overview subpage: general information about the app and basic steps to use it
if menu == "Overview":
    st.subheader("Welcome to EMERY 👋")
    st.write(" ")
    st.write("**EMERY** (**E**xplainable **M**achine l**E**arning fo**R** p**Y**basket) is a user-friendly interactive platform that allows the visualisation and exploration"
             " of results obtained from the pyBasket pipeline. The main goal of EMERY is to aid in the investigation of"
             " the omics profiles of the samples or patients involved in clinical trials and gain biological insights"
             " into the resistance patterns showed to a drug by different subgroups of patients. Also, interpretable Machine Learning"
             " methods can be used to identify features (genes) that have an important impact on the drug resistance rates."
             " EMERY can be used to interpret the results of the pyBasket pipeline and enhance the design of basket trials by identifying subgroups "
             "of patients within a disease type that are more likely to respond to a treatment.")
    st.write(" ")
    st.write("##### Basket trials")
    st.write("Basket trials are a type of early phase clinical trial design that tests how well a new drug or treatment works in"
             " patients that are grouped based on the same molecular mutation or biomarker within a type of disease."
             " This innovative approach is more effective and have a reduced cost, as it requires fewer patients to be enrolled "
             "for a study to be carried out and allows "
             "intractable or rare cancer patients to also be included in general clinical trials."
             " As they allow for multiple diseases to be studied, it facilitates faster drug development for specific "
             "subpopulations of patients and allows for a more personalised treatment approach. ")
    st.write(" ")
    st.write("##### pyBasket")
    url = "https://glasgowcompbio.github.io/pyBasket/"
    st.write("pyBasket is a two-stage approach previously developed by other researchers that incorporates omics data into clinical trial response prediction from patients."
             " In the first stage, patients are clustered based on their omics profile using K-means clustering. "
             " These clusters assignments are used within a hierarchical Bayesian model along with response rates from different types of disease to "
             " estimate overall response rates and predict interaction terms between each type of disease and cluster. Further information about the pyBasket model can be found here: [link](%s)" % url)
    st.write("---")
    with st.expander("##### Basic steps to use EMERY"):
        st.write(" ")
        col11, col12 = st.columns((2,2))
        with col11:
            st.write("##### :one: Upload the data")
            st.write("Navigate to the _Data Upload_ subpage in this same page. Select a pickle file (.p) with results obtained from the pyBasket pipeline to be uploaded. You can get some example results"
                     " file from the pyBasket pipeline's [OneDrive folder](%s)." %files)
            st.write(" ")
            st.write("##### :two: Check file information")
            st.write(
                "Navigate to the _Trial Information_ subpage in this same page to find information about the file uploaded.")
            st.write(" ")
            st.write("##### :three: Check drug information")
            st.write(
                "In the same _Trial Information_ subpage, further information about the drug tested in the trial being analysed can also be found. "
                "This includes direct links to the drug's pages in sites like DrugBank, Wikipedia or PubMed. ")
            st.write(" ")
        with col12:
            st.write("##### :four: Explore the data")
            st.write(
                "Navigate to the _Data Exploration_ page located in the left sidebar. Find more general information about the samples, the predicted treatment response rates and "
                "perform Dimensionality reduction, find the Prototypical sample of each group or do Differential Expression Analysis.")
            st.write(" ")
            st.write("##### :five: Select and analyse a disease-cluster interaction")
            st.write("In the left sidebar, navigate to the _Disease-Cluster_ _Interactions_ page to get an overview of all possible cluster"
                     " and disease combinations. Select a cluster number and a disease type in the sidebar. Samples that fall in this interaction will be selected. "
                     "Navigate to the _Selected interaction_ subpage to further explore results from the samples included in the selected disease-cluster interaction.")
            st.write(" ")
            st.write("##### :six: Use interpretable ML methods ")
            st.write(
                "Navigate to the _Feature_ _Importance_ _Analysis_ subpage located in the left sidebar to use interpretable Machine Learning methods and explore important transcripts.")
            st.write(" ")

#Upload data subpage
if menu == "Data Upload":
    st.subheader("Please upload your data")
    st.write("Click the box to select a pickle file with results from the pyBasket pipeline to be uploaded. Other results files can be sourced from the pyBasket project's [OneDrive folder](%files).")
    input_file = st.file_uploader('Upload your data file (pickle format .p)', type='p')
    file_name = ""
    if input_file is not None:
        with tempfile.NamedTemporaryFile(suffix=".py") as tmp_file:
            tmp_file.write(input_file.getvalue())
            tmp_file.flush()
            file_name = tmp_file.name
            st.success('The file was successfully uploaded!', icon="✅")
            save_data = readPickle(tmp_file.name)
            data = Data(save_data,input_file.name)
            analysis = Analysis(save_data,input_file.name)
            data.setData()
            st.session_state["data"] = data
            st.session_state["saved data"] = save_data
            st.session_state["File Name"] = input_file.name
            st.session_state["Analysis"] = analysis
    else:
        st.write(file_name)
    if "File Name" in st.session_state:
        st.success('Analysis ready', icon="✅")
        st.info("Current file: {}".format(st.session_state['File Name']))

#File and drug information subpage
if menu == "Trial information":
    st.subheader("File information")
    st.write("Once a file of results is uploaded, general information about the file is displayed here.")
    if "data" in st.session_state:
        st.session_state["data"].fileInfo()
        st.write(" ")
        st.subheader("Drug information")
        name = st.session_state["File Name"].split('_')
        drug = name[2]
        st.write("#### The drug used in this clinical trial was: **{}**".format(drug))
        st.write("Further information about the drug used in the pyBasket analysis.")
        accession_num = "0"
        if drug == "Erlotinib":
            accession_num = "DB00530"
        elif drug == "Docetaxel":
            accession_num = "DB01248"
        st.write("#### _DrugBank_")
        st.write("_DrugBank_ is a freely accesible online database containing information on drugs and drugs targets. It combines chemical,"
                 " pharmacological and pharmaceutical data with information about drug target (sequence, structure, pathway, etc.)"
                 " Further pharmacological and chemical information about {} can be found in the link to _DrugBank_ below.".format(drug))
        st.button('Open DrugBank', on_click=openDrugBank, args=(accession_num,))
        st.write("#### _PubMed_")
        st.write("_PubMed_ is a free full-text archive of biomedical and life sciences print and electronic journal literature,"
                 " and supports current biomedical and healthcare research and practice. Other research and journal literature related to "
                 "{} can be found in the link to _PubMed_ below.".format(drug))
        st.button('Open PubMed',on_click=openPubMed, args=(drug,))

        st.write("#### _Wikipedia_")
        st.write("_Wikipedia_ is a free online accessible encyclopedia where more general non-technical information about {} can be found.".format(drug))
        st.button('Open Wikipedia', on_click=openWikipedia, args=(drug,))
    else:
        st.warning("No data found. Please upload results in the 'Data Upload' tab")
