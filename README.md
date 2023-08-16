# EMERY: A Graphic User Interface for the pyBasket model

This is the repository for EMERY, a graphic-user interface developed for the pyBasket pipeline as part of my MSc Bioinformatics dissertation. Original pyBasket [repository](https://glasgowcompbio.github.io/pyBasket/)

**Abstract**
>Basket trials have become a popular clinical trial design as an alternative to traditional disease-specific trials. In this innovative design, patients within disease groups are arranged based on a shared molecular biomarker or alteration, facilitating the identification of the most effective targeted treatment for a subgroup of patients. As they focus on molecular profiles or biomarkers instead of a specific type of cancer, basket trials allow for the inclusion of intractable or rare cancer patients in the clinical trials, reducing the number of these needed, increasing efficiency and compatibility with a personalised medicine approach and reducing the cost. Based on this, previous researchers have developed pyBasket, a novel two-stage approach that incorporates transcriptomics data into a basket trial design to enhance response rates predictions. However, the complexity of the Machine Learning (ML) and statistical techniques that are part of the pyBasket pipeline complicates the investigation of the results obtained in a clinical setting. To address this downside, EMERY (Explainable Machine lEarning foR pYBasket) is a user-friendly interface that has been developed to explore and further investigate the omics profiles of patients and different subgroups, as well as the results from the pyBasket pipeline. The interface presents a more accessible and easier way to gain biological insights into the estimated response rates, to identify specific subgroups of patients that respond differently to a given treatment and aid in the design of enhanced basket trials. As ML methods and results are commonly not well reported and understood in a clinical context, EMERY has been equipped with several interpretable ML methods, data visualisation and statistical analyses techniques that can be applied. Besides the facilitation of interpretation and evaluation of results from the pyBasket pipeline, these are also aimed to identify interesting transcripts that are driving the treatment response predictions. The utility of EMERY was demonstrated with three applications on a dataset from the Genomics Drug Sensitivity in Cancer (GDSC) project and one of the drugs available, Erlotinib. First, a cancer group which is not normally treated with the given drug was shown to be predicted with a high probability of response. The second was the identification of two molecularly-defined subgroups of patients within the same cancer type that were predicted to respond differently to Erlotinib. Finally, the epidermal growth factor receptor (EGFR) transcript, which is specifically targeted by Erlotinib, was identified to be one of the most influential features in the response rates prediction.

**Steps to use EMERY**

***1. Download the required file(s):***

Download the necessary pyBasket pipeline results saved in a pickle format from the [OneDrive folder](https://gla-my.sharepoint.com/:f:/g/personal/ronan_daly_glasgow_ac_uk/Eod_I6-9hDtCgJ1CmKdBJCAB66sciwg58zlxDHD2fgtsMw?e=0MA2gb).

***2. Set up the virtual environment***

#### Option 1: Install dependencies using Pipenv and activate virtual environment
   1. Install pipenv: https://pipenv.readthedocs.io)
   2. Clone the repository in the desired location of your own machine: git clone 'clone link repository'
   3. In the cloned directory, run `$ pipenv install`. (If this doesn't work, use first: `$ sudo -H pip install -U pipenv`
   4. Enter virtual environment using `$ pipenv shell`.

#### Option 2: Manage dependencies using Anaconda Python:
   1. Install Anaconda Python (https://www.anaconda.com/products/individual).
   2. In the cloned Github repository, run `$ conda env create --file environment.yml`.
   3. Enter virtual environment using `$ conda activate pyBasketApp`.
 
***3. Run the app***

 1. Move to the "app" directory and from the command line run: `$ sh runApp.sh`
 2. The app will open in a new page in the web browser. Go to "Upload Data" in the Home page to upload the pickle file with pyBasket pipeline results.
