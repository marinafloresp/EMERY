# MScDiss: A Graphic User Interface for the pyBasket model

This is the repository for the pyBasket app developed as part of my MSc Bioinformatics dissertation. Original pyBasket [repository](https://glasgowcompbio.github.io/pyBasket/)

**Abstract**
>Basket trials have become a popular clinical trial design as an alternative to traditional disease-specific trials. This innovative strategy groups patients based on a shared molecular biomarker or alteration, irrespectively of the tumour’s anatomical location. As basket trials allow drugs that target specific molecular subtypes to be tested across different types of cancer, they can identify treatments that could be effective for subgroups of patients. This design is specially compatible with personalised medicine and is proven to be more cost-effective and facilitate the inclusion of rare cancers or cancers with rare genetic alterations. As basket trials are usually only focused on the response to treatment and overlook the underlying molecular information, previous researchers have developed pyBasket. This novel two-stage approach incorporates omics data into a basket trial design with cancer cell-lines to enhance response rates predictions. However, the complexity of the Machine Learning (ML) and modelling techniques that are part of the pyBasket pipeline might complicate the investigation of the results in a clinical setting. Here, a user-friendly interface to explore and further investigate the omics profiles of patients and the model’s results is presented. This represents a more accessible and easier way to gain biological insights into the estimated response rates and to identify specific groups of patients that respond differently to a given treatment. With clinical researchers with little or none ML expertise in mind, the pyBasket app has been equipped with several data visualisation techniques, statistical analyses and  interpretable ML methods. These are aimed to facilitate the report, interpretation and evaluation of results from the pipeline, as well as to identify potential interesting features that are driving the model’s predictions. With a case study, it is demonstrated how the pyBasket app can be used to explore results obtained by applying the pipeline on a dataset from the Genomics Drug Sensitivity in Cancer (GDSC) project and one of the drugs available, Erlotinib. Results enabled the identification of a basket and omics profile cluster combination that was predicted to likely respond to the treatment. The implementation of an interactive interface provides pyBasket with numerous advantages in comparison to existing alternative methods.

**Steps to use the pyBasket app**

To run the app, follow the following steps:

***1. Download the required file(s):***
1. Download a file with pyBasket model results from the "data" folder [here](https://gla-my.sharepoint.com/:f:/g/personal/ronan_daly_glasgow_ac_uk/Eod_I6-9hDtCgJ1CmKdBJCAB1LafkJ1UK3-Opp7UdQp1_Q?e=5SehkI)
2. This file is uploaded in the app.
   
***2. Set up the virtual environment***

#### Install dependencies using Pipenv and activate virtual environment
   1. Install pipenv: https://pipenv.readthedocs.io)
   2. Clone the repository in the desired location of your own machine: git clone 'clone link repository'
   3. In the cloned directory, run `$ pipenv install`. (If this doesn't work, use first: `$ sudo -H pip install -U pipenv`
   4. Enter virtual environment using `$ pipenv shell`.

#### Manage dependencies using Anaconda Python:
   1. Install Anaconda Python (https://www.anaconda.com/products/individual).
   2. In the cloned Github repo, run `$ conda env create --file environment.yml`.
   3. Enter virtual environment using `$ conda activate pyBasketApp`.
 
***3. Run the app***

 From the "app" directory, run in the command line: `$ sh runApp.sh`
