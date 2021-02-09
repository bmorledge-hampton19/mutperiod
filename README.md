# nucperiod
##### A Hybrid Python and R toolset for characterizing nucleosome mutational periodicities.
***
## Table of Contents
1. [Quickstart Guide](#quickstart-guide)
2. [Installation Guide](#installation-guide)
3. [Input Files and Formats](#input-files-and-formats)
4. [The Primary Data Pipeline](#the-primary-data-pipeline)
5. [Interpreting Results](#interpreting-results)
6. [A Representative Example](#a-representative-example)
7. [Acknowledgements](#acknowledgements)
***
## Quickstart Guide
#### 1. Install nucperiod 
Install nucperiod through apt using the following two commands (available on Ubuntu version 20.04, Focal Fossa):  
  `sudo add-apt-repository ppa:ben-morledge-hampton/nucperiod`  
  `sudo apt-get update`

#### 2. Set up the Data Directory
After installing, run the following command:  
  `nucperiod parseICGC`  
This should open up a dialog to choose a directory to store data files in.  Choose a directory.  
After choosing a directory, you should quit out of the following dialog to obtain the necessary data to run your analysis  

#### 3. Obtain Genome and Nucleosome Positioning Data
From the following link, download the zipped hg19 directory:  
<https://notAValidLinkYet.com>  
Navigate to the directory you chose in step 2 and unzip the contents of the hg19 file into the "nucperiod_data/\_\_external_data" directory.  
*Note:  Alternatively, you may use a different genome or nucleosome positioning file, but you must make sure the latter is formatted correctly, as detailed in [Section 3](#input-files-and-formats)*

#### 4. Obtain Mutation Data
Go to the [ICGC data portal](https://dcc.icgc.org/releases) to obtain mutation data for use in nucperiod.
Download any "simple_somatic_mutation" file with whole genome sequencing data.  
Place the gzipped file into a new directory (other than the "\_\_external_data" directory) within the "nucperiod_data" directory
*Careful:  nucperiod only maps mutations originating from whole genome sequencing data.  Exome sequencing data will not be carried through the pipeline, potentially resulting in blank output files.*  

*Note:  Alternatively, you may use any bed formatted mutation data for analysis, but you must may need to alter the data format slightly to be recognized by nucperiod as a CustomBed formatted file, as detailed in [Section 3](#input-files-and-formats)*

#### 5. Parse Input Data
If using data from ICGC, run the following command:  
  `nucperiod parseICGC`  
Otherwise, if you are using custom bed input, run this command:  
  `nucperiod parseBed`  
Fill out the resulting dialogs using the files obtained in steps 3 and 4.

#### 6. Perform Periodicity Analysis
Run the following command:
  `nucperiod mainPipeline`  
Select the directory you created within the "nucperiod_data" directory in step 4.  
Select the desired normalization method and search radius.  

After main pipeline has finished running, run this next command:  
  `nucperiod periodicityAnalysis`  
Once again, select the directory you created to find the nucleosome mutation counts files.  
Select an output file to store the results of the analysis.  
In the "Main Group" portion of the dialog, select the normalization and radius options that you used when running the main pipeline.  
All other options may be left unaltered.

*Note:  Both .rda and .tsv formats are supported as output, but only the .rda format supports figure generation in the next step.*

#### 7. Visualizing Results
To visualize results, run the following command:  
  `nucperiod generateFigures`  
Select the files you want to visualize and the output format for the figure(s).  

#### 8. Other Features
Nucperiod also supports stratification of input data by various conditions and comparison of the periodicities of each of the resulting cohorts.  
For more information, see the sections below.

***
## Installation Guide
some text  
some text  
some text  
some text  
some text  
some text  
some text  
some text  
some text  
***
## Input Files and Formats
some text  
some text  
some text  
some text  
some text  
some text  
some text  
some text  
some text  
***
## The Primary Data Pipeline
some text  
some text  
some text  
some text  
some text  
some text  
some text  
some text  
some text  
***
## Interpreting Results
some text  
some text  
some text  
some text  
some text  
some text  
some text  
some text  
some text  
***
## A Representative Example
some text  
some text  
some text  
some text  
some text  
some text  
some text  
some text  
some text  
***
## Acknowledgements
some text  
some text  
some text  
some text  
some text  
some text  
some text  
some text  
some text  
