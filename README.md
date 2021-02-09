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
  `sudo apt update`  
  `sudo apt install nucperiod`  

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
Easy installation of nucperiod can occur through the ppa at <https://launchpad.net/~ben-morledge-hampton/+archive/ubuntu/nucperiod>  
To install through this ppa, run the following commands:  
  `sudo add-apt-repository ppa:ben-morledge-hampton/nucperiod`  
  `sudo apt update`  
  `sudo apt install nucperiod`  

Currently, this installation method for nucperiod is only available on Ubuntu version 20.04, Focal Fossa, due to a dependency on a Python install of at least version 3.7.  
However, installation on other linux distributions is certainly possible through manual installation of the Python and R packages provided in this repository.  
If you believe a specific linux distribution should be supported by the ppa, but isn't, please contact me at b.morledge-hampton@wsu.edu

***
## Input Files and Formats
#### Directory Structure
All data files should be stored within the nucperiod_data directory.  

Genome data should be stored in the "\_\_external_data" directory under a sub-directory with the same name as the corresponding genome fasta file.  
e.g. "hg19.fa" should be stored in the "nucperiod_data\\\_\_external_data\\hg19" directory.  

Nucleosome positioning data should be stored in a sub-directory under the corresponding genome directory and should be named after the corresponding nucleosome positioning file.  
e.g. "MNase_nuc_pos.bed" should be stored in the "nucperiod_data\\\_\_external_data\\hg19\\MNase_nuc_pos" directory.  

Each individual mutation input file should be stored in its own directory under the "nucperiod_data" directory.  Nested directories are allowed.  
nucperiod populates these directories with all other files generated during the analysis.  

#### Genome Data
All genome information should be given in standard fasta format with chromosome identifiers as headers.

#### Nucleosome Positioning Data
All nucleosome positioning data 

#### Mutation Data
nucperiod supports two primary input formats for mutation data:  
First, data downloaded directly from the [ICGC data portal](https://dcc.icgc.org/releases) can be easily parsed using the following terminal command:  
`nucperiod parseICGC`  

Data from any other format should be converted to the specialized bed format recognized by nucperiod.  
This format is a variation on the standardized bed format and contains 6-7 tab separated data columns (with the 7th being optional).  
The columns should be formatted as follows:  
##### Column 1
Chromosome identifier.  (e.g. "chr1")  Should match the identifiers used in the corresponding genome fasta file.
##### Column 2
0 based mutation start position
##### Column 3
1 based mutation end position
##### Column 4
The base(s) in the reference genome at this position.  
If set to ".", the base(s) will be auto-acquired using the genome fasta file.  
Use the "\*" character to indicate an insertion between the two bases given in columns 2 and 3.  
##### Column 5
The base(s) that the position(s) were mutated to.  
Use the "\*" character indicates a deletion of the base(s) given in columns 2 and 3  
Use the string "OTHER" to indicate any other lesion or feature  
##### Column 6
The strand the mutation/alteration occurred in.  Single-base substitution mutations are flipped if necessary so that they occur in the pyrimidine-containing strand.  
If set to ".", the strand is determined from the genome file, if possible (not an insertion).
##### Column 7
The chort the tumor belongs to.  e.g. a donor ID or tumor type.  
This column is technically optional but is required for stratifying data in future steps.  
If any cohort designations are given, ALL entries must have designations.  
Use the "." character in this column to used to avoid assigning an entry to another cohort without breaking the above rule.

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
