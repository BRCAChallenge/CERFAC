# CERFAC
## Cloud Enabled, Rigorous, Functional Assay Calibration

Multiplexed assays of variant effect (MAVEs) and computational tools have the potential to help classify many variants that are discovered during clinical genetic sequencing. However, predictions from assays and tools are not direct measures of human health and require careful calibration and validation prior to clinical use.

Calibration of functional assays results in thresholds above or below which variants are predicted to be associated with disease or not. Calibration often relies on the ability of the assay to correctly classify variants of known effect (for example, to correctly classify variants which have been classified in ClinVar). However, it is known that this method introduces circularities and inflated assessments of performance; for example, existing ClinVar classifications might have been informed by the results of similar functional assays. Rather, if the calibration is validated on sets of clinical and control observations, such as case-control data, it is likely that many of the observed variants will not yet be classified, and are thus less subject to problems with circularities; this approach provides a more objective means to assess the performance of the assay on novel variants.

This strategy relies on identifying a group of variants that are known to be pathogenic and calculating the effect size (in this case an odds ratio, OR) of these variants in the case and control data. This provides a “standard candle” by which other groups of variants, such as those predicted to be associated with disease by the assay, can be evaluated.

Groups of variants with an OR close to or equal to 1 are consistent with benign variants (equally likely to be observed in cases and controls). Groups of variants with an elevated OR close to that of the “standard candle” would be consistent with variants that are associated with the disease in question (more likely to be observed in cases than in controls). Groups of variants that fall between an OR of 1 and the OR of the “standard candle” group will need to be carefully interpreted within the clinical context of the disease.

This workspace provides a resource for users to validate functional assay calibrations. It requires the following user input:

1. Functional assay scores and corresponding calibration thresholds
2. Clinical observational data for gene/disease of interest
3. Clinically informed ORs thresholds for gene/disease of interest


# Overview

This page contains:

1. A summary of the workflows and analyses available in this workspace, and estimates for the cloud costs required for their execution.  In the parlance of Terra documentation, “workflows” refer to scripts with ordered steps written in WDL (Workflow Description Language), while "Analyses" refer to Jupyter notebooks.

2. Step-by-step instructions for their execution

This page contains the step-by-step instructions first, followed by the technical summary information.

# Step-by-step instructions
## Before you begin: Create your own editable copy (clone) of this WORKSPACE


In order to run the workflows/notebooks in this tutorial, you need to have edit and execute permissions in the workspace. Cloning the workspace makes an exact copy with you as owner.
In the upper right of this page, click the circle with three dots, then click Clone. 

1. Click on the round circle with three dots in the upper right corner of this page:

![](https://github.com/BRCAChallenge/CERFAC/blob/main/doc_images/click_clone_upper_right.png?raw=true)

2. Select "Clone" from the dropdown menu

3. Rename your new workspace something memorable

Note: It may help to write down or memorize the name of your workspace

4. Choose your billing project (could be free credit) from the dropdown menu

5. Click the "Clone Workspace" button to make your own copy


Let's walk through an example of how to run these workflows. We'll be using BRCA1 as our example. 

## Step 1: Set up your Data Table with your gene of interest. 
Within this workspace, we have provided a sample table to keep track of our calibrations and to organize the results. The name of the table is “Configurations”.

Navigate to the DATA tab. You may wish to open it in a new tab and refer back to this documentation. 

Under TABLES, click on **Configurations**. If the workspace was cloned correctly, the Configurations data table will contain one row of workflow configuration data, with a sample ID of "BRCA1_template".  The steps here will show you how to reproduce the contents of this row.

![](https://github.com/BRCAChallenge/CERFAC/blob/main/doc_images/1_data_tab_configurations_boxed.png)


In this example, we are going to be adding a row for the gene BRCA1 to our sample data table. 

At the top left of the data table, click Edit > Add row. 

![](https://github.com/BRCAChallenge/CERFAC/blob/main/doc_images/Edit_add_row.png?raw=true)


Fill out the pop-up: The first box will be our sample name. In this case, we will add "BRCA1_run". Additonally, we will click on the > symbol next to GENE_NAME and type "BRCA1". It is important that you type the correct HGNC gene name. 

![](https://github.com/BRCAChallenge/CERFAC/blob/main/doc_images/Add_new_row_BRCA1.png?raw=true)

Leave the other inputs blank. 

Click "Add" in the blue box at the bottom right. 

## Step 2: Upload clinical data file and functional assay data file to workspace data. 

It is somewhat complicated to make a file available to a workflow in Terra, but doable. 

**First**, we must upload the file to our workspace Google Cloud Bucket in Terra...

Still in the DATA tab, navigate to OTHER DATA > Files. Click on Files. 

![](https://github.com/BRCAChallenge/CERFAC/blob/main/doc_images/click_on_files_data_section.png?raw=true)

A list of files in the workspace should appear in blue text (do not worry if you have different files from the screenshot). Click the blue Upload button in the top right corner. It will prompt you to add files locally from your computer. Make sure you add BOTH your clinical/case data file AND your functional assay data file!

![](https://github.com/BRCAChallenge/CERFAC/blob/main/doc_images/upload_to_Files_green_arrow.png?raw=true)

...**and then** we must add those files to the data table. 

After uploading a file to Files, it will appear in the list of files. If you mouse over the filename, a blue clipboard icon should appear. Click on that to copy the file's URL to your computer's clipboard. Keep in mind whether you are copying the name of the clinical or functional assay data file. 

![](https://github.com/BRCAChallenge/CERFAC/blob/main/doc_images/copy_file_url_to_clipboard_files_arrow.png?raw=true)


Navigate back to the sample data table. In the BRCA1 row, we mouse over to the "functional_assay_file" column. When you mouse over the box, a small blue pencil should appear. Click on it to edit, and paste the file name you copied earlier into the pop-up box. Click save. 

![](https://github.com/BRCAChallenge/CERFAC/blob/main/doc_images/editing_value_icon_add_file_name.png?raw=true)

Now we've got to do the other file. Go back to the Files section, and copy the URL for the clinical data file, then paste it in the relevant row and column of the sample Data Table. 

NOTE: Make sure you are adding the correct files to the correct column! Example: Your functional assay data file should be added to the "functional_assay_file" column, not the "clinical_data_file" column. 

## Step 3: Generate ClinVar calibation variants by running the workflow. 

Now that the Data Table has the necessary files, we can run the workflows. We're going to continue to use BRCA1 as an example. 

It is easiest to start the workflows directly from the sample Data Table. 

Check the box next to the row you wish to run, in this case BRCA1. Then click Open With...

![](https://github.com/BRCAChallenge/CERFAC/blob/main/doc_images/check_box_inBRCA1_row_open_with.png?raw=true)

In the pop-up that appears, click Workflow. 

![](https://github.com/BRCAChallenge/CERFAC/blob/main/doc_images/open_with_wdlworkflow.png?raw=true)

Then click 1-get_clinvar_variants.


![](https://github.com/BRCAChallenge/CERFAC/blob/main/doc_images/workflows_select_get_cvvars.png?raw=true)

The Workflow configuration page should appear. 

The user options are boxed in blue. 

![](https://github.com/BRCAChallenge/CERFAC/blob/main/doc_images/workflow_clinvar_variants_top_page.png?raw=true)


The only variable inputs that need to be configured is the GENE_NAME. Click in the input value box, and scroll to "this.GENE_NAME".  That tells the workflow to look in the GENE_NAME column of our data table. This should only need to be configured the first time you run a workflow. You can leave the other inputs blank. 

![](https://github.com/BRCAChallenge/CERFAC/blob/main/doc_images/clinvar_workflow_thisgenename.png?raw=true)


Then click the **save** button on the right. Then click Run Analysis. You can add a comment in the pop-up if you like, and then click launch. The workflow should begin running!

![](https://github.com/BRCAChallenge/CERFAC/blob/main/doc_images/Workflow_run_analysis.png?raw=true)



## Step 4: Generate gnomAD calibration variants by running the 2-get_gnomad_variants workflow. 

Repeat step 3, this time selecting the "2-get_gnomad_variants" workflow. Because of the size of the gnomAD database, this takes the longest to run. The processing time, and the memory needed to run the workflow, will depend in part on the length of the gene of interest. The workflow is optimized to allocate the minimum amount of memory needed based on the length of the gene. 



## Step 5: Merge calibration variants using the 3-merge_clinical_functional_data workflow.

Repeat step 3, this time selecting the "3-merge_clinical_functional_data" workflow.
This merges the calibration variants, the functional scores file you provide, and the clinical observational data file you provide. 


Because variant nomenclature varies widely, this workflow interacts with an API to obtain the VRS IDs of  the variants within each file, and then uses these VRS IDs to perform a merge. The files can be in any of the following formats: TSV, CSV, or tab separated TXT. 
However, only **HGVS variants** (cDNA/coding or genomic) *with valid accessions* OR **gnomAD VCF** style variants are supported. **The column with variant IDs must be the first column in your files.** The workflow will fail if that is not the case. 

 ## Step 6: Calibrating the assay: R-Jupyter Notebook Analysis
To calibrate the variants, use the Jupyter notebook analysis in the workspace. 

Make sure the environment is set to use R!

To perform this analysis you will need:
- The output (file) of the three workflows above, which combines the calibration variants the workflow generates for you, the functional assay file you provide, and the clinical observational data file you provide. 
- Functional threshold score cutoffs for the functional assay
- Odds ratio threshold
- The population frequency of the relevant disease
- Total control and case counts (from your clinical observational data)



## Troubleshooting the workflows
The workflow is optimized to allocate the minimum amount of memory needed based on the length of the gene. However, if the workflow fails due to an out of memory error, you can manually input a higher amount of memory. 

Checking the pre-emptible workflow option block lowers the costs of running a workflow, but there is a chance that your workflow will be preempted-and fail- without giving you a clear error message. 


Re-running the workflow with more GB > shouldn't be necessary

### Formatting errors that cause the workflow to abort
* The first line of the files must be column names. The column names can be whatever you want, but they must be there. 
* The HGVS variants must contain an accession. A gene name is not sufficient. 
* A colon (:) must separate the accession from the "g." or "c.". If it is missing in the first variant, the workflow will abort. 
* gnomAD VCF IDs must have dashes separating the fields, not colons or underscores. 

 ### Errors that cause the workflow to skip to the next variant
These errors won't interrupt the workflow, but they will skip to the next variant without generating a VRS_ID. 
* An HGVS variant contains "inv" for inversion. The API doesn't currently support inversions. 
* An HGVS variant contains bases (ACTG) after "dup" or "del". Some organizations use this format, but it is not supported by the API. 
* Variants containing bases other than ACTG. 
* HGVS coding variants in introns or upstream of a gene (containing + or -). Intronic HGVS coding variants are not supported (they are okay in VCF format). 
* HGVS variants containing asterisks. 

Variants that are skipped will be printed to stdout in the Google Cloud environment with the other workflow outputs. 


# Placeholder




# Workflows


In the parlance of Terra documentation, “workflows” refer to scripts with ordered steps written in WDL (Workflow Description Language), while "Analyses" refer to Jupyter notebooks.

This workspace contains 3 WDL workflows and 1 Jupyter notebook analysis.

The 3 workflows are used to generate calibration variants for the gene of interest from gnomAD and ClinVar variants.

The 3rd workflow requires 2 user-provided files:

- Functional assay data for the gene of interest including scores.

- Clinical observation data pertaining to the gene of interest.

The output will be a tsv file that the user can download and edit.

The Jupyter notebook contains the actual calibration.

To perform Jupyter notebook analysis you will need:

- The output (file) of the three workflows above, which combines the calibration variants the workflow generates for you, the functional assay file you provide, and the clinical observational data file you provide.

- Functional threshold score cutoffs for the functional assay

- Odds ratio threshold

- The population frequency of the relevant disease

- Total control and case counts (from your clinical observational data)


## 1-get_clinvar_variants
**What does it do?**

This WDL pipeline finds ClinVar variants for your gene of interest using the NCBI's [Entrez Direct (EDirect)](https://www.ncbi.nlm.nih.gov/books/NBK179288/) version of the [E-utilities program](https://www.ncbi.nlm.nih.gov/books/NBK25497/). It processes this data and returns  classifications, associated diseases, comments, functional assays results, variant effects, and submission information. ...

**What data does it require as input?**

It requires only the [HGNC gene name](https://www.genenames.org/tools/search/#!/?query=&rows=20&start=0&filter=document_type:gene) (gene symbol). Make sure you are using the most recent HGNC approved name. If unsure, you can check the gene name [here](https://www.genenames.org/tools/multi-symbol-checker/). 

**What does it return as output?**

It outputs a csv file containing rows for each submitter per variant. In other words, if a variant has multiple submissions in ClinVar, it will have multiple rows. (User feedback welcome.)

**Configuration notes**

The workflow is written in WDL1.0, and the tasks are written in a combination of bash and python. 
See [Terra documentation on how to run a  workflow](https://support.terra.bio/hc/en-us/articles/360036379771-Overview-Running-workflows-in-Terra). 


**Time and cost estimates**


For more information about controlling Cloud costs, see this article.

## 2-get_gnomad_variants
**What does it do?**

This WDL pipeline accesses the gnomAD database using the Hail python package and selects the portion of the database pertaining to your gene of interest. It returns variant information from both the exome and genome databases including frequency of variants by ancestry. 
It also merges variants with the ClinVar file. 

**What data does it require as input?**

It requires the gene name and the ClinVar variants file which is an output of 1-get_clinvar_variants. 

**What does it return as output?**
A csv file with the gnomAD and ClinVar variants combined for the gene of interest. Variant IDs are in the HGVS DNA coding variant form based on the MANE Select transcript. 


**Configuration notes**

The workflow is written in WDL1.0, and the tasks are written in a combination of bash and python. 


**Time and cost estimates**
Below is an example of the time and cost for running the workflow...

This workflow is the most computationally and memory-intensive part of the pipeline. The amount of memory needed will depend on the length of the gene transcript. 



## 3-merge_clinical_functional_data

**What does it do?**

This WDL pipeline merges the calibration output file with the user-provided clinical and functional assay data files. 

**What data does it require as input?**

It requires the previous two workflows have been run successfully and the user-provided clinical and functional assay data files. These two files should have the **variant ID in HGVS DNA coding variant format as the first column**, based on the **MANE select transcript**. However do not include the transcript ID or gene name in the variant ID. 
Examples: 

c.884G>A

c.897+9T>C

c.919_920del

[Here is a link describing the HGVS format for DNA coding transcripts.](https://hgvs-nomenclature.org/stable/recommendations/general/) 

**What does it return as output?**

A csv file with all variant information merged. It will be used as input for the Jupyter notebook analyses. 

## Analyses
(jupyter notebook stuff here)

## Time and Cost Information
*Under construction*

https://support.terra.bio/hc/en-us/articles/6123082826651-Overview-Costs-and-billing-in-Terra-GCP

### Memory used: 

Each of the 5 tasks in 1-get_clinvar_variants uses 10 GB of memory, but users may optionally input a smaller or larger amount if necessary when configuring the workflow. 

Memory allocation is more complicated for the 2-get_gnomad_variants workflow. The amount of memory required is larger for longer genes. 
* The workflow starts at 15 GB or memory and automatically allocates more memory for larger genes: 

* an additional 30 GB for genes longer than 1 million bases.

* and an additional 45 GB for genes longer than 2 million bases (for a total of 90 GB). 

#### Examples: 

* BRCA1: uses about 1 GB of memory for 1-get_clinvar_variants and  uses about 7 GB of memory for the 2-get_gnomad_variants workflow. 

* SRY: a short gene, uses about 1 GB of memory for 1-get_clinvar_variants and X GB for  2-get_gnomad_variants workflow

* DMD, a gene over 2 million bases long, uses about 1 GB of memory for 1-get_clinvar_variants and X GB for  2-get_gnomad_variants workflow.

| Gene              | Gene length | Step | Time | Memory | Cost |
| :----------------: | :---------------: | :---------------: | :---------------: | :---------------: | :---------------: |
| SRY           |  827 bases   | 1-get_clinvar_variants | 15 min | 1 GB | $X |
| SRY           |  827 bases   | 2-get_gnomad_variants | 11 min | 2.5 GB | $X |
| BRCA1        | 126,032  bases   | 1-get_clinvar_variants | 33 min | 1.5 GB | $X |
| BRCA1        | 126,032 bases   | 2-get_gnomad_variants | 12 min | 7 GB | $X |
| SOX5    | 1,033,146 bases   | 1-get_clinvar_variants | ? min | ? GB | $X |
| SOX5    | 1,033,146 bases   | 2-get_gnomad_variants | ? min | ? GB | $X |
| DMD    | 2,220,166 bases   | 1-get_clinvar_variants | 12 min | 1 GB | $X |
| DMD    | 2,220,166 bases   | 2-get_gnomad_variants | 24 min | 83 GB | $X |

The cost of 3-merge_clinical_functional_data will depend on the files the user submits.


Cost also depends on whether the workflow is pre-emptible or not. Pre-emptible workflows are cheaper but may fail without warning. [See article on managing cloud costs. ](https://support.terra.bio/hc/en-us/sections/360006459511-Managing-Cloud-costs). See [this article on estimating cloud costs](https://support.terra.bio/hc/en-us/articles/360029772212-Controlling-Google-Cloud-costs-sample-use-cases). 


Authors: archeney@ucsc.edu, mcline@ucsc.edu

Date: August 8, 2025
