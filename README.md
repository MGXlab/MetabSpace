# MetabSpace

This MetabSpace directory contains codes, data, and instructions for using the "chemical characteristics vector" approach for capturing the chemical space of biomes and mapping them using LC-MS/MS data. 
The main article related to MetabSpace is "*insertwhenarticleout*"


## General description 

The following workflow describes the general idea or the whole metabolomics data analysis approach, with Step 3-4 including the "chemical characteristics vector" part. 

**Step 1** is LC-MS/MS peak extraction which can be done using software like MZmine, MS-DIAL or R-based package patRoon.

**Step 2** is molecular fingerprint and compound class prediction using SIRIUS software.

**Step 3** is our developed approach to describe the ratio of compounds in the sample with specific chemical moiety.

**Step 4** illustrates how with this approach we are now able to compare the chemical space of compounds more efficiently. 

![github_metabspace](https://github.com/user-attachments/assets/9ce45b7c-b3a6-49c2-9e64-d5f5a8a82c86)


## Code

CCV_article.R  -> Codes for manuscript "". 
Functions_SIRIUS_DataAnalysis.R -> gathers functions to get data from SIRIUS calculations. Getting annotation tables and confidence scores into one table. 

## Data

## Literature
