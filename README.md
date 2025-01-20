# MetabSpace

This MetabSpace directory contains codes, data, and instructions for using the "chemical characteristics vector" approach for capturing the chemical space of biomes and mapping them using LC-MS/MS data. 
The main article related to MetabSpace is "*insertwhenarticleout*"

Data for the article is in Zenodo: "*[Add zenodo doi](https://doi.org/10.5281/zenodo.14506250)*"


## General description 

The following workflow describes the general idea or the whole metabolomics data analysis approach, with Step 3-4 including the "chemical characteristics vector" part. 

**Step 1** is LC-MS/MS peak extraction which can be done using software like MZmine, MS-DIAL or R-based package patRoon.

**Step 2** is molecular fingerprint and compound class prediction using SIRIUS software.

**Step 3** is our developed approach to describe the ratio of compounds in the sample with specific chemical moiety.

**Step 4** illustrates how with this approach we are now able to compare the chemical space of compounds more efficiently. 

![github_metabspace](https://github.com/user-attachments/assets/9ce45b7c-b3a6-49c2-9e64-d5f5a8a82c86)


## Code

CCV_article.R  -> Codes for gathering SIRIUS chemical characteristics and calculating averaged CCVs.
CCV_article_figures.R  -> Code for creating Figures for the article and re-analyzing the data. Uses data from Zenodo.

Functions_SIRIUS_DataAnalysis.R -> gathers functions to get data from SIRIUS calculations. Getting annotation tables and confidence scores into one table. 

## Data

Most data is available in Zenodo. For additional data, please contact Pilleriin Peets (pilleriin.peets@gmail.com, pilleriin.peets@ut.ee)

## Literature

For the article Peets et al. 202*, metabolomics data from the Earth Microbiome Project was used:

*Shaffer, J.P., Nothias, LF., Thompson, L.R. et al. Standardized multi-omics of Earth’s microbiomes reveals microbial and metabolite diversity. Nat Microbiol 7, 2128–2150 (2022). https://doi.org/10.1038/s41564-022-01266-x*
