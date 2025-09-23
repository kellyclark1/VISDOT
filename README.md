# NeuroTri2-VISDOT

This repository contains the data and code for [NeuroTri2-VISDOT](bhojlab.shinyapps.io/NeuroTri2-VISDOT).

## Files

#### extract_Bhaduri_data_summary.R

This R file includes the code needed to extract the summary of the publicly available data from [Bhaduri et al (2021)](https://www.nature.com/articles/s41586-021-03910-8). The summarized data file produced with this script is all_genes_second_trimester_sc_data.csv, and it is held in the ShinyApp.

#### all_genes_second_trimester_sc_data.csv.zip

This zip file contains the data summary extracted from the publicly available data from [Bhaduri et al (2021)](https://www.nature.com/articles/s41586-021-03910-8) using the above R code.

#### app.R

This R file contains the code needed to produce the ShinyApp.

#### NeuroTri2-VISDOT_documentation.pdf

This file contains the user guide for the online ShinyApp. For most users, this is what you will use. <br>

## How to use

If you are interested in using the ShinyApp to plot your genes of interest, see this site: [NeuroTri2-VISDOT](bhojlab.shinyapps.io/NeuroTri2-VISDOT).

If you want a guide on how to use the ShinyApp, see NeuroTri2-VISDOT_documentation.pdf.

If you are interested in using our code to extract a different data summary from [Bhaduri et al (2021)](https://www.nature.com/articles/s41586-021-03910-8) or from a different dataset, you can simply download and modify extract_Bhaduri_data_summary.R. Note that this single cell dataset (and others) is VERY large- you will likely need cloud computing/high-performance computing.

If you want to run VISDOT locally or use our ShinyApp code for your own project, you'll want to set up a directory for the code and data, which we'll call "visdot." Your folders should be organized as below:

    visdot

        app.R

        data

            all_genes_second_trimester_sc_data.csv

Note that you will need to unzip the csv data file. We recommend opening the script in RStudio to run.


## Cite

If you use this code or data, please cite:<br>

Bhaduri, A., Sandoval-Espinosa, C., Otero-Garcia, M. et al. An atlas of cortical arealization identifies dynamic molecular signatures. Nature 598, 200–204 (2021). https://doi.org/10.1038/s41586-021-03910-8 <br>

Clark KJ, Lubin EE, Gonzalez EM, Sangree AK, Layo-Carris DE, Durham EL, Ahrens-Nicklas RC, Nomakuchi TT, Bhoj EJ. NeuroTri2-VISDOT: An open-access tool to harness the power of second trimester human single cell data to inform models of Mendelian neurodevelopmental disorders. bioRxiv [Preprint]. 2024 Feb 4:2024.02.01.578438. doi: 10.1101/2024.02.01.578438. PMID: 38352329; PMCID: PMC10862881.
