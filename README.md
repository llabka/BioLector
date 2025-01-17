# BioLector data analysis

This repository contains a script for analyzing output data from BioLector microbioreactors and a supporting package "htmf" that provides functions to facilitate this analysis (reading, normalizing, plotting, and extracting features from growth profiles). The BioLector is a high-throughput microbioreactor system in microtiter plate format used for monitoring e.g. microbial growth (biomass), pH, dissolved oxygen, and fluorescence. The tool provided in this repository allows for the efficient analysis and visualization of data generated from BioLector experiments including calculation of the growth rate in exponential phase.

**Repository Contents**

- src/: Source folder containing BL2_plot_and_mu_script.R (the R script for analyzing BioLector data) and the supporting package, htmf, providing functions necessary for data analysis.
- Data/: Folder containing example data files used for demonstration and testing purposes.
- design_file.txt: design file is a table with information on each well required to run the script
- README.md: This file, providing an overview of the repository contents and instructions for usage.

**Usage**

The output file from BioLector runs is CSV format file containing all data. By running the R script, it will process the data in the .CSV file and generate an output file "Info" with the results from the run, an output file "Growthrates" with the analysed maximum growth rate (µmax) for each strain (as listed in the design file), and an image file with plots of the growth curves for each strain. 

For each experiment a .CSV file with BioLector data and a design file with updated strain and condition info is required as input. Feel free to customize the sections, add more details, or make adjustments based on your specific needs. 

***Growth rate calculation:***
The Growth rate (µmax) is estimated by first log-transforming the growth data, fitting a smoothing spline, and then using linear regression on selected growth phases based on the derivative (rate of change) of the log-transformed data. The slope og the regression line within each growth phase gives the growth rate. Overall, this method is robust for identifying periods of exponential growth and estimating the corresponding rates, even in noisy or complex growth data. For statistical confidence and goodness-of-fit info for the regression, the p-value and R-squared values are provided. Always ensure the calculated growth rates make biological sense by inspecting plotted data. 

***Data:***
This repository was created in conjunction with our recent publication, where the data provided in the "data" folder was analyzed using the included script. The same analysis can be reproduced with this script, as detailed in the publication. For more information see "Characterizing heterologous protein burden in Komagataella phaffii". 
DOI: 10.5281/zenodo.13625781

<img src="https://github.com/user-attachments/assets/cc096bc4-2b67-4081-a8b5-462e737f6857" width="200" height="100">
