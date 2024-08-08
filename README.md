# BioLector data analysis

This repository contains a script for analyzing output data from BioLector microbioreactors and a supporting package "htmf" that provides functionalities to facilitate this analysis. The BioLector is a high-throughput microbioreactor system in microtiter plate format used for monitoring e.g. microbial growth (biomass), pH, dissolved oxygen, and fluorescence. The tool provided in this repository allows for the efficient analysis and visualization of data generated from BioLector experiments.

Repository Contents

- BL2_plot_and_mu_script.R: Is the R script for analyzing BioLector data.
- htmf_0.5.1.tar.gz: zip file with the supporting package providing functions necessary for data analysis.
- data/: Folder containing example data files used for demonstration and testing purposes.
- design_file.txt: design file is a table with information on each well required to run the script
- README.md: This file, providing an overview of the repository contents and instructions for usage.

Usage

The output file from BioLector runs is CSV format file containing all data. The R script will process the data in the CSV file and generate an output file "Info" with the results from the run, an output file "Growthrates" with the analysed maximum growth rate (µmax) for each strain (as listed in the design file), and an image file with plots of the growth curves for each strain. 

For each experiment a .CSV file with BioLector data and a design file with updated strain and condition info is required. Feel free to customize the sections, add more details, or make adjustments based on your specific needs. 

<img src="https://github.com/user-attachments/assets/cc096bc4-2b67-4081-a8b5-462e737f6857" width="200" height="100">
