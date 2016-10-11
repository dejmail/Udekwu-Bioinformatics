# Udekwu-Bioinformatics
Bioinformatics scripts for the Udekwu lab - Stockholm University

***find_primers_demultiplex_and_cut.py*** - Using the metadata specification of the QIIME project, this script locates the forward and reverse primers in a FastQ file, and assigns the read to the relevant file (based on the sample name). It outputs a list of files, including files where only Forward, only Reverse or No primers have been found.