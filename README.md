# Finding proteins with overrepresentation of dipeptide motifs.

This tool aims to find proteins which exhibit a statistically significant overrepresentation of motifs consisting of two amino acids. The are also scripts to visualize the aminoacid content in found proteins, compare the proteins to their orthologs using multiple-sequence alignment and visualize the output and sort proteins into ranked lists. Currently in this repository the scripts are set to analyse the proteome of Arabidopsis thaliana, but with minor modifications any organisms' proteome can be used.

![Program output](https://user-images.githubusercontent.com/18538056/89332595-9366af80-d68b-11ea-9f04-ec21c58795c9.png)

The main calculations are done by using calculate_zscore.py, which outputs 210 files in json format - one file for each combination of two amino acids. 
For each protein 4 parameters are calculated:

* **count:** the number of dipeptides found in this protein
* **length:** the length of the protein
* **zscore:** the z-score, number of standard deviations that the protein differs from a protein with the mean amount of calculated dipeptides
* **ratio:** ratio of the protein's content of dipeptides to the number of dipeptides that appear on average in all proteins

Full results for Arabidopsis thaliana are available in json format in the full_fixed_results directory.

## Screenshots of results of the analysis of A. thaliana proteome.

Users are able to generate HTML files with a proteins parameters, as well as a multiple-sequence alignment with examined motifs higlighted.

![Output file html](https://i.imgur.com/ZSzm6cJ.png)
![Output file html](https://i.imgur.com/OPU7Uc8.png)

Out of 937 found proteins with a statistically significant overrepresentation of motifs, most exhibit overrepresentation of only one motif, but some exhibit overrepresentation of up to 91 different dipeptide motifs.

![Output file html](https://i.imgur.com/PZHX5lf.png)

A histogram grouping proteins by the count of the WG/GW motif.

![Output file html](https://i.imgur.com/nifuFo7.png)
