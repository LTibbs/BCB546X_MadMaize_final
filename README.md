# BCB546X_MadMaize_final
Dec. 3, 2018

Final project for BCB546X, MadMaize project group: Laura Tibbs, Ben Cortes, Qi Mu, Jialu Wei, and Andy Herr

This was based on the paper Blanca, José, et al. "Genomic variation in tomato, from wild ancestors to contemporary breeding accessions." BMC genomics 16.1 (2015): 257. This paper can be found at https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1444-1 .

## Overview of repository and project:

### Home directory:
Contains files that give an overview of the project (markdown format), the "custom scripts" mentioned by Dr. Blanca in his paper and provided to us via email, and the final presentation. Also contains the following subdirectories:

### Code directory:
Contains all code used in this project.

#### Processing_Data directory: 
Contains code used to pre-process data as described in the Methods of Blanca, et al. Also contains intermediate files from this process. The data produced by this code was stored in the `code` directory and used as input files for the other analyses in this project. Code written and annotated by Laura Tibbs and Ben Cortes.

#### PCA directory: 
Contains code used to process data for and run PCA analysis, then to make figures for this. The output figures are found in the subdirectory `figures`. Also contains intermediate files from this process; these were associated with the Eigensoft v3.0 program, so the files are found in `eigensoft`. This analysis corresponds to Fig 1-2 in Blanca, et al. Code written by Laura Tibbs.

### Classification_Comparison directory:
Contains code used to process classification data so that the Passport and Genetic Classifications could be compared between samples. To do this the data was converted into a pandas data frame, filtered and sorted. Counts were taken on how many times a pare of classifications occurred. This processed data was then used to create a plot graph showing the linear comparison between the two classifications. Written and coded by Andy Herr

#### Neighbor_Network directory: 
Contains code (R markdown file) used to process data for and run Neighbor network analysis. The output files and intermediate files can be found in the sub-directory `Figure5_Output`. Explicit annotations can be found inside the R markdown code file. The R code was used to generate the stats for neighbor network, and the actual figure was produced by interface software SplitsTree4 (v 4.14.6) as indicated by the paper. A similar figure can be generated in R, however, without proper available package, the quality was low so I did not include it in results. This analysis corresponds to Fig 5 (A and B) in Blanca, et al. Code written and annotated by Qi Mu.

#### Phylogenetic_Analysis directory: 
Contains code used to process data for phylogenetic analysis, then to make figures for this. The output figures can be found in the subdirectory `figures`. Software TASSEL was used for clustering with Neighbor Joining method instead of Beast as they did; input and output files of it can be found in `TASSEL`. This analysis corresponds to Fig 6 in Blanca, et al. Code written by Jialu Wei.

#### Rarefaction_Analysis directory: 
Contains code used to process data for rarefaction analysis, then to make figures for this. The output figures can be found in the subdirectory `figures`. ADZE1.0 software was used for analysis as indicated by the paper; input and output files of it can be found in `ADZE`. This analysis corresponds to Fig 7 in Blanca, et al. Code written by Jialu Wei.

#### Weight_Shape directory:
Contains code used to replicate Figure 8. This involved first organizing accessions by species then by group, then plotting the standardized percentage of ancestral and derived alleles in a group for each weight and shape locus, and finally calculating binomial confidence intervals. Code is annotated and writen in the file README.Rmd. Written and annotated by Ben Cortes

### Data directory: 
Contains original data from the tables in Blanca, et al. (`12864_2015_1444_MOESM1_ESM.txt` /`12864_2015_1444_MOESM1_ESM.csv`/ `Suppl_Table_1.csv` and `Suppl_Table_2.csv`) as well as a paper it referenced (Sim, Sung-Chur, et al. "Development of a large SNP genotyping array and generation of high-density genetic maps in tomato." PloS one 7.7 (2012): e40563; found at https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0040563) (these files are `S1_SNPs_Sim_etal.xlsx` and `S8_Map_Sim_etal.csv`. Also contains data processed using a script in the `code/Processing_Data` directory (`final_data.csv` and `rarefaction_LD_data.csv`).
