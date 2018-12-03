# Overview of original paper: Blanca, et al. 2015

# Overview of analysis and comparison to original results

## Pre-processing data: Laura Tibbs and Ben Cortes
The code associated with this section is found in `code/Processing_Data`. 
### BEN--do you want to finish this section? what we did and how it compares?

## Principal Component Analysis: Laura Tibbs
The annotated code associated with this section is found in `code/PCA`, and the resulting figures are found in `code/PCA/figures`. The analysis corresponds to Fig. 1-2 in Blanca, et al. 

See `1_Tibbs final assignment writeup.md` for a detailed overview of this section, and the files `2_Eigensoft_geno_format.Rmd` and `4_PCA_plotting_code.Rmd` (as well as the associated `.html` files) for detailed overviews of steps within this section.

In general, this section ran principal component analysis on both the full dataset (after pre-processing, see above) and on subsets of it that represented subpopulations. The steps were: 

1. Eigensoft 3.0
	- The paper said that they used Eigensoft version 3.0 to perform the PCA. 
	- In this step, I downloaded and installed Eigensoft 3.0. As I describe in the markdown file, I downloaded this to `hpc-class` and did not add it to GitHub because, among other reasons, the program files were quite large and caused us to go over our storage limit on GitHub if I tried to add them.

2. Format data for Eigensoft 3.0
	-	Eigensoft 3.0 required input files that were formatted in a custom format. I used R to manually convert from the file type of the data we had to that required for Eigensoft.
	-	This step was accomplished and described in `2_Eigensoft_geno_format.Rmd`

3. Run PCA in Eigensoft:
	- In this step, I ran the PCA using Eigensoft on the `hpc-class` server.
	- First, I had to move the files I prepared to the directory where the Eigensoft program was.
	- Then, I used the `ind2pheno.perl` script from Eigensoft's `CONVERTF` directory to create additional input files for the PCA.
	- I used the `example.perl` file as a template to write `3_Eigensoft_PCA.perl` and ran this to carry out the PCA on the data prepared above.
	- I used the Eigensoft program's native `ploteig` function to plot the results.
	- Finally, I moved the output files and associated graphs to the directory associated with the GitHub repository.

4. Plot Figures 1-2 in R
	- I used R to read in and tidy the output files from the Eigensoft PC analysis.
	- Then, I used these files to produce scatterplots of the PCA to replicate Fig 1-2 in Blanca, et al. At this step, I found that the signs of the PCs were flipped in several of my analyses compared to those in the paper; however, because the interpretation of a PCA remains the same when all its signs are flipped, I flipped the signs  in my own results where necessary to reproduce the published figures. 
	- Finally, I re-plotted the data to recreate the appearance of the plots in the paper (colors, theme, etc.) as nearly as possible.
	- See `4_PCA_plotting_code` for more details.

5. Tidy repository:
	- I moved intermediate files produced by the analysis to tidy up the directory.

When compared with the results in the published paper, Fig 1-2 appeared quite similar. The overall layout of the groups in each of the PC graphs was quite close to those in the paper, though they did not match exactly. However, the percent of variation explained by each of the first two principal components was higher in my analysis than in their published paper in all cases. Any discrepancies are probably due to differences in data used; after data processing, they used 952 accessions and 2313 markers, while our processed data included 950 accessions and 2959 markers despite our best efforts to replicate their processing methods (described above).