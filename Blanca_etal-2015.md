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

## Rarefaction Analysis: Jialu Wei
The annotated code associated with this section is found in `code/Rarefaction_Analysis`, and the resulting figures are found in `code/Rarefaction_Analysis/figures`. The analysis corresponds to Fig. 7 in Blanca, et al. 

In general, this section ran rarefaction analysis with two sets of markers. The first set included one marker every 0.1 cM (2959 SNPs, SNP data can be found in `data/final_data.csv`) and the second set included 6471 SNPs (SNP data can be found in `data/rarefaction_ld_data.csv`) after removing monomorphic SNPs and with 10% missing. ADZE1.0 was used for analysis and R was used for visualizing. The steps were:

1. ADZE1.0
	- The paper said that they used ADZE1.0 to conduct rarefaction analysis. 
	- I downloaded ADZE1.0 to my personal computer in UNIX environment, since it's more convenient for future analysis with command lines to adjust parameters.  
	- Read manual and run examples insides software to know how it works, the format of input files and how to set parameters depending on our data.

2. Format data for ADZE1.0
	- ADZE1.0 required input files were formatted in a custom format. I used R to manually convert from the file type of the data we had to that required for ADZE1.0.
	- Each sample was labeled with its genetic subgroup name.
		- Six genetic subgroup were defined by authors for analysis, but they didn't mention in detail what each subgroup included in terms of the information from passport files. So, I manually inferred and grouped them based on the information in passport doc (column `group2` basically).
		- Combined SNP file with `group2` column from passport doc with corresponding sample names. 
		- Added new column of genetic subgroup as I earlier defined manually.
		- Formatted each sample with 2 lines because they are diploid. Each line indicated one allele and they were all numeric instead of ATCG.  

3. Run rarefaction analysis in ADZE1.0:
	- Place formatted input files to the same folder with ADZE1.0. 
	- Set parameters based on our own data. 
		-	Basically, they need the number of loci, number of individuals, max number of within-group individuals (MAX_G). These can be obtained by inspecing data in R.
		-	One critical parameter here is the threshold for missing data. I started from 1 and examined the output files to see if MAX_G was reached after calculation.
		-	 Manually adjusted parameter and rerun by command line `./ADZE1.0 parameter.txt -t 0.1`, in which `-t` are used for missing data rate setting
		-	 Repeating the last step until MAX_G was reached in our final output. 
	- Two sets of markers were analysis separately in same way and differing in parameters.
	- Output files were generated and I placed four of them to our `./ADZE` repository, which would be used later for visualizing.

4. Plot Figures 7 in R
	- R was used to read in and tidy the output files from the ADZE1.0 rarefaction analysis.
	- ggplot2 package was used for data visualizing.
	- Colors, theme and figure placement were adjusted based on their original ones as much as possible.
	- Saving figures as `Figure7.png` and `Figure7.pdf`.

5. Tidy repository:
	- Intermediate files were moved to `./ADZE` folder and final figures were moved to `./figures` .
	- All codes used above were saved in `Code_for_Figure7.Rmd` and `Code_for_Figrue7.html`, which can be found under `Rarefaction_Analysis` repository.


Compared our results with those in the published paper, they look generally similar.  

1. SLC Andean has the highest number of alleles per loci and SP also has high one.
2. SLL processing and SLL vintage have low number of alleles per loci.
3. Speaking of frequecy of private alleles, SP group shows obvious much higher result while the other groups all have similar frequencies which are lower than 0.05. This makes sense considering that SP group is the most ancient group with has high diversity. 
  
Howerver, there are some differences which are also obvious. 
1. The most difference comes from SLC_non_Andean group, which has similar result in figureA,C with SLL vintage group, but in our result, it turns out to have much higher number of alleles per loci. The reason for this could be mainly due to the grouping issue in my opinion, which I mentioned above that I grouped the individuals manually. Their might be some heterozygous individuals misgroupped to SLC_non_Andean group.
2. The x and y axis scale are a little different between the results. Our individual number are less than theirs while the number of alleles per loci has higher level. I considered two reasons for this. First of all, of course, our filtered SNP data are not exactly the same as theirs, especially when filtered by position. Secondly, there was one critical step in rarefaction analysis when we adjusted the parameter of missing data rate. This step affected the final individual numbers a lot when I ran the program. They didn't mention in detail about the parameter setting in the paper, so I could not guarantee that I used the same threshold as they did. Plus, in terms of the number of individuals, the outputs of ADZE use individual as line, which is half of the number of samples; in my figures, I use individuals to mean samples, which makes more sense to me. In this way, if their plot used directly the number from outputs (without dividing by 2), we can have some differences in scale. 

To sum up, the result of this section are reproducible.
