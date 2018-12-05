# Overview of original paper: Blanca, et al. 2015

 In agronomy, comparative genomics is used to elucidate the many questions of how domestication events shaped the genotypes and phenotypes of crop species. What phenotypes did the ancestor possess, and what phenotypes were selected for during domestication?
 How do the modern lineages differ from their ancestor? When and where did the domestications occur, and what is the overall phylogenetic relationship between modern lines of cultivated crops and their wild-type ancestors? 
 In this paper, the authors sought to answer these questions as they pertained to tomatoes and to provide a description of the evolutionary history and domestication of the domesticated tomato, Solanum lycopersicum, from its closest wild ancestor, Solanum pimpinellifolium (SP).
 There are two botanical varieties of S. lycopersicum, S. l. cerasiforme (SLC) and S. l. lycopersicum (SLL), and it is hypothesized that SLC was domesticated from SP before SLL originated as a further domestication of SLC. 
 The authors sought to confirm this hypothesis, as well as to determine the location of domestication events in tomato and the ancestral characteristics of genes which affect the weight and shape of the tomato fruit. 
 
 In this paper, the authors began by genotyping a wide variety of 1008 tomato accessions (952 of which were unique) using the Tomato Infinium Array. This array probed for 8784 biallelic SNPs, 7720 of which passed quality control and were used for the raw dataset.
 SNPs which would not have been helpful for differentiating the genomes of the 1008 accessions were removed. The removed SNPs included those with greater than 10% missing data and those with a major allelic frequency of 0.95. 
 Additionally, SNPs which mapped too closely together to be distinguishable (separated by a distance of less than 0.1 cM) were filtered out, resulting in a final dataset of 2313 SNPs. Duplicate accessions were also corrected, and 8 with inconsistent data were removed as well.
 
 For each of the three species, principal component analyses were conducted using the final dataset to separate the accessions into genetic groups. Further PCAs broke these genetic groups down into subgroups. 
 By comparing the characteristics of these groups and subgroups to the geographical metadata of the accessions, the authors determined that classifying the accessions into groups and subgroups matched the results of classification of the accessions based on native location. 
 The relationships between these groups and subgroups were further elucidated through the construction of neighbor networks and a phylogenetic tree. These figures further demonstrated the relatedness and evolutionary history of the tomato accessions.
 GPS and climate data for the native location of each accession was also analyzed, further corresponding with the aforementioned PCA and phylogenetic tree data. The authors also determined heterozygosity and diversity for the accessions, comparing them between the ancestral and domesticated populations.
 As expected, domestication of the wild-type resulted in the higher heterozygosity and diversity of the wild type decreasing, resulting in significantly more homozygosity in domesticated accessions. Rarefaction and LD analysis were used to confirm the absence of bias in the heterozygosity/diversity analysis.
 Finally, the authors genotyped 6 genes involved in fruit weight and shape using PCR. Comparing the frequency of the derived allele to the ancestral allele in each group demonstrated that, in the most domesticated groups, alleles became fixed for domesticated fruit and weight traits, shifting from the ancestral trait. 
 
 The authors' analysis supported the currently accepted hypothesis that SP is the closest wild ancestor of domestic tomatoes, as evidenced by the phylogenetic tree and neighbor networks, higher genetic diversity, and ancestral fruit weight and shape genotypes. Additionally, the authors conclude that the first tomato domestication occurred in S. America.
 Supporting this conclusion, N. Ecuadorian SP accessions appear to be the closest relatives to the SLC group in terms of phylogenetic data, being especially closely related to the more ancestral SLC accessions found in N. Ecuador. However, the N. Ecuador SLC group has significantly higher diversity than the N. Ecuador SP group. 
 Thus, the authors admit that the transition from SP to SLC during domestication is probably more complicated, possibly resulting from secondary gene flow. More confidently, however, it appears that the Peru SLC group is highly related to the most ancestral SLL group.
 Therefore, the authors conclude that a second domestication event delineating SLC from SLL occurred in Mesoamerica. 

# Overview of analysis and comparison to original results

## Pre-processing data: Laura Tibbs and Ben Cortes
This step was done in R. The code associated with this section 'FinalREADME.Rmd' is found in `code/Processing_Data`. 

Starting with a raw dataset of 7720 SNPs, the authors stated that for all analyses, SNPs where >10% of the accessions had missing data were removed. SNPs where the major allelic frequency was >0.95 were removed as well. 
Additionally, for all analyses except for the rarefaction and LD analyses, SNPs that were separated by a distance of less than 0.1 cM were removed. During the PCA analyses, the authors noticed that some accessions were duplicated or had inconsistent data, and these were rectified as well.

1. Imported and transposed raw dataset (Suppl_Table_2.csv in the Data directory).
	- Raw dataset had SNPs as columns. It is much easier to filter rows of data, so the dataset was transposed.

2. Removed SNPs with >10% missing data.
	- Counted NAs per row and removed the row if >10% of the elements were NAs.
	- Authors state that 240 SNPs were removed in this way (7480 remaining). We removed 251 (7469 remaining).

3. Removed SNPs with >0.95 allelic frequency.
	- Used unite to concatenate all the elements in a row into a single string, then removed NAs with str_replace_all.
	- Counted number of occurrences of each base. Divided largest count by the total length of the concatenated string. If result was >0.95, SNP was removed. 
	- Authors state that 1137 SNPs were removed in this way (6343 remaining). We removed 998 (6471 remaining). 

4. Created consensus map.
	- According to the authors, "genetic distances were based on the genetic maps of Sim et al." (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0040563)." However, that reference has eight different genetic maps, and the authors don't state which one(s) they used.
		1. The genetic map (Table S8) seemed the most complete, so it was used as the basis of our code. Table 8 contains the pre-existing genetic maps of two panels of tomato crosses.
		2. The authors would have had to create a consensus map from this genetic map, but they don't state how they did. We used the package LPmerge (Jeffrey Endelman, "Merging Linkage Maps by Linear Programming, 2018). 
	- Sorted through Table S8 and formatted it into 2 tables, one for the genetic map of each tomato cross.
		1. Separated the 2 tables out by chromosome, creating a list of chromosomes. Each chromosome in turn was a list of SNPs and genetic map information for that chromosome.
	- Ran LPmerge with the parameters that it create four consensus maps for each chromosome.
		1. Manually inspected and chose consensus maps based on the lowest mean RMSE, and broke ties using the lowest sd for RMSE (https://potatobreeding.cals.wisc.edu/wp-content/uploads/sites/161/2014/01/LPmerge_tutorial.pdf).
		2. Used rbind to merge the consensus maps for each chromosome into one consensus map for the entire genome with 5296 SNPs.

5. Rectified discrepancies between SNPs listed in raw dataset and those in the consensus map.
	- 4206 SNPs were shared between the raw dataset and the consensus map.
		1. Additionally, the consensus map had 1090 unique SNPs, and the raw dataset had 2265 unique SNPS.
	- The 1090 unique SNPs in the consensus map were actually SNPs in the raw dataset, but they were given different names by by Sim et al.
		1. Names and naming scheme in consensus map SNPs were changed to match raw dataset.
	- However, because we were using the most complete tomato genetic maps, the 2265 unique SNPs in the raw dataset were determined to simply be SNPs that had never been mapped.
		1. Because the authors stated that "SNPs that mapped closer than 0.1 cM were removed", it was decided that the unmapped SNPs would be retained as they were never mapped and so could not be known to be separated by less than 0.1 cM.

6. Removed SNPs separated by less than 0.1 cM
	- The authors did not state how they decided which SNPs to remove.
		1. Therefore, we first calculated the distances between each SNP.
		2. Kept the first SNP found and marked all SNPs within the next 0.1 cM for removal.
	- The algorithm we used does not account for many consecutive SNPs separated by less than 0.1 cM but that span a total distance of 0.1 cM or greater.
		1. To account for this, we searched for SNPs separated by less than 0.1 cM but more than 0 cM.
		2. Only 17 of these SNPs were found, and after manual inspection, 5 SNPs where the above problematic scenario occurred were removed. 
	- Finally, filtered out all SNPs marked for removal. 
	- Authors state that 4030 SNPs were removed (processed dataset of 2313 SNPs). We removed 3510 (processed dataset of 2961 SNPs).

7. Transposed the processed dataset.
	- Had to remove duplicate accessions, so accessions were transposed to rows.
		
8. Removed duplicate accessions (lines that were genotyped twice) using passport information to decide which of each duplicate to remove.
	- Imported passport information (12864_2015_1444_MOESM1_ESM.txt in the Data directory).
	- Changed accession names in the passport information to match accession names in the raw dataset.
	- According to the authors, an accession which was duplicated was removed completely "unless it was clear based on the passport information, which genotype was correct."
	- What exactly this means could vary from the provided example depending on the situation, so was not amenable to coding. We decided that the best thing to do would be to go through and inspect the duplicate data manually after arranging it. The following workflow was used in Excel.
		1. The file used (`duplicate_accessions.csv`) was generated in R (below) by filtering for duplicated accessions in the passport.data tibble.
		2. The accession duplicates were checked for inconsistencies in the columns that contain information from the PCA; that is, the species, group1, and group2 columns.
		3. If there were no inconsistencies found, one sample was chosen at random (using a random number generator) for that accession to keep. The other samples for the accession were marked to delete.
		4. If inconsistencies were found, the passport and geographic data columns were assessed to determine which accession was "correct." We took "conflicting classification" to indicate that a given sample was not the correct sample for that accession. However, if it was not clear which sample was correct, we removed the accession entirely.
		5. Also removed the entire accession if the passport information between the duplicates was inconsistent.
		6. We saved the file with the lines to delete and keep marked as `duplicate_accessions_new.csv`.
	- Authors state they removed 8 duplicate accessions with inconsistent information and had 952 unique accessions in the final dataset. We also removed 8 duplicate accessions with inconsistent information, but our final number of unique accessions was 950.
	- Also removed duplicate accessions for the rarefaction and LD analysis data.
	
Final dataset
	- Ours - 950 accessions, 2961 SNPs
	- Authors - 952 accessions, 2313 SNPs

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

## Neighbor Network Analysis: Qi Mu  
This part is associate with the reconstruction of Figure 5 in Blanca et al. 2015 paper. The annotated code can be found in `code/Neighbor_Network` folder. The code is written in R markdown file. The R project and knitted R markdown files are also included. Please note that, knitted R markdown file cannot be viewed properly directly from Github, they have to cloned or downloaded first. The output files are stored in `code/Neighbor_Network/Figure5_Output` folder. Two intermediate files (trimmed individuals and groups `Geno_group_filtered.csv` and pairwise genetic distances `FstMatrix.nxs` are also included in this folder. The R markdown code itself is well annotated, so here I will summarized the meaning of figure 5, reasons I used this method, overall flow, some technical details, and reproducibility of the original paper.   

1. Explanation of Neighbor Network Analysis (Neighbor-Net)  
	Neighbor-Net is first presented in a paper from Bryant and Moulton at *Mol. Biol. Evol.* 2004. It is a distance based method for constructing phylogentic network that is based on the Neighbor-Joining (NJ) algorithm (Bryant and Moulton, 2004). So instead of using sequence data, the Neighor-Net phylogentic network uses genetics distance data. 

2. Genetic distance matrix  
	Genetic distance is a measure of the genetic divergence between species or between populations within a species. In order to construct the Netghbor-Net, we will need the genetic distance matrix computed among each pair of groups. Every group (composite of varied numbers of individuals) will have a distance with every other group, yielding a distance matrix for all pairwised the groups. There are several method to calculate genetic distances. I will introdce two methods that are related to this study.  
	- The fixation index (Fst) is a measure of the population differentiation due to genetic structure. A high Fst value between two populations means the allele fixation is high for each population, that they might have inbreeding for a long time and show low levels of breeding with one another. A low Fst value between two populations means that the allele fixation is low, meaning they have high level of inter-breeding. Fst is the population distance I used in my analysis.  
	- Gst is another common term to measure genetic differentiation. However, there are some drawbacks using Gst based on `Jost 2008 paper`. So instead, the `Blanca et al. 2015` paper used an improved term called true Differentiation estimator (Dest), which can estimate the diversity (genetic distance) among population unbiasedly.  

  
3. The script problem  
	The paper has used custom python scripts for processing and analyzing data. They also wrote a custom python library to compute genetic distance using the algorithm of Dest. The scripts and library were not included in the paper, so I have asked the authors to share them with us. The script they sent us is an integrete script that combined everything that they carried out in this paper, and without much annotations. I have tried to undertand their script for a couple of weeks, and figured out that most of the functions were based on Python 2 (keep running into bugs with Python 3). The script and the library is a complicated system, there were more than 20 individual python scripts, and the variables and packages are all cross-referenced in different scripts. It took about 3-5 steps locate the variables by tracing back several scripts to find the variable, class, or fuction's definitaions. It took me several nights to figure out part of the the variables and functions, and I have stuck at the input file format, which requirs a metadata on the genotype file. I have asked again from the author on the format of input data file, he said he couldn't recall what he did back then. The system was outdated and they have swicthed to a new system after this paper. He said he would run into probelms like I did if he had to redo the analysis. Hence, I decide to not use their python code and instead, finding similar pacckages in R and generate figures close to theirs as much as possible. 

4.  Overall work flow  
	I had moved on to R after the unsuccessful python decoding. Here is the work flow I carried out:  
	
    1) Data processing (tidyverse). The genotype file (`Suppl_Table_1`) is combined with grouping information (`Suppl_Table_2`). Individuals were filtered so that groups with less than 5 individuals were removed, as well as the mixture groups which contained individuals with unspecific grouping information. This step reduced the number of individuals from 950 to 760 and reduced genetic subgroups from 52 to 36.  
    	
    2) Data formatting. Packages `adegenet` and `hierfstat` were used in these steps. To work with the data, we need to convert the data set to a `genind` object using package `adegent`. The `genind` object can then easily be converted into a `hierfstat` (package hierfstat) object. The `genind` is used to store individual genotypes, `hierfstat` is used to estimate hierarchical F-statistics for any number of hierarchical levels using the method described in `Yang, R.C. (1998) Estimating hierarchical F-statistics. Evolution 52(4):950-956`. It also contains functions allowing to test the significance of population differentiation.  
    	
    3) I also calculated and visualized the observed and expected heterozygosity, which was not required to make the figure. But it helps me understand the data better. I observed there was significant differences between expected and observed heterozygosities. The observed heterozygosity is significantly lower than the expected, meaning different groups are quite homogeneous among themselves. I also calculated the Fst based on each locus, this is the fixation index on each locus, not population-wised differentiation Fst.  
    	
    4) Pairwise Fst matrix. This generated the pairwise genetic distance (Fst) across all the loci among different genetic subgroups. The estimation method is based on Nei, M. (1987). This is really the main part of this analysis. Since it calculates genetic distance between every pair of groups across all loci, requires large amount of computations, which took about **35 min** on my computer.  
    	
    5) Convert to nexus file format. In order to use the Fst matrix for neighbor-net building, it has to be in certain format. So I converted the Fst genetic distance matrix into nexus format.  
    	
    6) Final visualization. In this step, I used the interface software SplitsTree4 to generate the final figure, as it was also used by the authors after their custom python scripts. The figure can also be viewed in R, but withour proper packages, the plot does not look informative. So I kept the final results from SplitsTree4. The Fst genetic distance matrix in nexus format was imported to SplitsTree4, and a final figure was produced. Coloring was done in SplitsTree4.  

5. Reproducibility  
	
    1) Technical aspect. Overall the technical aspects regarding to certain computational method was not reproducible. Especially difficult when it comes to the custom scripts. Even though I have the original code, it was almost impossible to reproduce because of minimal annotations, changed version, system alterations. The critial calculation of genetic distance is based on Dest in this paper, using costumed script (library). This library is overall complicated and poorly annotated, it was not developed for publich use. So is was almost impossible to implement by others without spending large amount of time on it. Also, there was no available R package that can calculate the specific Dest they used, so the reproducibility is low in this aspect. However, I understand why they chose to wirte a custom library for Dest, this seems to be the most accurate way of calculating population distance based on literature. The lack of related R packages made this difficult so they had to generate custom scripts.It would be great if they could develop this package for public use.
    	
    2) General aspect. Although I could not use their specific parameters to generate the neighbor-net, the figures I produces are similar enough to theirs to draw their conclusions. There might be small parts that differed, but general trends aligned well. The domestication history should be correct. Thus, the overall results were robust.

## Phylogenetic Analysis: Jialu Wei
The annotated code associated with this section is found in `code/Phylogenetic_Analysis`, and the resulting figures are found in `code/Phylogenetic_Analysis/figures`. The analysis corresponds to Fig.6 in Blanca, et al. 

Based on the published paper, they conducted phylogenetic analysis with SNAPP package in Beast software, which was build up based on a recently developed method with finite-sites model likelihood algorithm. However, this software is sort of complicated to play with which require specific java environment; I could not be able to figure it out. Then I looked for similar packages in R but it turned out that there was no package using the same algorithm with SNAPP. To make things easier, I decided to use the traditional Neighbor Joining Clustering method for phylogenetic analysis. The steps are:

1. Beast2 and SNAPP
	- I downloaded Beast2 in UNIX environment and loaded SNAPP package inside it.
	- Read manual and figured out the input format.
		- I found that the examples provided by BEAST were all single character for each site coded with ATCG, while our data were diallelic. To format our data correctly, I looked into the python custum script provided by the authors and recoded our data with same rules.
		- SNAPP has high computational demands, so they selected one accession per genetic subgroup in analysis.
		- I manually selected the accessions based on their figure from total SNP datasat in R.
		- Also, the software requirs input file to be .nex or .xml format. R package was used to convert that as well.
	- However, after struggling for a long time with the preparing my data and figuring out the SNAPP instruction. I realized I couldn't run the .xml file with error popping out. (even with their examples)
	- At last, I had to decide leaving this software and looked for other way instead with limited time left.

2. Run phylogenetic analysis in TASSEL
	- I firstly tried to look for package with same algorithm in R but had no result.
	- Considering the time spending in the last step, TASSEL with friendly interface was a easy way to do the phylogenetic analysis.
	- Some data preparation had already done from first step. Here R was used again to have .fasta file (`figure6_sample.fasta`) for TASSEL.
	- Neighbor Joining method was used in TASSEL for clustering.
	- The output of analysis was saved as `figure6_tree.nwk`, which included the clustering information like nodes and branch.length.

4. Plot Figures 6 in R
	- R ggtree pachage was used to read in and tidy the output files from the TASSEL.
	- `figure6_tree.nwk` was read in as basic tree. More adjustments were done to be alike with their results as much as possible.
		- branch placemnet was adjusted to be the same as theirs with flip() function.
		- each branch was labelled with species name followed by sample name.
		- clustered groups were highlighted with different colors, and group names were placed aside.
	- Save figures as `Figure6.png` and `Figure6.pdf`. 
5. Tidy repository:
	- Intermediate files were moved to `./TASSEL` folder and final figures were moved to `./figures` .
	- All codes used above were saved in `Figure6.Rmd` and `Figrue6.html`, which can be found under `Phylogenetic_Analysis` repository.

Comparing our result with theirs:  
1. Individuals in group SLC and SP were clustered together respectively. (well replicated)
2. Peruvian SP showed more basal status than others in our result; in their result, Peruvian SP was the basal group of the red-fruited species. (kind of different but not conflicted)
3. Ecuadorian SP was phylogenetically closest to SLC. (Well replicated)
4. SLC Ecuador1 was basal to the entire SLC group. (Well replicated)
5. SG clustered very close to Ecuadorian SP in their result, while in ours it was basal of Ecuadorian SLC and SP groups. (different but still not conflicted with PCA result)

Considering that we used different algorithms for analysis, the above results were acceptable. So, this section is reproducible as well.

## Rarefaction Analysis: Jialu Wei
The annotated code associated with this section is found in `code/Rarefaction_Analysis`, and the resulting figures are found in `code/Rarefaction_Analysis/figures`. The analysis corresponds to Fig. 7 in Blanca, et al. 

This section ran rarefaction analysis with two sets of markers. The first set included one marker every 0.1 cM (2959 SNPs, SNP data can be found in `data/final_data.csv`) and the second set included 6471 SNPs (SNP data can be found in `data/rarefaction_ld_data.csv`) after removing monomorphic SNPs and with 10% missing. ADZE1.0 was used for analysis and R was used for visualizing. The steps were:

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


Comparing our results with those in the published paper, they look generally similar.  

1. SLC Andean has the highest number of alleles per loci and SP also has high one.
2. SLL processing and SLL vintage have low number of alleles per loci.
3. Speaking of frequecy of private alleles, SP group shows obvious much higher result while the other groups all have similar frequencies which are lower than 0.05. This makes sense considering that SP group is the most ancient group with has high diversity. 
  
Howerver, there are some differences which are also obvious. 
1. The most difference comes from SLC_non_Andean group, which has similar result in figureA,C with SLL vintage group, but in our result, it turns out to have much higher number of alleles per loci. The reason for this could be mainly due to the grouping issue in my opinion, which I mentioned above that I grouped the individuals manually. Their might be some heterozygous individuals misgroupped to SLC_non_Andean group.
2. The x and y axis scale are a little different between the results. Our individual number are less than theirs while the number of alleles per loci has higher level. I considered two reasons for this. First of all, of course, our filtered SNP data are not exactly the same as theirs, especially when filtered by position. Secondly, there was one critical step in rarefaction analysis when we adjusted the parameter of missing data rate. This step affected the final individual numbers a lot when I ran the program. They didn't mention in detail about the parameter setting in the paper, so I could not guarantee that I used the same threshold as they did. Plus, in terms of the number of individuals, the outputs of ADZE use individual as line, which is half of the number of samples; in my figures, I use individuals to mean samples, which makes more sense to me. In this way, if their plot used directly the number from outputs (without dividing by 2), we can have some differences in scale. 

To sum up, the result of this section are reproducible.

>>>>>>> 6a41fad7d5d376c55b5a564983022f08c0464108
