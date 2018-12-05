# Overview of original paper: Blanca, et al. 2015

 In agronomy, comparative genomics is used to elucidate the many questions of how domestication events shaped the genotypes and phenotypes of crop species. What phenotypes did the ancestor possess, and what phenotypes were selected for during domestication?
 How do the modern lineages differ from their ancestor? When and where did the domestications occur, and what is the overall phylogenetic relationship between modern lines of cultivated crops and their wild-type ancestors? 
 In this paper, the authors sought to answer these questions as they pertained to tomatoes and to provide a description of the evolutionary history and domestication of the domesticated tomato, Solanum lycopersicum, from its closest wild ancestor, Solanum pimpinellifolium (SP).
 There are two botanical varieties of S. lycopersicum, S. l. cerasiforme (SLC) and S. l. lycopersicum (SLL), and it is hypothesized that SLC was domesticated from SP before SLL originated as a further domestication of SLC. 
 The authors sought to confirm this hypothesis, as well as to determine the location of domestication events in tomato and the ancestral characteristics of genes which affect the weight and shape of the tomato fruit. 
 
 In this paper, the authors began by genotyping a wide variety of 1008 tomato accessions (952 of which were unique) using the Tomato Infinium Array. This array probed for 8784 biallelic SNPs, 7720 of which passed quality control and were used for the raw dataset.
 SNPs which would not have been helpful for differentiating the genomes of the 1008 accessions were removed. The removed SNPs included those with greater than 10% missing data and those with a major allelic frequency of 0.95. 
 Additionally, SNPs which mapped too closely together to be distinguishable (0.1 cM or closer) were filtered out, resulting in a final dataset of 2313 SNPs. Duplicate accessions were also corrected, and 8 with inconsistent data were removed as well.
 
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

## Neighbor Network Analysis: Qi Mu  
This part is associate with the reconstruction of Figure 5 in Blanca et al. 2015 paper. The annotated code can be found in `code/Neighbor_Network` folder. The code is written in R markdown file. The R project and knitted R markdown files are also included. Please note that, knitted R markdown file cannot be viewed properly directly from Github, they have to cloned or downloaded first. The output files are stored in `code/Neighbor_Network/Figure5_Output` folder. Two intermediate files are also included in this folder. The R markdown code itself is well annotated, so here I will summarized the meaning of figure 5, reasons I used this method, overall flow, some technical details, and reproducibility of the original paper.  

1. Explanation of Neighbor Network Analysis (Neighbor-Net)  
	Neighbor-Net is first presented in a paper from Bryant and Moulton at *Mol. Biol. Evol.* 2004. It is a distance based method for constructing phylogentic network that is based on the Neighbor-Joining (NJ) algorithm (Bryant and Moulton, 2004). So instead of using sequence data, the Neighor-Net phylogentic network 

2. Genetic distance matrix  
	Genetic distance is a measure of the genetic divergence between species or between populations within a species. In order to construct the Netghbor-Net, we will need the genetic distance matrix computed by among each pair of groups. Every group (composite of varied numbers of individuals) will have a distance between every other group, yield a distance matrix for all the groups. There are several method to calculate genetic distances.  
	- The fixation index (Fst) is a measure of the population differentiation due to genetic structure. A high Fst value between two populations means the allele fixation is high for each population, that they might have inbreeding for a long time and show low levels of breeding with one another. A low Fst value between two populations means that the allele fixation is low, meaning they have high level of inter-breeding.  
	- Gst is another common term to measure genetic differentiation. However, there are some drawbacks using Gst based on `Jost 2008 paper`. So instead, the `Blanca et al. 2015` paper used a improved term called ture Differentiation estimator Dest, which can unbiased estimate the diversity (genetic distance) among population.  
  
3. The script problem  
	The paper has used custom python scripts for processing data. They also wrote a custom python library to compute genetic distance using the algorithm of Dest. The scripts and library were not included in the paper, so I have asked the author to share them with us. The script they sent us is an integrete script that combined everything that has carried out in this paper, without any proper annotations. I have tried to undertand their script for a couple of weeks, and figured out that most of the functions were based on Python 2 (keep running into bugs with Python 3). The script and the library is a complicated system, there were more than 20 individual python scripts, and the variables and packages are all cross-referencing in different scripts. It took about 3-5 steps to go through the scripts and locate the variable, class, or fuction's definitaions. It took my several nights to trace back the variables and functions, and I have stuck at the input file format, which requirs a metadata on the genotype file. I have asked again from the author on the format of input data file, he said he couldn't recall what he did then. The system was outdated and they have swicthed to a new system after this paper. He said he would run into probelms like I did if he had to redo the analysis. Hence, I decide to not using their python code and instead finding similar pacckages in R and generate figures close to theirs as much as possible. 

4.  Overall work flow  
	I had moved on to R after the unsuccessful python decoding. Here is the work flow I carried out:  
	
    1) Data processing (tidyverse). The genotype file (Suppl_Table_1) is combined with grouping information (Suppl_Table_2). Individuals were filtered so that groups with less than 5 individuals were removed, as well as the mixture groups which contains individuals with unspecific grouping information. This step reduced the number of individuals from 1008 to 760 and reduced genetic subgroups from 52 to 36.  
    	
    2) Data formatting. Packages `adegenet` and `hierfstat` were used in this steps. To work with the data, we need to convert the data set to a `genind` object.The result can then be converted to a `genind` object (from package adegent). The `genind` object can then easily be converted into a `hierfstat` (package hierfstat) object.`genind` is used to store individual genotypes, `hierfstat` is used to estimate hierarchical F-statistics for any number of hierarchical levels using the method described in Yang, R.C. (1998) Estimating hierarchical F-statistics. Evolution 52(4):950-956. It also contains functions allowing to test the significance of population differentiation.  
    	
    3) I also calculated and visualized the observed and expected heterozygosity, which is not required in making the figure. This helps me understand the data better. I observed there is significant differences between expected and observed heterozygosity. The observed heterozygosity is significantly lower than the expected, meaning different groups are quite homogeneous among themselves. I also calculated the Fst based on each locus, this is the fixation index on each locus, not population-wised differentiation Fst.  
    	
    4) Pairwise Fst matrix. This generated the pairwise genetic distance (Fst) across all the loci among different genetic subgroups. The estimation method is based on Nei, M. (1987). This is really the main part of this analysis. Since it calculates genetic distance between every pair of groups across all loci, the computation took about 35 min on my computer.  
    	
    5) Convert to nexus file format. In order to use the Fst matrix for neighbor-net building, it has to be certain format. So I converted the Fst genetic distance matrix into nexus format.  
    	
    6) Final visualization. In this step, I used the interface software SplitsTree4 to generate the final figure, as it was also used by the authors after their custom python scripts. The figure can also be viewed in R, but withour proper packages, the plot does not look informative. The Fst genetic distance matrix in nexus format was imported to SplitsTree4, and a final figure was produced. Coloring was done in SplitsTree4.  

5. Reproducibility  
	
    1) Technical aspect. Overall the technical aspects regarding to certain computation method was not reproducible. Especially difficult when it comes to custom scripts. Even if I have the original code, it was almost impossible to reproduce because of minimal annotations, version changed, system alterations. The critial calculation of genetic distance is based on Dest in this paper, using costumed script, which was not able to reproduce. Also, there was no available R package that can calculate the specific Dest they used, so the reproducibility is low in this aspect. However, I understand why they choose to wirte a custom library for Dest, this seems to be the most accurate way of calculating population distance based on literature. However, the lack of related R package makes this difficult so they have to generate custom script.If they could develop this package for other people to use, that will be great.  
    	
    2) General aspect. Although I could not use their specific parameter to generate the neighbor-net, the figures I produces are similar enough to theirs to draw their conclusion. There might be small details differed, but general trends aligned well. The domestication history should be correct. Thus, the overall results were robust.

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

