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

 
