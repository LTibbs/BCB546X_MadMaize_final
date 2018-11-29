# PCA in Eigensoft 3.0: Fig. 1-2
Laura Tibbs

## Intro:
All of the following Unix code was run on the hpc-class server. The R code was run using: `R version 3.4.4 (2018-03-15), Platform: x86_64-w64-mingw32/x64 (64-bit), Running under: Windows >= 8 x64 (build 9200)`.

## Step 1: Install Eigensoft 3.0
Eigensoft 3.0 is available through a GitHub repository. I ran the following code on the hpc-class server to download and then make and test the program:

	#clone the repository into the directory that itself contains BCB546X_MadMaize_final on the hpc-class:
	# I did not clone into the BCB546X_MadMaize_final in GitHub because the Eigensoft program is quite large and took up too much space in our directory
	git clone https://github.com/argriffing/eigensoft.git 
	
	# make the program (following instructions in the README for the most recent Eigensoft version, found at https://github.com/DReichLab/EIG)
	cd eigensoft/src
	make
	make install

	# run the example to check that the program works (following instructions in Eigensoft's README):
	cd ../EIGENSTRAT
	perl example.perl

## Step 2: Format data for Eigensoft 3.0: `2_Eigensoft_geno_format.Rmd`
I used R to format the SNP data for the Eigensoft 3.0 program. This code is found in `2_Eigensoft_geno_format.Rmd` and the associated `html` file. See that file for code and workflow for this step.

## Step 3: Run Eigensoft: `3_Eigensoft_PCA.perl`
First, copy the output files from Step 2 to the directory `eigensoft/EIGENSTRAT` in order to use the Eigensoft program on them.

	#in the BCB546X_MadMaize_final/code/PCA directory, run:
	cp full.eigensoft.indiv ../../../eigensoft/EIGENSTRAT/
	cp full.eigensoft.geno ../../../eigensoft/EIGENSTRAT/
	cp SP.geno ../../../eigensoft/EIGENSTRAT/
	cp SP.indiv ../../../eigensoft/EIGENSTRAT/	
	cp SLC.geno ../../../eigensoft/EIGENSTRAT/
	cp SLC.indiv ../../../eigensoft/EIGENSTRAT/
	cp SLL.geno ../../../eigensoft/EIGENSTRAT/
	cp SLL.indiv ../../../eigensoft/EIGENSTRAT/
	cp Fig2d.geno ../../../eigensoft/EIGENSTRAT/
	cp Fig2d.indiv ../../../eigensoft/EIGENSTRAT/
	cp eigensoft.snp ../../../eigensoft/EIGENSTRAT/

Use Eigensoft's CONVERTF functionality to create `.pheno` files from the `.indiv` files:

	#in the eigensoft/CONVERTF directory, run:
	./ind2pheno.perl ../EIGENSTRAT/full.eigensoft.indiv ../EIGENSTRAT/full.eigensoft.pheno
	./ind2pheno.perl ../EIGENSTRAT/SP.indiv ../EIGENSTRAT/SP.pheno
	./ind2pheno.perl ../EIGENSTRAT/SLC.indiv ../EIGENSTRAT/SLC.pheno
	./ind2pheno.perl ../EIGENSTRAT/Fig2d.indiv ../EIGENSTRAT/Fig2d.pheno

This does output the warning `WARNING: no cases` and `WARNING: no controls` but that is fine because this paper did not have accessions labelled as cases or controls; instead, we labelled them according to species or group. According to the documentation, the program will still work without these labels.
	
Then, using the `example.perl` file as a template, I wrote a `.perl` script to run Eigensoft on the full and subset data (that is, `full.eigensoft`, `SP`, `SLC`, `SLL`, and `Fig2d`). This is in the file `3_Eigensoft_PCA.perl`.	Run this code as follows:

	# move the script to the necessary directory; in the BCB546X_MadMaize_final/code/PCA directory, run:
	cp full.eigensoft.indiv ../../../eigensoft/EIGENSTRAT/
	
	# in the eigensoft/EIGENSTRAT directory, run:
	perl 3_Eigensoft_PCA.perl

The error "OOPs bad phenotype ..." is thrown by this process, but all required output is created (the error occurs after the PCA itself), so can move on despite this error.

Now, use Eigensoft's native `ploteig` function to plot these PC graphs. However, from information on a known issue on GitHub (https://github.com/DReichLab/EIG/issues/13), I first had to use vi to comment out the call to `fixgreen` within the `ploteig` file. Then, I could run ploteig as follows:

	
	module load gnuplot/5.2.0-py2-fnj3nkc # load the gnuplot module required by `ploteig`
	
	# in the eigensoft/EIGENSTRAT directory:	
	
	# plot the full dataset: -i gives the input file (`.pca.evec` output by `3_Eigensoft_PCA.perl` from above, -c gives the PC range desired to plot (here, PC 1 and PC 2), -p gives the list of all potential categories into which individuals are grouped in the input file, and -o gives the output file name
	perl ../bin/ploteig -i full.eigensoft.pca.evec -c 1:2 -p SP:SG:SLC:SLL:mixture:SN -x -o Fig1_Eigensoft_ploteig.xtxt 

	
	# plot for Fig 2:
	perl ../bin/ploteig -i SP.pca.evec -c 1:2 -p SP-Montane-1:SP-Peru-2:SP-Peru-9:SP-Ecuador-1:SP-Ecuador-2:SP-Peru-1:SP-Peru-4:SP-Montane-2:SP-Peru-3:SP-Peru-8:SP-Ecuador-3:SP-Peru-5:SP-Peru-6:SP-mixture:SP-non-Andean:SP-Peru- -x -o Fig2-SP-ploteig.xtxt 
	perl ../bin/ploteig -i SLC.pca.evec -c 1:2 -p SLC-Sinaloa:SLC-Mesoamerica:SLC-mixture:SLC-world:SLC-Ecuador-3:SLC-Asia:SLC-Peru-1:SLC-SP-Peru:SLC-Ecuador-1:SLC-Peru-2:SLC-Costa-Rica:SLC-Ecuador-2:SLC-Colombia:SLC-1:SLC-vintage:SLC-Peru-3:SLC-LA2135:SLC-LA2792:SLC-Wva106 -x -o Fig2-SLC-ploteig.xtxt 
	perl ../bin/ploteig -i SLL.pca.evec -c 1:2 -p SLL-Mesoamerica:SLL-vintage-fresh:SLL-vintage-1:SLL-processing-1-2:SLL-processing-1-1:SLL-1:SLL-early-breed:SLL-fresh-1:SLL-vintage-2:SLL-processing-2:SLL-processing-1-3:SLL-fresh-2 -x -o Fig2-SLL-ploteig.xtxt 
	perl ../bin/ploteig -i Fig2d.pca.evec -c 1:2 -p SLL-Mesoamerica:SLC-Mesoamerica:SLC-world:SLC-Ecuador-3:SLC-Peru-1:SLC-Peru-2:SLC-Costa-Rica:SLC-Colombia:SLC-vintage -x -o Fig2-Fig2d-ploteig.xtxt 

Now that output files from the PCA and associated graphs have been made, move them to GitHub.

	# move output files to GitHub; run in eigensoft/EIGENSTRAT directory:
	cp full.eigensoft.pca.evec ../../BCB546X_MadMaize_final/code/PCA/eigensoft/
	cp SP.pca.evec ../../BCB546X_MadMaize_final/code/PCA/eigensoft/
	cp SLC.pca.evec ../../BCB546X_MadMaize_final/code/PCA/eigensoft/
	cp SLL.pca.evec ../../BCB546X_MadMaize_final/code/PCA/eigensoft/
	cp Fig2d.pca.evec ../../BCB546X_MadMaize_final/code/PCA/eigensoft/
	cp Fig1_Eigensoft_ploteig.pdf ../../BCB546X_MadMaize_final/code/PCA/figures/
	cp Fig2-SP-ploteig.pdf ../../BCB546X_MadMaize_final/code/PCA/figures/
	cp Fig2-SLC-ploteig.pdf ../../BCB546X_MadMaize_final/code/PCA/figures/
	cp Fig2-SLL-ploteig.pdf ../../BCB546X_MadMaize_final/code/PCA/figures/
	cp Fig2-Fig2d-ploteig.pdf ../../BCB546X_MadMaize_final/code/PCA/figures/
		
	
## Step 4: Plot PCs in R: `4_PCA_plotting_code.Rmd`

The graphs produced by Eigensoft's default plotting function(`ploteig`) were similar overall to those in the paper, but had some differences. In order to more closely replicate the figures, I used R on the same input files as used to create these graphs in Eigensoft's `ploteig`. The code for this section is found in `4_PCA_plotting_code.Rmd` and the associated `html` file.

## Step 5: Tidy files in repository:

Move figures to a directory, Eigensoft input into the Eigensoft directory, any other intermediate files???

Thoughts--maybe I copy `example.perl`, edit it, and go from there?? Or is it easier just to rename all my files?? I don't see where to rename the example.pheno.

###
trying to plot: in EIGENSTRAT,

	 perl ../bin/ploteig -i example.pca.evec -c 1:2 -p SLC:SLL:mixture:SP:SN:SG:SLL:conflicting -x -o Fig1_Eigensoft_plotting.xtxt

But this causes--

	SLC:SLL:mixture:SP:SN:SG:SLL:conflicting -x -o test_example.xtxt
	## number of fields: 5
	SLC:SLL:mixture:SP:SN:SG:SLL:conflicting
	sh: gnuplot: command not found
	Can't exec "fixgreen": No such file or directory at ../bin/ploteig line 118, <FF> line 1009.

https://github.com/DReichLab/EIG/issues/13 said to comment out the fixgreen, which I did, but still had the gnuplot problem.

Do: `module load gnuplot/5.2.0-py2-fnj3nkc`

