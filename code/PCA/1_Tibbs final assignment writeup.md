# PCA in Eigensoft 3.0: Fig. 1-2
Laura Tibbs

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

## Step 2: Format data for Eigensoft 3.0
I used R to format the SNP data for the Eigensoft 3.0 program. This code is found in `2_Eigensoft_geno_format.Rmd`. See that file for code and workflow for this step.

## Step 3: Run Eigensoft
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

This does output the warning `WARNING: no cases` and `WARNING: no controls` but that is fine because this paper did not have accessions labelled as cases or controls. According to the documentation, the program will still work without these labels.
	
Then, using the `example.perl` file as a template, I wrote a `.perl` script to run Eigensoft on the full and subset data (that is, `full.eigensoft`, `SP`, `SLC`, `SLL`, and `Fig2d`). This is in the file `3_Eigensoft_PCA.perl`.	Run this code as follows:

	# move the script to the necessary directory; in the BCB546X_MadMaize_final/code/PCA directory, run:
	cp full.eigensoft.indiv ../../../eigensoft/EIGENSTRAT/
	
	# in the eigensoft/EIGENSTRAT directory, run:
	perl 3_Eigensoft_PCA.perl

The error "OOPs bad phenotype" is thrown by this process, but all required output is created (the error occurs after the PCA itself), so can move on despite this error.

Now, use Eigensoft's native `ploteig` function to plot these PC graphs. However, from information on a known issue on GitHub (https://github.com/DReichLab/EIG/issues/13), I first had to use vi to comment out the call to `fixgreen` within the `ploteig` file. Then, I could run ploteig as follows:

	# in the eigensoft/EIGENSTRAT directory:
	module load gnuplot/5.2.0-py2-fnj3nkc # load the gnuplot module required by `ploteig`
	
	# plot the full dataset: -i gives the input file (`.pca.evec` output by `3_Eigensoft_PCA.perl` from above, -c gives the PC range desired to plot (here, PC 1 and PC 2), -p gives the list of all potential categories into which individuals are grouped in the input file
	perl ../bin/ploteig -i full.eigensoft.pca.evec -c 1:2 -p SLC:SLL:mixture:SP:SN:SG:SLL:conflicting -x -o Fig1_Eigensoft_plotting.xtxt 
	# move output files for the full dataset to GitHub:
	cp full.eigensoft.pca ../../../BCB546X_MadMaize_final/code/PCA
	cp full.eigensoft.pca.evec ../../../BCB546X_MadMaize_final/code/PCA

#### TO DO: make real .snp files for each dataset; re-run with these; make and move all output files/plots to github
	
	
## Step 4: Plot PCs in R:


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

