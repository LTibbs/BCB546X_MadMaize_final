# PCA in Eigensoft 3.0: Fig. 1-2
Laura Tibbs

## Step 1: Install Eigensoft 3.0
Eigensoft 3.0 is available through a GitHub repository. I ran the following code on the hpc-class server to download and then make and test the program:

	#clone the repository (in the directory BCB546X_MadMaize_final/code/PCA):
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

Use Eigensoft's CONVERTF functionality to create `.pheno` files from the `.indiv` files:

	#in the eigensoft/CONVERTF directory, run:
	./ind2pheno.perl ../EIGENSTRAT/full.eigensoft.indiv ../EIGENSTRAT/full.eigensoft.pheno
	./ind2pheno.perl ../EIGENSTRAT/SP.indiv ../EIGENSTRAT/SP.pheno
	./ind2pheno.perl ../EIGENSTRAT/SLC.indiv ../EIGENSTRAT/SLC.pheno
	./ind2pheno.perl ../EIGENSTRAT/Fig2d.indiv ../EIGENSTRAT/Fig2d.pheno

Then, change all the names in the script (or change the file names themselves?? need to decide...) and run:

	perl example.perl

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

