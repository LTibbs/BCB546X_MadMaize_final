# PCA in Eigensoft 3.0: Fig. 1-2 and S1-3
Laura Tibbs

## Step 1: Install Eigensoft 3.0
Eigensoft 3.0 is available through a GitHub repository. I ran the following code on the hpc-class server to download and then make the program:

	#clone the repository:
	git clone https://github.com/argriffing/eigensoft.git 
	
	# make the program (following instructions in the README for the most recent Eigensoft version, found at https://github.com/DReichLab/EIG)
	cd eigensoft/src
	make
	make install

	# run the example to check that the program works (following instructions in eigensoft's README):
	cd ../EIGENSTRAT
	perl example.perl

## Step 2: Format data for Eigensoft 3.0
I used R to format the 

## Step 3: Run Eigensoft
First, use Eigensoft's CONVERTF functionality to create a `.pheno` file from the `.indiv` file:

	#in the eigensoft/CONVERTF directory, use this command:
	./ind2pheno.perl ../EIGENSTRAT/full.eigensoft.indiv ../EIGENSTRAT/example.pheno

Then, change all the names in the script (or change the file names themselves?? need to decide...) and run:

	perl example.perl

Thoughts--maybe I copy `example.perl`, edit it, and go from there?? Or is it easier just to rename all my files?? I don't see where to rename the example.pheno.
