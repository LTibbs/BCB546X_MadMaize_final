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