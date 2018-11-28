#!/usr/bin/perl

$ENV{'PATH'} = "../bin:$ENV{'PATH'}"; 
# MUST put smartpca bin directory in path for smartpca.perl to work

# This is the code from example.perl, to use as a template:

#$command = "smartpca.perl";
#$command .= " -i example.geno ";
#$command .= " -a example.snp ";
#$command .= " -b example.ind " ;
#$command .= " -k 2 ";
#$command .= " -o example.pca ";
#$command .= " -p example.plot ";
#$command .= " -e example.eval ";
#$command .= " -l example.log ";
#$command .= " -m 5 ";
#$command .= " -t 2 ";
#$command .= " -s 6.0 ";
#print("$command\n");
#system("$command");

#$command = "smarteigenstrat.perl "; 
#$command .= " -i example.geno ";
#$command .= " -a example.snp ";
#$command .= " -b example.ind ";
#$command .= " -p example.pca ";
#$command .= " -k 1 ";
#$command .= " -o example.chisq ";
#$command .= " -l example.log ";
#print("$command\n");
#system("$command");

#$command = "gc.perl example.chisq example.chisq.GC";
#print("$command\n");
#system("$command");

# Now, update the file names in the above template code to run Eigensoft on the full dataset:
# commented commands are commented out in order to use the default settings when the paper did not
# provide any specific settings to use.

$command = "smartpca.perl";
$command .= " -i full.eigensoft.geno ";
$command .= " -a full.eigensoft.snp ";
$command .= " -b full.eigensoft.indiv " ;
# $command .= " -k 2 "; 
$command .= " -o full.eigensoft.pca ";
$command .= " -p full.eigensoft.plot ";
$command .= " -e full.eigensoft.eval ";
$command .= " -l full.eigensoft.log ";
#$command .= " -m 5 ";
#$command .= " -t 2 ";
#$command .= " -s 6.0 ";
print("$command\n");
system("$command");

$command = "smarteigenstrat.perl "; 
$command .= " -i full.eigensoft.geno ";
$command .= " -a full.eigensoft.snp ";
$command .= " -b full.eigensoft.indiv ";
$command .= " -p full.eigensoft.pca ";
#$command .= " -k 1 ";
$command .= " -o full.eigensoft.chisq ";
$command .= " -l full.eigensoft.log ";
print("$command\n");
system("$command");

$command = "gc.perl full.eigensoft.chisq full.eigensoft.chisq.GC";
print("$command\n");
system("$command");