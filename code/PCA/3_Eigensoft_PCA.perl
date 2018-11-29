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
$command .= " -a eigensoft.snp ";
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
$command .= " -a eigensoft.snp ";
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

# Again, edit file names to run Eigensoft on the SP dataset:

$command = "smartpca.perl";
$command .= " -i SP.eigensoft.geno ";
$command .= " -a eigensoft.snp ";
$command .= " -b SP.eigensoft.indiv " ;
# $command .= " -k 2 "; 
$command .= " -o SP.eigensoft.pca ";
$command .= " -p SP.eigensoft.plot ";
$command .= " -e SP.eigensoft.eval ";
$command .= " -l SP.eigensoft.log ";
#$command .= " -m 5 ";
#$command .= " -t 2 ";
#$command .= " -s 6.0 ";
print("$command\n");
system("$command");

$command = "smarteigenstrat.perl "; 
$command .= " -i SP.eigensoft.geno ";
$command .= " -a eigensoft.snp ";
$command .= " -b SP.eigensoft.indiv ";
$command .= " -p SP.eigensoft.pca ";
#$command .= " -k 1 ";
$command .= " -o SP.eigensoft.chisq ";
$command .= " -l SP.eigensoft.log ";
print("$command\n");
system("$command");

$command = "gc.perl SP.eigensoft.chisq SP.eigensoft.chisq.GC";
print("$command\n");
system("$command");

# Again, edit file names to run Eigensoft on the SLC dataset:

$command = "smartpca.perl";
$command .= " -i SLC.eigensoft.geno ";
$command .= " -a eigensoft.snp ";
$command .= " -b SLC.eigensoft.indiv " ;
# $command .= " -k 2 "; 
$command .= " -o SLC.eigensoft.pca ";
$command .= " -p SLC.eigensoft.plot ";
$command .= " -e SLC.eigensoft.eval ";
$command .= " -l SLC.eigensoft.log ";
#$command .= " -m 5 ";
#$command .= " -t 2 ";
#$command .= " -s 6.0 ";
print("$command\n");
system("$command");

$command = "smarteigenstrat.perl "; 
$command .= " -i SLC.eigensoft.geno ";
$command .= " -a eigensoft.snp ";
$command .= " -b SLC.eigensoft.indiv ";
$command .= " -p SLC.eigensoft.pca ";
#$command .= " -k 1 ";
$command .= " -o SLC.eigensoft.chisq ";
$command .= " -l SLC.eigensoft.log ";
print("$command\n");
system("$command");

$command = "gc.perl SLC.eigensoft.chisq SLC.eigensoft.chisq.GC";
print("$command\n");
system("$command");

# Again, edit file names to run Eigensoft on the SLL dataset:

$command = "smartpca.perl";
$command .= " -i SLL.eigensoft.geno ";
$command .= " -a eigensoft.snp ";
$command .= " -b SLL.eigensoft.indiv " ;
# $command .= " -k 2 "; 
$command .= " -o SLL.eigensoft.pca ";
$command .= " -p SLL.eigensoft.plot ";
$command .= " -e SLL.eigensoft.eval ";
$command .= " -l SLL.eigensoft.log ";
#$command .= " -m 5 ";
#$command .= " -t 2 ";
#$command .= " -s 6.0 ";
print("$command\n");
system("$command");

$command = "smarteigenstrat.perl "; 
$command .= " -i SLL.eigensoft.geno ";
$command .= " -a eigensoft.snp ";
$command .= " -b SLL.eigensoft.indiv ";
$command .= " -p SLL.eigensoft.pca ";
#$command .= " -k 1 ";
$command .= " -o SLL.eigensoft.chisq ";
$command .= " -l SLL.eigensoft.log ";
print("$command\n");
system("$command");

$command = "gc.perl SLL.eigensoft.chisq SLL.eigensoft.chisq.GC";
print("$command\n");
system("$command");

# Again, edit file names to run Eigensoft on the Fig2d dataset:

$command = "smartpca.perl";
$command .= " -i Fig2d.eigensoft.geno ";
$command .= " -a eigensoft.snp ";
$command .= " -b Fig2d.eigensoft.indiv " ;
# $command .= " -k 2 "; 
$command .= " -o Fig2d.eigensoft.pca ";
$command .= " -p Fig2d.eigensoft.plot ";
$command .= " -e Fig2d.eigensoft.eval ";
$command .= " -l Fig2d.eigensoft.log ";
#$command .= " -m 5 ";
#$command .= " -t 2 ";
#$command .= " -s 6.0 ";
print("$command\n");
system("$command");

$command = "smarteigenstrat.perl "; 
$command .= " -i Fig2d.eigensoft.geno ";
$command .= " -a eigensoft.snp ";
$command .= " -b Fig2d.eigensoft.indiv ";
$command .= " -p Fig2d.eigensoft.pca ";
#$command .= " -k 1 ";
$command .= " -o Fig2d.eigensoft.chisq ";
$command .= " -l Fig2d.eigensoft.log ";
print("$command\n");
system("$command");

$command = "gc.perl Fig2d.eigensoft.chisq Fig2d.eigensoft.chisq.GC";
print("$command\n");
system("$command");