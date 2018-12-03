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

# don't need this part for this paper, so comment out:
#$command = "gc.perl full.eigensoft.chisq full.eigensoft.chisq.GC";
#print("$command\n");
#system("$command");

# Again, edit file names to run Eigensoft on the SP dataset:

$command = "smartpca.perl";
$command .= " -i SP.geno ";
$command .= " -a eigensoft.snp ";
$command .= " -b SP.indiv " ;
# $command .= " -k 2 "; 
$command .= " -o SP.pca ";
$command .= " -p SP.plot ";
$command .= " -e SP.eval ";
$command .= " -l SP.log ";
#$command .= " -m 5 ";
#$command .= " -t 2 ";
#$command .= " -s 6.0 ";
print("$command\n");
system("$command");

$command = "smarteigenstrat.perl "; 
$command .= " -i SP.geno ";
$command .= " -a eigensoft.snp ";
$command .= " -b SP.indiv ";
$command .= " -p SP.pca ";
#$command .= " -k 1 ";
$command .= " -o SP.chisq ";
$command .= " -l SP.log ";
print("$command\n");
system("$command");

# don't need this part for this paper, so comment out:
#$command = "gc.perl full.eigensoft.chisq full.eigensoft.chisq.GC";
#print("$command\n");
#system("$command");

# Again, edit file names to run Eigensoft on the SLC dataset:

$command = "smartpca.perl";
$command .= " -i SLC.geno ";
$command .= " -a eigensoft.snp ";
$command .= " -b SLC.indiv " ;
# $command .= " -k 2 "; 
$command .= " -o SLC.pca ";
$command .= " -p SLC.plot ";
$command .= " -e SLC.eval ";
$command .= " -l SLC.log ";
#$command .= " -m 5 ";
#$command .= " -t 2 ";
#$command .= " -s 6.0 ";
print("$command\n");
system("$command");

$command = "smarteigenstrat.perl "; 
$command .= " -i SLC.geno ";
$command .= " -a eigensoft.snp ";
$command .= " -b SLC.indiv ";
$command .= " -p SLC.pca ";
#$command .= " -k 1 ";
$command .= " -o SLC.chisq ";
$command .= " -l SLC.log ";
print("$command\n");
system("$command");

# don't need this part for this paper, so comment out:
#$command = "gc.perl full.eigensoft.chisq full.eigensoft.chisq.GC";
#print("$command\n");
#system("$command");

# Again, edit file names to run Eigensoft on the SLL dataset:

$command = "smartpca.perl";
$command .= " -i SLL.geno ";
$command .= " -a eigensoft.snp ";
$command .= " -b SLL.indiv " ;
# $command .= " -k 2 "; 
$command .= " -o SLL.pca ";
$command .= " -p SLL.plot ";
$command .= " -e SLL.eval ";
$command .= " -l SLL.log ";
#$command .= " -m 5 ";
#$command .= " -t 2 ";
#$command .= " -s 6.0 ";
print("$command\n");
system("$command");

$command = "smarteigenstrat.perl "; 
$command .= " -i SLL.geno ";
$command .= " -a eigensoft.snp ";
$command .= " -b SLL.indiv ";
$command .= " -p SLL.pca ";
#$command .= " -k 1 ";
$command .= " -o SLL.chisq ";
$command .= " -l SLL.log ";
print("$command\n");
system("$command");

# don't need this part for this paper, so comment out:
#$command = "gc.perl full.eigensoft.chisq full.eigensoft.chisq.GC";
#print("$command\n");
#system("$command");

# Again, edit file names to run Eigensoft on the Fig2d dataset:

$command = "smartpca.perl";
$command .= " -i Fig2d.geno ";
$command .= " -a eigensoft.snp ";
$command .= " -b Fig2d.indiv " ;
# $command .= " -k 2 "; 
$command .= " -o Fig2d.pca ";
$command .= " -p Fig2d.plot ";
$command .= " -e Fig2d.eval ";
$command .= " -l Fig2d.log ";
#$command .= " -m 5 ";
#$command .= " -t 2 ";
#$command .= " -s 6.0 ";
print("$command\n");
system("$command");

$command = "smarteigenstrat.perl "; 
$command .= " -i Fig2d.geno ";
$command .= " -a eigensoft.snp ";
$command .= " -b Fig2d.indiv ";
$command .= " -p Fig2d.pca ";
#$command .= " -k 1 ";
$command .= " -o Fig2d.chisq ";
$command .= " -l Fig2d.log ";
print("$command\n");
system("$command");

# don't need this part for this paper, so comment out:
#$command = "gc.perl full.eigensoft.chisq full.eigensoft.chisq.GC";
#print("$command\n");
#system("$command");