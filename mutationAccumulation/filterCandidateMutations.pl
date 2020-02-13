use strict;
use warnings;
use Getopt::Long;

## This perl script takes a list of candidate de novo mutations identified between 
## two samples, the bam files of those two samples, and the reference genome;
## and will split the de novo mutations into 3 groups: "good", "bad" and "triallelic". 

## The bam files should be filtered e.g. for base and mapping qualities, strand bias, etc...
## and samtools should be available somewhere on the $PATH (samtools is used to generate
## mpileup information for each candidate de novo mutation)

## Triallelic sites: position has >2 alleles across the 2 samples
## Good mutations: position has 2 different alleles across the 2 samples + 1 sample has only 1 allele
## Bad mutations: position has 2 different alleles across the 2 sampes + both samples have the 2 alleles

# If you use this script please cite Fazalova & Nevado "Low spontaneous mutation rate and 
# Pleistocene radiation of pea aphids", submitted, Molecular Biology and Evolution

# Author: Bruno Nevado

my $usage = "example usage: filterCandidateMutations.pl -reference Scaffold1.fa -bam1 perMA.bam -bam2 posMA.bam -candidates candidateMutations.txt\n";
my $noisy = 0; 
my $reference = "";
my $candidates = "";
my $bam1 = "";
my $bam2 = "";

GetOptions (  "reference=s" => \$reference,
              "candidates=s"   => \$candidates,
              "bam1=s"  => \$bam1,
              "bam2=s" => \$bam2 )
  or die("Error in command line arguments\n$usage");

my @BAMS = ($bam1, $bam2);


my %iupac;
$iupac{'A'} = ["A","A"];
$iupac{'a'} = ["a","a"];
$iupac{'C'} = ["C","C"];
$iupac{'c'} = ["c","c"];
$iupac{'G'} = ["G","G"];
$iupac{'g'} = ["g","g"];
$iupac{'T'} = ["T","T"];
$iupac{'t'} = ["t","t"];
$iupac{'R'} = ["A","G"];
$iupac{'r'} = ["a","g"];
$iupac{'Y'} = ["C","T"];
$iupac{'y'} = ["c","t"];
$iupac{'S'} = ["G","C"];
$iupac{'s'} = ["g","c"];
$iupac{'W'} = ["A","T"];
$iupac{'w'} = ["a","t"];
$iupac{'K'} = ["G","T"];
$iupac{'k'} = ["g","t"];
$iupac{'M'} = ["A","C"];
$iupac{'m'} = ["a","c"];



open my $fh, '<', $candidates or die "ERROR: cant open for reading file $candidates: $!\n";

my @goodMutations;
my @badMutations;
my @triMutations;
my $c = 0;
while(<$fh>){
	$c++;
	print "\xd", "Processing line: $c";
	chomp;
	my ($scaff, $pos, $genos) = split /\s+/, $_;
	my @genotypes = split //, $genos;
	if($noisy) {print "scaff $scaff, pos $pos, geno $genos, alleles " ;}
	my %alleles = ();
	foreach my $i (@genotypes){
		
		my $base1 = $iupac{$i}[0];
		my $base2 = $iupac{$i}[1];
		if(!exists ($alleles{$base1})){$alleles{$base1} = 1;}
		if(!exists ($alleles{$base2})){$alleles{$base2} = 1;}
	}

    if($noisy) {foreach my $k (keys %alleles){print $k ," ";}}
    if($noisy) {print "\n";}

    # skip if SNP has more than 2 alleles
    if(scalar keys %alleles > 2){
    	push @triMutations, $_;
    	next;
    }
    
    # else check whether all samples have both alleles
    my $nsamplesWithBoth = 0;
    my @alleles = keys %alleles;
    foreach my $bam (@BAMS){
    	my $cmd = "samtools mpileup -r $reference $bam -r$scaff:$pos-$pos 2> /dev/null";
    		my $samout = `$cmd`;
       		my @fields = split /\t/,$samout;
    		if($noisy){print "--$samout-- SCALAR IS ", scalar @fields, "\n";}
    		my @basesInMpileup = split //, uc $fields[4];
    		if  ($noisy){print "BASES IN MPILEUP @basesInMpileup\n\n";}
    		if( $alleles[0] ~~ @basesInMpileup && $alleles[1] ~~ @basesInMpileup ){
    			$nsamplesWithBoth++;
    		}

    }
    if($nsamplesWithBoth < 1 ){ die "ERROR: alleles not found in any sample for line $_\n"; }
    if($nsamplesWithBoth == 1){ push @goodMutations, $_; }
    if($nsamplesWithBoth > 1){ push @badMutations, $_; }
}

print "\nChecked ", scalar @goodMutations + scalar @badMutations + scalar @triMutations, ". Pass: ", scalar @goodMutations , ", Bad: ", scalar @badMutations, ", Tri: ", scalar @triMutations, "\n";

open my $fh1, '>', $candidates.".good.txt" or die "ERROR: cant open for writing file $candidates.good.txt: $!\n";
open my $fh2, '>', $candidates.".bad.txt" or die "ERROR: cant open for writing file $candidates.bad.txt: $!\n";
open my $fh3, '>', $candidates.".tri.txt" or die "ERROR: cant open for writing file $candidates.tri.txt: $!\n";
foreach (@goodMutations){
	print $fh1 $_, "\n";
}

foreach (@badMutations){
	print $fh2 $_, "\n";
}

foreach (@triMutations){
	print $fh3 $_, "\n";
}
close($fh1);close($fh2);close($fh3);
