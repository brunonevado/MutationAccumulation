use strict;
use warnings;
use Carp;

# Author: Bruno Nevado

my $cds = "cdss.gff";
my $folder = "A";
my $outfolder = "A_genes";

# needs get4foldSites (https://github.com/brunonevado/get4foldSites) 
# and concatenateFasta (https://github.com/brunonevado/concatenateFasta) in PATH

## read cds info
open my $fh, '<', $cds or die "ERROR: cant open cds file $cds: $!\n";

my %cdsinfo; # {scaf}{start} = 'geneName'
my %genes;	# {geneName} = @files

while(<$fh>){
	my @fields = split /\s+/,$_;
	my $scaf = $fields[0];
	my $start = $fields[3];
	$fields[8] =~ m/.+GeneID:(.+?),.+/;
	my $geneID=$1;
	$cdsinfo{$scaf}{$start} = $geneID;
	
}




opendir(my $dh, $folder) || die "ERROR: can't opendir $folder: $!";
my @files = grep { /.+\.fas/ && -f "$folder/$_" } readdir($dh);
closedir $dh;

print "Found ", scalar @files ," in folder $folder\n"; 
my $found = 0;
my $notfound = 0;
foreach my $f (@files){
	my @fields = split /\./, $f;
	my @fields2 = split /_/,$fields[1];
	my $scaffold = $fields[0];
	my $start = $fields2[0];
	my $thisfound = 0;
	for(my $i = $start -2 ; $i < $start + 2; $i++){
		if(exists $cdsinfo{$scaffold}{$i}){
			$thisfound = 1;
			my $gene = $cdsinfo{$scaffold}{$i};
			if(exists($genes{$gene})){
				push @{$genes{$gene}}, $f;
				}else{
				@{$genes{$gene}} = ($f); 
				}
			}
		}
	if($thisfound){$found++;}else{$notfound++;}
}

print "Found $found files, not found $notfound files\n";

print "Found ", scalar keys %genes , " genes\n";

my %nexons;
my $count = 0;
foreach my $g (keys %genes){
	$count++;
	my $nex = scalar @{$genes{$g}};
	$nexons{$nex}++;
	open my $fho, '>', "tempMergeCDS" or die "ERROR: cant open for writing tempMergeCDS: $!\n";
	foreach my $f ( @{$genes{$g}}){
		print $fho "$folder/$f\n";
		}
	close($fho);
	my $cmd =" concatenateFasta -missChar N -files tempMergeCDS -outfile $outfolder/gene.$g.fas 2> tempMergeCDS.log";
	
	system($cmd);
	open my $fh02, '>', "tempMergeCDS" or die "ERROR: cant open for writing tempMergeCDS: $!\n";
	print $fh02 "$outfolder/gene.$g.fas\n";
	close($fh02);
	my $cmd4fold = " get4foldSites -infile tempMergeCDS -outfile $outfolder/gene.$g.4fold.fas -iupac 1 -verbose 0";
	system($cmd4fold);
	print "\xd",  "Finished $count genes";
	}

foreach my $n (sort {$a <=> $b} keys %nexons){
	print $n, " ", $nexons{$n}, "\n";
	}




