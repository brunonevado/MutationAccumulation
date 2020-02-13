use strict;
use warnings;

# =======================================================================
# few things that are not checked:
# fasta file should be desinterleaved
# fasta file should contain only upper case characters
# assignment should be to 'X' or 'A', not checked 
# =======================================================================

my $verbose = 1;

my $infile_XA = "X_A_Scaffold.txt";
my $infile_genome = "genome.fasta"; # has to be desinterleaved!
my $infile_cds = "cdss.gff";

my $outfileX = "out.X.blastlike";
my $outfileA = "out.A.blastlike";


my %comp;
$comp{'A'} = 'T';
$comp{'T'} = 'A';
$comp{'C'} = 'G';
$comp{'G'} = 'C';

my %isStopCodon;
$isStopCodon{'TAG'} = 1;
$isStopCodon{'TAA'} = 1;
$isStopCodon{'TGA'} = 1;


sub getScaffold {
	# scaffold is name of scafold to get
	# infile is assembly (must be desinterleaved)
	# infile 2 is information on X/A, with 4 fields "scaffold	chromosome	start	end"
    my ($scaffold, $infile, $infile2 ) = @_;
	
	# =======================================================================
	# read genome
	# =======================================================================

	my %genome = ();
	my %genome_ass = ();
	my $temp_seqname;
	print "Reading genome from $infile looking for scaffold $scaffold...\n" if $verbose;

    open my $fhgen, '<', $infile or die "ERROR: cant open genome file $infile: $!\n";

  	while(<$fhgen>){
    	chomp;
    	if(s/>//){
    		$temp_seqname = $_;

    	}
    	elsif($temp_seqname eq $scaffold){
    		$genome{$temp_seqname} = [split //, $_];
    		my @temp_array = ('N') x length($_);
    		$genome_ass{$temp_seqname} = \@temp_array;
    		print "Scaffold found, ",  length($_) , " bp long\n" if $verbose ;
			last;
    	}
	}
	close($fhgen);

	# =======================================================================
	## get information on X vs A into %xa{scaffold}
	# =======================================================================

	print "Reading X/A assignment from $infile2... \n";
	open my $fhxa, '<', $infile2 or die "ERROR: cant find infile $infile2: $!\n";
    my $counter;
	while(<$fhxa>){
		next if /start\tend/;
		chomp;
		my @fields = split /\t/,$_;
		next unless $fields[0] eq $scaffold;
		$counter++;
		for my $i($fields[2] .. $fields[3]){
			$genome_ass{$fields[0]}[$i-1] = $fields[1];
		}
	} 
	close($fhxa);
	print "Found ", $counter, " assigned regions\n...\n" if $verbose;

	return %genome_ass;
    }




# =======================================================================
# And now get CDS
# =======================================================================


open my $fhcds, '<', $infile_cds or die "ERROR: cant open for reading gff file $infile_cds: $!\n";

open my $fhoX, '>', $outfileX or die "ERROR: Cant open for output file $outfileX: $!\n";
open my $fhoA, '>', $outfileA or die "ERROR: Cant open for output file $outfileA: $!\n";

my $current_scaffold = "";
my %assignment = ();

while(<$fhcds>){
	my @fields = split /\t/,$_;
	next unless $fields[2] eq "CDS";
    my $scaffold = $fields[0];
	my $start = ($fields[6] eq "-") ? $fields[4] : $fields[3];
	my $stop = ($fields[6] eq "+") ? $fields[4] : $fields[3];
	my $reverse = ($fields[6] eq "+") ? 0 : 1;
	my $add = $fields[7];

	if( $scaffold ne $current_scaffold){
		%assignment = getScaffold ( $scaffold, $infile_genome, $infile_XA );
		$current_scaffold = $scaffold;
	}

    my @thisCDS = ();
    if($reverse){
    	# adjust for frame position, down on rev strand
    	$start -= $add;
    }
    else{
    	# adjust for frame position, up on for strand
		$start += $add;
	}
	# now check where to print
    my $ass = $assignment{$scaffold}[$start-1];
    next if $ass eq 'N';
    my $printEntry = 1;
    for my $i ($start..$stop){
    	if( $assignment{$scaffold}[$i-1] ne $ass){
    		$printEntry = 0;
    		last;

    	}

    }
    next unless $printEntry;
    if($ass eq 'A'){
   		print $fhoA $scaffold,"\t", $start , "\t", $stop , "\n";
    	}else{
   		print $fhoX $scaffold,"\t", $start , "\t", $stop , "\n";

    	}
	
}





