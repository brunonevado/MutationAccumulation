use strict;
use warnings;

# =======================================================================
# This script should be ran after getCDS_1.pl
# It will read the resulting "blastlike" file (either X or A, given as argument)
#     and print out fasta files for each CDS.
# It expects fasta files named scaffold.fas in folder infolder
# Scaffold files should be desinterleaved and be aligned (doesnt check)
# =======================================================================

# Author: Bruno Nevado

my $verbose = 1;

my $infile = "out.A.blastlike";
my $infolder = "../FASTA/";
my $outfolder = "./A";

my %comp;
$comp{'A'} = 'T';
$comp{'T'} = 'A';
$comp{'C'} = 'G';
$comp{'G'} = 'C';

$comp{'R'} = 'Y';
$comp{'Y'} = 'R';
$comp{'S'} = 'S';
$comp{'W'} = 'W';
$comp{'K'} = 'M';
$comp{'M'} = 'K';

$comp{'N'} = 'N';



my %isStopCodon;
$isStopCodon{'TAG'} = 1;
$isStopCodon{'TAA'} = 1;
$isStopCodon{'TGA'} = 1;


# =======================================================================
# getFasta
# =======================================================================

sub readFasta {
	# returns fasta file as @AoA plus @names
	# from folder $scaffold.fas in folder $infolder
	my ($scaffold, $folder ) = @_;

	my @names = ();
	my @AoA = ();
	my $temp_seqname = "";
	print "Reading scaffold $scaffold from folder $folder ...\n" if $verbose;

    open my $fhgen, '<', $folder."/".$scaffold.".fas" or die "ERROR: cant open genome file $scaffold in folder $infolder: $!\n";

  	while(<$fhgen>){
    	chomp;
    	if(s/>//){
    		s/>//;
    		$temp_seqname = $_;
    		push @names, $_;

    	}
    	else{
    		push @AoA, [split //, $_];
        	}
	}
	close($fhgen);
	print "  done, ", scalar @AoA, " sequences, ",  scalar @{$AoA[0]},  " bp long\n" if $verbose;

	return (\@names, \@AoA);
}




# =======================================================================
# And now get CDS
# =======================================================================


open my $fhcds, '<', $infile or die "ERROR: cant open for reading infile $infile: $!\n";


my $current_scaffold = "";
my %assignment = ();
my $ref_names;
my $ref_AoA;
while(<$fhcds>){
	chomp;
	my @fields = split /\t/,$_;
    my $scaffold = $fields[0];
	my $start = ($fields[1] > $fields[2]) ? $fields[2] : $fields[1];
	my $stop = ($fields[1] < $fields[2]) ? $fields[2] : $fields[1];
	my $reverse = ($fields[1] < $fields[2]) ? 0 : 1;
	

	print "Looking for $scaffold from $start to $stop, is reversed: $reverse\n" if $verbose > 1;

	if( $scaffold ne $current_scaffold){
		($ref_names, $ref_AoA) = readFasta($scaffold, $infolder );
		$current_scaffold = $scaffold;
	}

	# print CDSs
	my $ofile = $outfolder."/".$scaffold.".".$start."_".$stop.".$reverse.fas";
    open my $fh, '>', $ofile or die "ERROR: cant open for writing file $ofile: $!\n";



	if($reverse){
    	for my $ind (0..scalar @{$ref_AoA} -1 ){
    		print $fh ">",$$ref_names[$ind], "\n";
    		for (my $i = $stop ; $i > $start; $i--){
				print $fh $comp{ $ref_AoA->[$ind]->[$i-1] };
				}
			# add gaps at end if needed
			my $length = $stop - $start ;
			while( $length % 3 ){
				print $fh '-';
				$length++;

			}
			print $fh "\n";
    	}
    }

    else{
		for my $ind (0..scalar @{$ref_AoA} -1){
			print $fh ">",$$ref_names[$ind], "\n";
			for my $i ($start..$stop){
				print $fh $ref_AoA->[$ind]->[$i-1] ;

			}
			# add gaps at end if needed
			my $length = $stop - $start + 1;
			while( $length % 3 ){
				print $fh '-';
				$length++;

			}			
	    	print $fh "\n";
		}
	}
	
	close($fh);
    	
}





