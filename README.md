# MutationAccumulation  
  
Scripts used in Fazalova & Nevado (submitted) "Low spontaneous mutation rate and Pleistocene radiation of pea aphids"  
  
  
Folders:  
  
divergenceEstimates  
  
  	Scripts used to assign genomic regions to X/Autosomes, identify CDSs, and obtain 4-fold degenerate sites.  
  	Includes:  
  	getCDS.pl (generates list of autosomal and X-linked genomic regions according to Jaquiery et al. 2019)  
  	getCDS_2.pl (takes list of intervals, aligned fasta files for each scaffold, and output CDSs for each gene)  
  	mergeExons.pl (takes information from 2 previous steps and generates single CDS for each gene by merging all exons)  
  
  
mutationAccumulation  
  
  	Scrips used to identify and filter putative de novo mutations.
  	getCandidateMutations (generate list of candidate mutation from vcfs of 2 samples. C++ program, use make to install then run it without arguments for further info)  
  	filterCandidateMutations.pl (takes output from getCandidateMutations and filters mutations using bam files from the 2 samples)  
  	