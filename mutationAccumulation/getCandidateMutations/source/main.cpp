/*

 Uses BIO++
 
 */

#include <iostream>
#include <fstream>      // std::ofstream

#include <Bpp/Seq/Alphabet.all> /* this includes all alphabets in one shot */
#include <Bpp/Seq/Container.all> /* this includes all containers */
#include <Bpp/Seq/Io.all> /* this includes all sequence readers and writers */
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Seq/SiteTools.h>

#include "args.h"
#include "vcf.h"

// 25022019: this code was modified from vcf2fas to identify candidate mutations from 2 vcf files
// 04022019: cleaned up code and added help section

void help(){
    std::cout << "#################################\n  getCandidateMutations 04022020 \n#################################" << std::endl << std::endl;
    std::cout << "Create list of candidate de novo mutations from genotype calls for " << std::endl;
    std::cout << "2 samples (first and last generation in the MA experiment)." << std::endl <<  std::endl;;
    std::cout << "Usage: "<< std::endl;
    std::cout << " getCandidateMutations -reference Scaffold1.fa -vcfs vcfs.txt -outfile out.txt" << std::endl<< std::endl;
    std::cout << "   -reference: reference genome used." << std::endl;
    std::cout << "   -vcfs: text file with path to the two vcf files to use." << std::endl;
    std::cout << "   -outfile: text file with positions of putative mutations." << std::endl<< std::endl;
    
    std::cout << "Requirements: " << std::endl;
    std::cout << "  BIO++" << std::endl;
    std::cout << "  VCF files must be obtained with SAMtools v1.0+, and filtered with bcftools " <<std::endl;
    std::cout << "    (should contain 'PASS' for confident variants)." << std::endl;
    std::cout << "  Reference homozygous calls should be present." << std::endl;
    std::cout << "Notes: INDELS are ignored." << std::endl;
    
    
}

int main(int argc, const char * argv[]) {
    
    bool outputFasta = false;
    int nSitesCoveredInBoth = 0; // might not be enough for really long genomes
    int nMutations = 0;
    
    sargs myargs;
    try{
        myargs = args::getargs(argc, argv, std::vector<std::string> {"reference","vcfs","outfile"}, std::vector<std::string> {}, std::vector<std::string>  {}, std::string {"contigs"}, std::string {}); }
    catch (std::string e){
        std::cout << " Args failed: " << e << std::endl;
        help();
        exit(1);
    }
    
    std::string reference = myargs.args_string.at(0);
    std::string infilesVCF = myargs.args_string.at(1);
    std::string outfileMut = myargs.args_string.at(2);
    std::string infileCONTIGS = (myargs.args_string_optional.size() > 0) ? myargs.args_string_optional.at(0) : "";
    
    // READ VCFS
    std::vector < vcf > vec_vcfs;
    std::ifstream fh_vcfs ( infilesVCF );
    if( !fh_vcfs.is_open() ){
        std::cerr << "ERROR: Unable to open for reading infile " << infilesVCF << std::endl;
        exit(1);
    }
    std::string cfile;
    while (getline(fh_vcfs, cfile)) {
        vcf avcf(cfile);
        
        try {
            avcf.readfile();
        }
        catch(std::string e){
            std::cerr << "ERROR READING VCF FILE " << cfile << " : " <<  e << std::endl;
            exit(1);
        }
        vec_vcfs.push_back(avcf);
        std::clog << "Finished parsing infile " << cfile << ", ambiguous/total lines: " << avcf.get_n_failed_lines() << "/" << avcf.get_total_variants_read() << std::endl;
    }
    
    if( vec_vcfs.size() != 2 ){
        std::cerr << "ERROR: Need exactely two vcf files in infile " << infilesVCF << std::endl;
        exit(1);
    }
    // READ REF
    bpp::Fasta fasReader(-1);
    fasReader.strictNames(true);
    bpp::OrderedSequenceContainer *reference_genome = fasReader.readSequences(reference, &bpp::AlphabetTools::DNA_ALPHABET);
    std::clog << "Read " << reference_genome->getNumberOfSequences() << " contigs form reference file " << reference << std::endl;
    
    std::ofstream fho;
    fho.open(outfileMut);
    if(!fho.is_open()){
        std::cerr << "ERROR: cant open for writing file " << outfileMut << std::endl;
        exit(1);
    }
    
    if( infileCONTIGS == "" ){
        for ( unsigned int icontig = 0 ; icontig < reference_genome->getNumberOfSequences() ; icontig++ ) {
            bpp::AlignedSequenceContainer bpp_contig ( &bpp::AlphabetTools::DNA_ALPHABET );
            
            try{
                for (unsigned int iind = 0; iind < vec_vcfs.size(); iind++) {
                    std::string temp = vec_vcfs.at(iind).make_fas(reference_genome->getSequencesNames().at(icontig), reference_genome->getSequence(icontig).toString());
                    bpp::Sequence *sequence = new bpp::BasicSequence( vec_vcfs.at(iind).get_ind_name() , temp, &bpp::AlphabetTools::DNA_ALPHABET);
                    bpp_contig.addSequence(*sequence);
                    delete  sequence;
                }
            }
            catch(...){
                std::cerr << "ERROR: problems doing contig " << reference_genome->getSequencesNames().at(icontig) << std::endl;
                exit(1);
            }
            
            if( outputFasta == true ){
                try{
                    fasReader.writeSequences(std::string( reference_genome->getSequencesNames().at(icontig) + ".fas"  ), bpp_contig );
                }
                catch(...){
                    std::cerr << "ERROR: problems writing to file contig " << reference_genome->getSequencesNames().at(icontig) << std::endl;
                    exit(1);
                }
            }
            else{
                bpp::SiteContainer *vsc2 = new bpp::VectorSiteContainer(bpp_contig);
                for( unsigned int isite = 0; isite < vsc2->getNumberOfSites() ; isite++ ){
                    //nSitesCoveredInBoth
                    //nMutations
                    int a = vsc2->getSite(isite).getValue(0);
                    int b = vsc2->getSite(isite).getValue(1);
                    
                    if(a < 14 && b < 14) {
                        nSitesCoveredInBoth++;
                       if( a != b ){
                            nMutations++;
                            fho << reference_genome->getSequencesNames().at(icontig) << "\t" << vsc2->getSite(isite).getPosition() << "\t"<< vsc2->getSite(isite).toString() << std::endl;
                        }
                    }
                }
            }
            if( icontig % 100 == 0 ){
                bpp::ApplicationTools::displayGauge(icontig/100, (reference_genome->getNumberOfSequences()/100)-1, '.' );
            }
        }
    }
    else{
        std::cerr << "Contigs option not available yet!" << std::endl;
        exit(1);
    }
    fho.close();
    std::clog  << std::endl << "Number of sites covered in both samples: " << nSitesCoveredInBoth << std::endl;
    std::clog << "Number of putative mutations: " << nMutations << std::endl;
    return 0;
}
