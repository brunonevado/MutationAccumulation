// Author: Bruno Nevado
//  Created by bnevado on 03/10/2014.
//

#include "vcf.h"


vcf::vcf ( std::string in ){
    infile = in;
    failed_parse = 0;
}

void vcf::readfile () {
    std::ifstream fh_infile;
    std::string cline;
    
    fh_infile.open(infile);

    if(!fh_infile.is_open()){  throw std::string ( "Unable to open infile vcf ");}
    
    
    while( getline(fh_infile, cline) ){
        if( strncmp(cline.c_str(), "##", 2) == 0){
            continue;
        }
        else if( strncmp(cline.c_str(), "#CHR", 3) == 0){
            std::vector <std::string > fields = msplit(cline, "\t");
            if( fields.size() != 10 ){ throw std::string("more than 1 individual in infile vcf?");  }
            ind_name = fields.at(9);
            continue;
        }
        else{
            std::vector <std::string > fields = msplit(cline, "\t");
            vcf_line aline;
            try{
                if ( fields.at(6) != "PASS" || fields.at(7).substr(0,5) == "INDEL" ) { continue; }
                // homozygous ref call
                else if ( fields.at(7).substr(0,3) == "END" ){
                    aline.start = std::stoi(fields.at(1));
                    aline.end = std::stoi( fields.at(7).substr(4)  );
                    aline.ref = true;
                    
                    if( contigs.count(fields.at(0)) == 0  ){
                        vcf_region new_region;
                        new_region.contig = fields.at(0);
                        new_region.lines.push_back(aline);
                        contigs.insert( std::pair <std::string, vcf_region> (fields.at(0), new_region));
                        
                    }
                    else{
                        contigs.at(fields.at(0)).lines.push_back(aline);
                    }
                    
                    
                    
                    continue;
                }
                // SNP call
                else{
                    //std::cout << cline << std::endl;;
                    
                    std::string genotype = get_genotype(fields);
                    
                    //std::cout << genotype << std::endl;
                    
                    
                    aline.start = std::stoi(fields.at(1));
                    aline.end = std::stoi( fields.at(1)  );
                    //std::cout << cline << std::endl;
                    aline.alt = toIUPAC(genotype);
                    aline.ref = false;
                    if( contigs.count(fields.at(0)) == 0  ){
                        vcf_region new_region;
                        new_region.contig = fields.at(0);
                        new_region.lines.push_back(aline);
                        contigs.insert( std::pair <std::string, vcf_region> (fields.at(0), new_region));
                        
                    }
                    else{
                        contigs.at(fields.at(0)).lines.push_back(aline);
                    }
                    continue;
                }
            }
            catch(...){
                //std::clog <<"Failed to parse VCF line: "  <<  cline ;
                failed_parse++;
                //throw "Failed to parse VCF line: " + cline ;
            }
            
        }
    }
    fh_infile.close();
}


std::string vcf::get_genotype ( const std::vector <std::string > & fields){

    // alternative alleles
    std::vector <std::string> alt_alleles = msplit( fields.at(4), "," );
    std::vector <std::string> all_alleles;
    all_alleles.push_back(fields.at(3));
    for(auto i:alt_alleles){
        all_alleles.push_back(i);
    }
    
    msplit( fields.at(4), "," );

   // check where GT is stored
    unsigned int gt_field=0;
    bool found_gt = false;
    std::vector <std::string> keys = msplit( fields.at(8), ":" );
    for(gt_field=0; gt_field < keys.size(); gt_field++){
        if( strncmp(keys.at(gt_field).c_str(), "GT", 2) == 0 ){
            found_gt = true;
            break;
        }
    }
    
    if(!found_gt){ throw std::string("Cannot find GT field"); }
    
    std::vector <std::string> values = msplit( fields.at(9), ":" );
    std::vector <std::string> gt_values = msplit( values.at(gt_field), "/" );
    if(gt_values.size() != 2){ throw std::string("GT field does not have 2 values separated by /"); }

    std::string geno = all_alleles.at(std::stoi(gt_values.at(0))) + all_alleles.at(std::stoi(gt_values.at(1)));
    
    return geno;
    
}


char vcf::toIUPAC (const std::string instr){
    std::string in = instr;
    std::transform(in.begin(), in.end(), in.begin(), ::tolower);
    if( in.at(0) == in.at(1) ){
        return in.at(0);
    }
    else if ( (in.at(0) == 'a' && in.at(1) == 'c') || (in.at(1) == 'a' && in.at(0) == 'c') ){
        return 'm';
    }
    else if ( (in.at(0) == 'a' && in.at(1) == 'g') || (in.at(1) == 'a' && in.at(0) == 'g') ){
        return 'r';
    }
    else if ( (in.at(0) == 'a' && in.at(1) == 't') || (in.at(1) == 'a' && in.at(0) == 't') ){
        return 'w';
    }
    else if ( (in.at(0) == 'c' && in.at(1) == 'g') || (in.at(1) == 'c' && in.at(0) == 'g') ){
        return 's';
    }
    else if ( (in.at(0) == 't' && in.at(1) == 'g') || (in.at(1) == 't' && in.at(0) == 'g') ){
        return 'k';
    }
    else if ( (in.at(0) == 'c' && in.at(1) == 't') || (in.at(1) == 'c' && in.at(0) == 't') ){
        return 'y';
    }
    else{
        std::cerr << "ERROR (turnIUPAC): no code available for " << in.at(0) << in.at(1) << "\n";
        exit(1);
    }
}

std::string vcf::make_fas ( std::string contig_name, std::string reference_seq ){
    if( contigs.count(contig_name) == 0){
        std::string res( reference_seq.size(), 'N' );
        return res;
    }
    else{
        std::string res( reference_seq.size(), 'N' );
        for ( unsigned int ivar = 0; ivar < contigs.at(contig_name).lines.size()  ; ivar++ ) {
            if( contigs.at(contig_name).lines.at(ivar).ref == true ){
                res.replace(contigs.at(contig_name).lines.at(ivar).start - 1
                            ,  contigs.at(contig_name).lines.at(ivar).end - contigs.at(contig_name).lines.at(ivar).start + 1
                            , reference_seq.substr( contigs.at(contig_name).lines.at(ivar).start - 1, contigs.at(contig_name).lines.at(ivar).end - contigs.at(contig_name).lines.at(ivar).start + 1 ) );
            }
            else{
                res.at( contigs.at(contig_name).lines.at(ivar).start - 1  ) = contigs.at(contig_name).lines.at(ivar).alt;
            }
        }
        return res;
    }
    
}

int vcf::get_total_variants_read(){
    int tot = 0;
    std::map< std::string, vcf_region >::iterator it = contigs.begin();
    for (it=contigs.begin(); it!=contigs.end(); ++it)
        tot += it->second.lines.size();
    return tot;
}


void vcf::speak(){
    std::cout << "Infile: " << infile << std::endl << "N contigs: " << contigs.size() << std::endl << "N vars per contig:" << std::endl;
    std::map< std::string, vcf_region >::iterator it = contigs.begin();
    for (it=contigs.begin(); it!=contigs.end(); ++it)
        std::cout << it->first << " : " << it->second.lines.size() << '\n';
}
