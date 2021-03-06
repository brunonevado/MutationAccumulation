// Author: Bruno Nevado
//  Created by bnevado on 03/10/2014.
//

#ifndef __vcf2fas__vcf__
#define __vcf2fas__vcf__

#include <stdio.h>

#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <cstring> // for strncmp
#include <algorithm>

#include "common.h"


struct vcf_line {
    int start, end;
    char alt;
    bool ref;
    
};

struct vcf_region {
    std::vector < vcf_line > lines;
    std::string contig;
};


class vcf {

private:
    int failed_parse;
    std::string infile;
    std::string ind_name;
    std::map < std::string, vcf_region > contigs;
    static std::string get_genotype( const std::vector <std::string > & fields);
public:
    vcf ( std::string infile );
    void readfile();
    std::string get_infile () { return infile; };
    std::string get_ind_name () { return ind_name; };
    static char toIUPAC (const std::string instr);
    int get_n_failed_lines ( ) { return failed_parse; };
    std::string make_fas ( std::string contig, std::string reference );
    int get_total_variants_read();
    void speak();
};




#endif /* defined(__vcf2fas__vcf__) */
