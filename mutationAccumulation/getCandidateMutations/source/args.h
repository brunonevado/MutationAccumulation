/*
 B. Nevado, 01/10/2014
*/



#ifndef __args__
#define __args__

#include <stdio.h>
#include <vector>
#include <string>
#include <iostream> 


struct sargs {
    std::vector <std::string > args_string;
    std::vector <int > args_int;
    std::vector <float > args_float;
    std::vector <bool> args_booleans;
    std::vector <std::string > args_string_optional;
    std::vector <int > args_int_optional;
    //std::vector <float > args_float_optional;
    //std::vector <bool> args_booleans_optional;
    
};


class args {
public:
    sargs static getargs ( int argc, const char * argv[],
                          const std::vector <std::string > & args_string,
                          const std::vector <std::string > & args_bool,
                          const std::vector <std::string > & args_int,
                          const std::string & args_string_opt,
                          const std::string & args_int_opt
           
                          );
    
    
};


#endif /* defined(__args__) */
