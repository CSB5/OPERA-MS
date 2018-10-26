#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string.h>
#include "contig_reader.h"

#include "sigma.h"

ContigReader::~ContigReader() {}

AllReader::AllReader(){}

long int AllReader::read(const char* contigs_file, ContigMap* contigs) {
    int length = 0;
    long int assembly_size = 0;
    std::string id; 

    std::ifstream contigs_fp(contigs_file);
    
    if (contigs_fp.is_open()) {
        while(true){

           std::string line;

            if(!std::getline(contigs_fp, line)){
                if (length !=0 && length >= Sigma::contig_len_thr){
                    contigs -> insert(std::make_pair(id, new Contig(id, length)));
                //    std::cerr << length << "\n";
                //    std::cerr << id <<"\n";
                }

                assembly_size += length;
                break;
            }

            if(line[0] != '>'){
                length += (int) line.size();
            }

            else{

                assembly_size += length;
                if (length != 0 && length >= Sigma::contig_len_thr){
                    contigs -> insert(std::make_pair(id, new Contig(id, length)));
                //    std::cerr << id << "\n";
                 //   std::cerr << length << "\n";

                }
                id = line.substr(1, line.find(' ') - 1);
                length = 0;   
             } 
        }
        contigs_fp.close();
    } 
    
    else {
        fprintf(stderr, "Error opening file: %s\n", contigs_file);
        exit(EXIT_FAILURE);
    }

    std::string assembly_size_name = Sigma::output_dir + "/assembly_size.dat";          
    FILE* assembly_size_file = fopen(assembly_size_name.c_str(), "w");
    if (assembly_size_file != NULL) {
        fprintf(assembly_size_file, "%ld\n", assembly_size);
        fclose(assembly_size_file);
    }
    
    
    return assembly_size;
   
}

void AllReader::get_assembly_size(const char* contigs_file){
    int length = 0;
    long int assembly_size = 0;
    long int assembly_nb_contig = 0;
    std::ifstream contigs_fp(contigs_file);
    
    if (contigs_fp.is_open()) {
        std::string line; 
      
         while(true){

           std::string line;

           if(!std::getline(contigs_fp, line)){
               assembly_size+= length;
               break;
           }

           if(line[0] != '>'){
               length += (int) line.size();
           }

           else{
               assembly_size += length;
               assembly_nb_contig++;
               
               length = 0;   
           } 
         }

        contigs_fp.close();
    } 
    
    else {
        fprintf(stderr, "Error opening file: %s\n", contigs_file);
        exit(EXIT_FAILURE);
    }

    std::string assembly_size_name = Sigma::output_dir + "/assembly_size.dat";          
    FILE* assembly_size_file = fopen(assembly_size_name.c_str(), "w");
    if (assembly_size_file != NULL) {
        fprintf(assembly_size_file, "%ld\n", assembly_size);
        fprintf(assembly_size_file, "%ld\n", assembly_nb_contig);
        fclose(assembly_size_file);
    }

    
    
    Sigma::total_assembly_size = assembly_size;
    Sigma::total_assembly_nb_contig = assembly_nb_contig;

}


SOAPdenovoReader::SOAPdenovoReader() {}

long int SOAPdenovoReader::read(const char* contigs_file, ContigMap* contigs) {
    char id[256];
    int length;
    long int assembly_size = 0;

    FILE* contigs_fp = fopen(contigs_file, "r");

    if (contigs_fp != NULL) {
        while (!feof(contigs_fp)) {
            // >[ID] length [LENGTH] cvg_[COVERAGE]_tip_[TIP]\n
            if (fscanf(contigs_fp, ">%s %*s %d %*s\n", id, &length) == 2) {
                assembly_size += length;
                if (length >= Sigma::contig_len_thr) {
                    contigs->insert(std::make_pair(id, new Contig(id, length)));
                }
            } else {
                if( fscanf(contigs_fp, "%*[^\n]\n") );
            }
        }

        fclose(contigs_fp);
    } else {
        fprintf(stderr, "Error opening file: %s\n", contigs_file);
        exit(EXIT_FAILURE);
    }

    std::string assembly_size_name = Sigma::output_dir + "/assembly_size.dat";          
    FILE* assembly_size_file = fopen(assembly_size_name.c_str(), "w");
    if (assembly_size_file != NULL) {
        fprintf(assembly_size_file, "%ld\n", assembly_size);
        fclose(assembly_size_file);
    }
    
    
    return assembly_size;
}

void SOAPdenovoReader::get_assembly_size(const char* contigs_file) {
    int length;
    long int assembly_size = 0;
    long int assembly_nb_contig = 0;
    
    FILE* contigs_fp = fopen(contigs_file, "r");

    if (contigs_fp != NULL) {
        while (!feof(contigs_fp)) {
            // >[ID] length [LENGTH] cvg_[COVERAGE]_tip_[TIP]\n
            if (fscanf(contigs_fp, ">%*s %*s %d %*s\n", &length) == 1) {
                assembly_size += length;
                assembly_nb_contig++;
            } else {
                if( fscanf(contigs_fp, "%*[^\n]\n") );
            }
        }

        fclose(contigs_fp);
    } else {
        fprintf(stderr, "Error opening file: %s\n", contigs_file);
        exit(EXIT_FAILURE);
    }

    std::string assembly_size_name = Sigma::output_dir + "/assembly_size.dat";          
    FILE* assembly_size_file = fopen(assembly_size_name.c_str(), "w");
    if (assembly_size_file != NULL) {
        fprintf(assembly_size_file, "%ld\n", assembly_size);
        fprintf(assembly_size_file, "%ld\n", assembly_nb_contig);
        fclose(assembly_size_file);
    }
    
    Sigma::total_assembly_size = assembly_size;
    Sigma::total_assembly_nb_contig = assembly_nb_contig;

    //return assembly_size;
}

RAYReader::RAYReader() {}

long int RAYReader::read(const char* contigs_file, ContigMap* contigs) {
    char id[256];
    int length;
    long int assembly_size = 0;

    FILE* contigs_fp = fopen(contigs_file, "r");

    if (contigs_fp != NULL) {
        while (!feof(contigs_fp)) {
            // >[ID] [LENGTH] nucleotides\n
            if (fscanf(contigs_fp, ">%s %d %*s\n", id, &length) == 2) {
                assembly_size += length;
                if (length >= Sigma::contig_len_thr) {
                    contigs->insert(std::make_pair(id, new Contig(id, length)));
                }
            } else {
                if( fscanf(contigs_fp, "%*[^\n]\n") );
            }
        }

        fclose(contigs_fp);
    } else {
        fprintf(stderr, "Error opening file: %s\n", contigs_file);
        exit(EXIT_FAILURE);
    }

    std::string assembly_size_name = Sigma::output_dir + "/assembly_size.dat";          
    FILE* assembly_size_file = fopen(assembly_size_name.c_str(), "w");
    if (assembly_size_file != NULL) {
        fprintf(assembly_size_file, "%ld\n", assembly_size);
        fclose(assembly_size_file);
    }
    
    
    return assembly_size;
}

void RAYReader::get_assembly_size(const char* contigs_file) {
    int length;
    long int assembly_size = 0;
    long int assembly_nb_contig = 0;
    
    FILE* contigs_fp = fopen(contigs_file, "r");

    if (contigs_fp != NULL) {
        while (!feof(contigs_fp)) {
            // >[ID] [LENGTH] nucleotides\n
            if (fscanf(contigs_fp, ">%*s %d %*s\n", &length) == 1) {
                assembly_size += length;
                assembly_nb_contig++;
            } else {
                if( fscanf(contigs_fp, "%*[^\n]\n") );
            }
        }

        fclose(contigs_fp);
    } else {
        fprintf(stderr, "Error opening file: %s\n", contigs_file);
        exit(EXIT_FAILURE);
    }

    std::string assembly_size_name = Sigma::output_dir + "/assembly_size.dat";          
    FILE* assembly_size_file = fopen(assembly_size_name.c_str(), "w");
    if (assembly_size_file != NULL) {
        fprintf(assembly_size_file, "%ld\n", assembly_size);
        fprintf(assembly_size_file, "%ld\n", assembly_nb_contig);
        fclose(assembly_size_file);
    }
    
    
    Sigma::total_assembly_size = assembly_size;
    Sigma::total_assembly_nb_contig = assembly_nb_contig;

    //return assembly_size;
}



VelvetReader::VelvetReader() {}

long int VelvetReader::read(const char* contigs_file, ContigMap* contigs) {
    char id[256];
    int length;
    long int assembly_size = 0;

    FILE* contigs_fp = fopen(contigs_file, "r");

    if (contigs_fp != NULL) {
        while (!feof(contigs_fp)) {
            // >NODE_[ID]_length_[LENGTH]_cov_[COVERAGE]\n
            if (fscanf(contigs_fp, ">%s\n", id) == 1 && sscanf(id, "%*[^_]_%*[^_]_%*[^_]_%d_%*s", &length) == 1) {
                assembly_size += length;
                if (length >= Sigma::contig_len_thr) {
                    contigs->insert(std::make_pair(id, new Contig(id, length)));
                }
            } else {
                if( fscanf(contigs_fp, "%*[^\n]\n") );
            }
        }

        fclose(contigs_fp);
    } else {
        fprintf(stderr, "Error opening file: %s\n", contigs_file);
        exit(EXIT_FAILURE);
    }


    std::string assembly_size_name = Sigma::output_dir + "/assembly_size.dat";          
    FILE* assembly_size_file = fopen(assembly_size_name.c_str(), "w");
    if (assembly_size_file != NULL) {
        fprintf(assembly_size_file, "%ld\n", assembly_size);
        fclose(assembly_size_file);
    }
    

    return assembly_size;
}

void VelvetReader::get_assembly_size(const char* contigs_file) {
    char id[256];
    int length;
    long int assembly_size = 0;
    long int assembly_nb_contig = 0;
    
    FILE* contigs_fp = fopen(contigs_file, "r");

    if (contigs_fp != NULL) {
        while (!feof(contigs_fp)) {
            // >NODE_[ID]_length_[LENGTH]_cov_[COVERAGE]\n
            if (fscanf(contigs_fp, ">%s\n", id) == 1 && sscanf(id, "%*[^_]_%*[^_]_%*[^_]_%d_%*s", &length) == 1) {
                assembly_size += length;
                assembly_nb_contig++;
            } else {
                if( fscanf(contigs_fp, "%*[^\n]\n") );
            }
        }

        fclose(contigs_fp);
    } else {
        fprintf(stderr, "Error opening file: %s\n", contigs_file);
        exit(EXIT_FAILURE);
    }


    std::string assembly_size_name = Sigma::output_dir + "/assembly_size.dat";          
    FILE* assembly_size_file = fopen(assembly_size_name.c_str(), "w");
    if (assembly_size_file != NULL) {
        fprintf(assembly_size_file, "%ld\n", assembly_size);
        fclose(assembly_size_file);
    }
    
    Sigma::total_assembly_size = assembly_size;
    Sigma::total_assembly_nb_contig = assembly_nb_contig;
    
    //return assembly_size;
}
