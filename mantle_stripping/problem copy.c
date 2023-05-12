
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"


int main(int argc, char* argv[]){
    
    for (int n=1; n<2; n++){
        char arch_file[150];
        char inp_file[150];
        
        char arch_ext[10] = ".bin";
        char inp_ext[10] = ".txt";
        
        char arch_pw[100] = "archive";
        char inp_pw[150] = "../../mantle_stripping_input/2step_mantle_stripping_input";
        
        char file_no[2];
        sprintf(file_no, "%d", n);
        
        sprintf(arch_file, "%s%s%s", arch_pw, file_no, arch_ext);
        sprintf(inp_file, "%s%s%s", inp_pw, file_no, inp_ext);
        printf("%s", arch_file);
        
        struct reb_simulationarchive* archive = reb_open_simulationarchive(arch_file);
        struct reb_simulation* sim = reb_create_simulation_from_simulationarchive(archive, 0);
        FILE* of_msi = fopen(inp_file,"a+");
        
        for (int i=1;i<sim->N-2;i++){
            struct reb_orbit o =  reb_tools_particle_to_orbit(sim->G, sim->particles[i], sim->particles[0]);
            fprintf(of_msi, "%u\t", sim->particles[i].hash);
            fprintf(of_msi, "%e\t", sim->particles[i].m);
            //double CMF = 0.3;
            //double CMF = -0.10959*(o.a-0.35)+0.5; // line equation = CMF = -0.10959(a-0.35) + 0.5   for max CMF = 0.5 and min CMF = 0.1
            //double CMF = 0.6-(0.071701907*pow(o.a-0.35, 3./2.)); //power law equation for max CMF = 0.6 and min CMF = 0.1
            //double CMF = 0.7*exp(-0.5331260682*(o.a-0.35)); //exponential equation for CMF: max CMF = 0.7 and min CMF = 0.1
            //double CMF = 0.6*exp(-0.4908930053*(o.a-0.35)); //exponential equation for CMF: max CMF = 0.6 and min CMF = 0.1
            //fprintf(of_msi, "%e\t", CMF); //core frac of body
            
            /*if (o.a < 1.0){
                   fprintf(of_msi, "%e\t", 0.7); //core frac of body
                   }   else if(1.0 <= o.a && o.a <= 3.0){
                   fprintf(of_msi, "%e\t", 0.5); //core frac of body
                   }   else if(o.a > 3.0){
                   fprintf(of_msi, "%e\t", 0.3); //core frac of body
                   }*/
            
            if (o.a <= 2.0){
            fprintf(of_msi, "%e\t", 0.4); //core frac of body
            }   else if(2.0 < o.a ){
            fprintf(of_msi, "%e\t", 0.2); //core frac of body
            }
            
            fprintf(of_msi, "%e\t", o.a);
            fprintf(of_msi, "%e\t", o.e);
            fprintf(of_msi, "\n");
        }
        
        fclose(of_msi);
        
    }
    
    //FILE* of_dmsi = fopen("../../mantle_stripping_input/nu_mantle_stripping_input4.txt","a+");
    //FILE* of_dmsi = fopen("../../mantle_stripping_input/nu_mantle_stripping_input6.txt","a+");
    
    //struct reb_simulationarchive* archive = reb_open_simulationarchive("archive4.bin");
    //struct reb_simulation* sim = reb_create_simulation_from_simulationarchive(archive, 0);
        
    /*for (int i=0;i<sim->N;i++){
        struct reb_orbit o =  reb_tools_particle_to_orbit(sim->G, sim->particles[i], sim->particles[0]);
        fprintf(of_dmsi, "%u\t", sim->particles[i].hash);
        fprintf(of_dmsi, "%e\t", sim->particles[i].m);
        fprintf(of_dmsi, "%e\t", 0.3); //core frac of body
        fprintf(of_dmsi, "%e\t", 0.7); //mantle frac of body
        fprintf(of_dmsi, "%e\t", o.a);
        fprintf(of_dmsi, "%e\t", o.e);
        fprintf(of_dmsi, "\n");
    }*/
    
    /*for (int i=0;i<sim->N;i++){
        struct reb_orbit o =  reb_tools_particle_to_orbit(sim->G, sim->particles[i], sim->particles[0]);
        fprintf(of_dmsi, "%u\t", sim->particles[i].hash);
        fprintf(of_dmsi, "%e\t", sim->particles[i].m);

        if (o.a < 1.0){
        fprintf(of_dmsi, "%e\t", 0.5); //core frac of body
        fprintf(of_dmsi, "%e\t", 0.5); //mantle frac of body
        }   else if(1.0 <= o.a && o.a <= 3.0){
        fprintf(of_dmsi, "%e\t", 0.3); //core frac of body
        fprintf(of_dmsi, "%e\t", 0.7); //mantle frac of body
        }   else if(o.a > 3.0){
        fprintf(of_dmsi, "%e\t", 0.1); //core frac of body
        fprintf(of_dmsi, "%e\t", 0.9); //mantle frac of body
        }
        
        fprintf(of_dmsi, "%e\t", o.a); //mantle frac of body
        fprintf(of_dmsi, "%e\t", o.e); //mantle frac of body
        fprintf(of_dmsi, "\n");
    }
    
    fclose(of_dmsi);*/

    /*FILE* of_fp = fopen("../../final_orbital_parameters/final_orbital_parameters5.txt","a+");
    
    struct reb_simulation* r = reb_create_simulation_from_simulationarchive(archive, -1);
    
    for (int i=0;i<r->N;i++){
        struct reb_orbit o =  reb_tools_particle_to_orbit(r->G, r->particles[i], r->particles[0]);
        fprintf(of_fp, "%u\t", r->particles[i].hash);
        fprintf(of_fp, "%e\t", r->particles[i].m);
        fprintf(of_fp, "%e\t", o.a);
        fprintf(of_fp, "%e\t", o.e);
        fprintf(of_fp, "\n");
        }
        
        fclose(of_fp);

        reb_free_simulation(r);*/
    }


