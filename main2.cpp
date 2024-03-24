// include header
#include "planet.hpp"
#include "system.hpp"
#include "header.hpp"

// useful line for debugging:
//std::cout << __FILE__ << __LINE__<< std::endl;

int main(){


int step_size_rel_to_orbit = 10000; //50000 
double h = 2*M_PI/step_size_rel_to_orbit;
//timestep h
//even if timestep is adaptive, this is the initial step

int num_orbits = 10;
double t_max = num_orbits*2*M_PI;  //10*2pi = 10 orbits

int iter = (int) t_max/h;



// planet parameters 
double m1 = 1;
double m2 = 1e-3;
double e = 0.5;
double a = 1;

// system parameters
double j = 0;


// printing parameters
double width = 7;
double prec = 7;

//@@ -28,9 +36,10 @@
System Sosy;
Sosy.initialize_kepler_orbit(e, a, m1 ,m2);
Sosy.coord_transf();
bool adaptive = true;
//fixed time_step or adaptive?
double eta = 0.001;
//constant for calculation of adaptive time step

std::string integrator = "Euler";
// open file .txt to save data
std::ofstream ofile;
//ofile.open("test3.txt");
std::string file_dir ="D:/drive dateien/Uni/Semester_7_WiSe23_24/Astro-F_Praktikum/Chaos/data/";
//change directory to your liking
//This line does not give an error if your adress is non-existent
std::string file_name = "2body" + integrator+ "_" + int_to_str(step_size_rel_to_orbit)+ "_steps_per_Orbit";
if (adaptive== true){
     file_name = file_name + "_adaptive";
}
std::string file_end = ".txt";
std::string file = file_dir + file_name +file_end;

ofile.open(file);
//initialize time t to zero
double t=0.;
// save to file data: t p1 x_1 y_1 z_1 vx_1 vy_1 vz_1 p2 x_2 y_2 z_2 vx_2 vy_2 vz_2 dt
for(int i = 0; t<t_max; i++){
     //break when we reach the given max time instead of #iterations, this opens up the adaptable time_step
     //the index i is probably redundant now
     ofile  << scientific_format(t, width, prec);

     for(int i_planet = 0 ; i_planet < Sosy.planets.size(); i_planet++){
          ofile   << " " << i_planet+1
                    << " " << scientific_format(Sosy.planets.at(i_planet).r(0), width, prec)
                    << " " << scientific_format(Sosy.planets.at(i_planet).r(1), width, prec)
                    << " " << scientific_format(Sosy.planets.at(i_planet).r(2), width, prec)
                    << " " << scientific_format(Sosy.planets.at(i_planet).v(0), width, prec)
                    << " " << scientific_format(Sosy.planets.at(i_planet).v(1), width, prec)
                    << " " << scientific_format(Sosy.planets.at(i_planet).v(2), width, prec);

    }
    ofile  << " "<< scientific_format(Sosy.compute_eccentricity(), width, prec);
    ofile  << " "<< scientific_format(Sosy.compute_energy(), width, prec);
    ofile  << " "<< scientific_format(Sosy.compute_semi_maj_ax(), width, prec);
    ofile  << " "<< scientific_format(norm(Sosy.compute_spec_ang_mom()), width, prec);
    ofile  << " "<< scientific_format(h, width, prec); //added timestep to the text file
    ofile  << std::endl;

     //calculate variable time_step if appropriate
     if(t!=0){
          if (adaptive== true){
               double new_h= Sosy.adaptive_time_step(eta, h);
               //calc. variable time step h
               if (new_h !=0){
                    if(new_h != eta/h){
                         h = new_h;
                         //only use it if it's non-zero, a_dot is non-zero(eta/h) and it's not the first step
                    }
               }
          }
     }
   


     
   // Evolve the system integrator of choice 
   if (integrator == "Euler"){
       Sosy.evolveEuler(h);
   }
   else if(integrator == "EulerCromer"){
        Sosy.evolveEulerCromer(h);
   }
   else if (integrator == "LeapFrog"){
        Sosy.evolveLeapFrog(h, i, iter);
   }
   else if (integrator == "RK4"){
        Sosy.evolve_RK4(h);
   }

   t= t+h;
   //track time explicitly
     

}
    // close file
    ofile.close();
    return 0;
}