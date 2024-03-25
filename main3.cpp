// include header
#include "planet.hpp"
#include "system.hpp"
#include "header.hpp"

// useful line for debugging:
//std::cout << __FILE__ << __LINE__<< std::endl;

int main(){


int step_size_rel_to_orbit = 10000; //50000 
double P = 2*M_PI; //Period P is 2pi for our settings
double h = P/step_size_rel_to_orbit;
//timestep h
//even if timestep is adaptive, this is the initial step

int num_orbits = 10;  //should be 
double t_max = num_orbits*P;  //10*2pi = 10 orbits

int iter = (int) t_max/h;



// planet parameters 
double m1 = 1;
double m2 = 1e-3;
double mu2 = m2/m1;   //mass ratio
double m3= 1e-5;
double mu3 = m3/m1;
double e = 0;
double a1 = 1;
double delta =0.01*a1;
double a2 = a1 + delta;


// system parameters
double j = 0; //spec. ang mom, not needed in this exercise


// printing parameters
double width = 7;
double prec = 7;

//@@ -28,9 +36,10 @@
System Sosy;
Sosy.initialize_kepler_orbit_3body(e, a1, a2, m1 , m2, m3);
Sosy.coord_transf();


double eta = 0.001;
//constant for calculation of adaptive time step, could be dependent on Period P

//fixed time_step or adaptive?
bool adaptive = true;
std::string integrator = "RK4";
//RK4 adaptive is not yet implemented
//pay close attention to the spelling of the string
// open file .txt to save data
std::ofstream ofile;
//ofile.open("test3.txt");
std::string file_dir ="D:/drive dateien/Uni/Semester_7_WiSe23_24/Astro-F_Praktikum/Chaos/data/";
//change directory to your liking
//This line does not give an error if your adress is non-existent
std::string file_name = "3body" + integrator + "_delta="+ scientific_format(delta, 1, 1);
if (adaptive== true){
     file_name = file_name + "_adaptive_eta=" +scientific_format(eta, 1, 1);
}
else{
    file_name = file_name + "_" + int_to_str(step_size_rel_to_orbit)+ "_steps_per_Orbit";
}
std::string file_end = ".txt";
std::string file = file_dir + file_name +file_end;

ofile.open(file);
//initialize time t to zero
double t=0.;
// save to file data: t #p1 x1 y1 z1 vx1 vy1 vz1 #p2 x2 y2 z2 vx2 vy2 vz_2 #p3 x3 y3 z3 vx3 vy3 vz3 p3 e1 e2 a1 a2 dt
//a1 is somehow always 0, a2 is not
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
    ofile  << " "<< scientific_format(Sosy.compute_eccentricity_3body(1), width, prec);  //for m2
    ofile  << " "<< scientific_format(Sosy.compute_eccentricity_3body(2), width, prec);  //m3
    //ofile  << " "<< scientific_format(Sosy.compute_energy(), width, prec);
    ofile  << " "<< scientific_format(Sosy.compute_semi_maj_ax_3body(1), width, prec);  //m2
    ofile  << " "<< scientific_format(Sosy.compute_semi_maj_ax_3body(2), width, prec); //m3
    //ofile  << " "<< scientific_format(norm(Sosy.compute_spec_ang_mom()), width, prec);
    ofile  << " "<< scientific_format(h, width, prec); //added timestep to the text file
    ofile  << std::endl;

     //calculate variable time_step if appropriate
     if(t!=0){
          if (adaptive== true){
              double new_h= Sosy.adaptive_time_step_diff_quot(eta, h);
              h = new_h;          
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
          if(t+h > t_max){
               Sosy.evolveLeapFrog(h, iter-1, iter);       
          }
          else{
               Sosy.evolveLeapFrog(h, i, iter);
          }
        //ignore last step? gets pretty complicated
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
