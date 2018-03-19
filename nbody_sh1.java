import java.io.*;
import java.util.Scanner;
/*
=============================================================================
*
* nbody_sh1.C: an N-body integrator with a shared but variable time step
* (the same for all particles but its size changing in time),
* using the Hermite integration scheme (4th order integrator).
*
* ref.: Hut, P., Makino, J. & McMillan, S., 1995,
* Astrophysical Journal Letters 443, L93-L96.
*
* note: in this first version, all functions are included in one file,
* without any use of a special library or header files.
*_____________________________________________________________________________
*
* USAGE: This is a JAVA program run from the command line. 
* The minimum required command line is:
* java nbody_sh1 input.dat
* 
* Additional options can be specified in the following order:
* 
* dt_param dt_tot dt_out dt_dia 
*
* usage: nbody_sh1 [dt_param: step_size_control_parameter]
* [dt_out: output_interval] [dt_tot: total_duration]
* [dt_dia: time interval between diagnostics]
*
* "step_size_control_parameter" is a coefficient determining the
* the size of the shared but variable time step for all particles
*
* "diagnostics_interval" is the time between output of diagnostics,
* in the form of kinetic, potential, and total energy; with the
* -x option, a dump of the internal particle data is made as well
*
* "output_interval" is the time between successive snapshot outputs
*
* "total_duration" is the integration time, until the program stops
*
* Input/output are read/written from/to the standard i/o streams.
* Since all options have sensible defaults, the simplest way to run
* the code is by only specifying the i/o files for the N-body
* snapshots:
*
* java nbody_sh1 data_in.dat > data_out.dat
*
* The diagnostic information is written to file error.txt
*
* Note: if any of the times specified in the -e, -o, or -t options are not an
* an integer multiple of "step", output will occur slightly later than
* predicted, after a full time step has been taken. And even if they
* are integer multiples, round-off error may induce one extra step.
*_____________________________________________________________________________
*
* External data format:
*
* The program expects input of a single snapshot of an N-body system,
* in the following format: the number of particles in the snapshot n;
* the time t; mass mi, position ri and velocity vi for each particle i,
* with position and velocity given through their three Cartesian
* coordinates, divided over separate lines as follows:
*
* n
* t
* m1 r1_x r1_y r1_z v1_x v1_y v1_z
* m2 r2_x r2_y r2_z v2_x v2_y v2_z
* ...
* mn rn_x rn_y rn_z vn_x vn_y vn_z
*
* Output of each snapshot is written according to the same format.
*
* Internal data format:
*
* The data for an N-body system is stored internally as a 1-dimensional
* array for the masses, and 2-dimensional arrays for the positions,
* velocities, accelerations and jerks of all particles.
* 
* In order for the internal data to be passed to and from a method whilst
* being changed (without resorting to global variables), 
* it is necessary to pass variables by reference which meant that t, epot and
* coll_time all had to be stored as single value arrays.
*_____________________________________________________________________________
*
* version 1: June 2016 David Reid (from v1 of the Art of Computational Physics by Piet Hut and Jun Makino.
*_____________________________________________________________________________
* Earlier versions of the nbody code include:
 *hermite6 uses a fourth order hermite integrator to find the new positions after each timestep.
 *
 *This code takes data from an INPUT FILE rather than having to change the initial conditions in the program
 *itself.
 *v6 changes from calculating collision time based upon velocity (since it could be 0) to a method
 * based on both velocity and acceleration (which is never 0 for an N>1 body simulation).
 * v5 includes a variable time step for more accurate integration during close encounters 
 * v4 IS MORE OBJECT ORIENTED, BREAKING THE CODE INTO MODULAR FUNCTIONS
 *COMMAND LINE ARGUMENT arg[0] is the .dat file .
 *Command line usage:
 *javac hermite6.java
 *
 */

public class nbody_sh1 {
    
    public static final int NDIM = 3;   //number of dimensions used
    
     
    
    public static void main(String[] args) throws IOException{
        
    System.setErr(new PrintStream("diagnostic.txt"));    //diagnostic output file    
    //COMMAND LINE ARGUMENTS FOR dt and t_end
    double dt_param = 0.03; //default value for dt_param
    double dt_dia = 10;      //time interval between diagnostic output
    double dt_out = 0.1;     //time interval between snapshot outputs
    double dt_tot = 10;    //default duration of integration
    boolean init_out = true;   //if true: snapshot output will start at t=0
                                // with an echo of the input snapshot
    boolean x_flag = false;     //if true: extra debugging diagnostics output
    
    //read options, must enter all earlier options to access later ones, can change later
    if (args.length > 1) {
       dt_param = Double.parseDouble(args[1]);  
       if (args.length > 2) {
           dt_tot = Double.parseDouble(args[2]);
        }
       if (args.length > 3) {
           dt_out = Double.parseDouble(args[3]);
        }
        if (args.length > 4) {
           dt_dia = Double.parseDouble(args[4]);
        }
    }
    
    
    //first 2 lines of input file are number of particles and start time t
    int n;     
    double t[] = new double[1];     
         
    // READ INPUT FROM DATA FILE
    if (args.length < 1) {
        System.out.println("Error: you must enter the datafile");
        System.out.println("Usage: java nbody_sh1 input.dat dt_param dt_tot dt_out dt_dia > ouput.dat");
        System.exit(1);
    }
    Scanner reader = new Scanner(new FileInputStream(args[0]));
    n = reader.nextInt();   //read n and t
    t[0] = Double.parseDouble(reader.next());
   
    double pos[][] = new double[n][NDIM]; //initialize vectors
    double vel[][] = new double[n][NDIM];
    double mass[] = new double[n];
    
  //set initial conditions from input file
   get_snapshot(mass,pos,vel,n,reader); 
  
   evolve(mass,pos,vel,n,t,dt_param,dt_dia,dt_out,dt_tot,init_out,x_flag);
  
}
    


 
 /*
 * write_diagnostics -- writes out diagnostics to the error stream:
 *                      if x_flag is true all particle positions, velocities,
 *                      mass, accelerations and jerks are also output.
 * Note Kinetic energy is calculated here whilst potential is done in the 
 * get_acc_jerk_pot_coll() function.
 * 
 * einit needs to be passed by reference as it will be changed and passed back to
 * the program. Therefore it is held as a single value array
 */
 public static void write_diagnostics(double[] mass, double[][] pos, double[][] vel,
                                        double[][] acc, double[][] jerk, int n, double[] t,
                                        double[] epot, int nsteps, double[] einit, boolean init_flag,
                                        boolean x_flag) {
    double ekin = 0;
    
    for (int i=0; i<n;i++) {
        for (int k=0;k<NDIM;k++) {
        ekin += 0.5*mass[i]*vel[i][k]*vel[i][k];
    }
}
     double etot = ekin + epot[0]; //total energy
     
     if(init_flag) {
         einit[0] = etot;
        }
  
    //calculate energy errors
    double e_errAbs = etot - einit[0];
    double e_errRel = e_errAbs/einit[0];
    //print out final energy and errors
    
    System.err.println("At time t = " + t[0] + " after " + nsteps + "steps"); 
    System.err.println("E_kin = " + ekin);
    System.err.println("E_pot = " + epot[0]);
    System.err.println("E_tot = " + etot);
    System.err.println("Absolute energy error = " + e_errAbs);
    System.err.println("Relative error = " + e_errRel);      
    
    if(x_flag) {
        System.err.println("For debugging, here is the internal data");
        for(int i = 0; i<n;i++) {
            int particleNum = i+1;
            System.err.println("Internal data for particle " + particleNum);
            System.err.print(mass[i]);
            for(int k=0;k<NDIM;k++) {
                System.err.print(" ");
                System.err.println(pos[i][k]);
            }
            for(int k=0;k<NDIM;k++) {
                System.err.print(" ");
                System.err.println(vel[i][k]);
        }
            for(int k=0;k<NDIM;k++) {
                System.err.print(" ");
                System.err.println(acc[i][k]);
    }
            for(int k=0;k<NDIM;k++) {
                System.err.print(" ");
                System.err.println(jerk[i][k]);
    }
}
}
}
 
/*
 * get_snapshot(m,r,v,n,t);
 * number of particles and time are read in during main method
 * This method reads in the particle data.
 * 
 */    
   public static void get_snapshot(double[] mass, double[][] pos, double[][] vel, int n, Scanner reader) {
       for(int i=0; i<n; i++) {
         mass[i] = Double.parseDouble(reader.next());
         for(int k= 0; k<NDIM; k++) {
             pos[i][k] = Double.parseDouble(reader.next());
            }
         for(int k= 0; k<NDIM; k++){ //mass, 3 position coords, 3 velocity coords
              vel[i][k] = Double.parseDouble(reader.next());
            }
        }
    }
    
/*
 * put_snapshot(double[] mass, double [][] pos, double[][] vel, int n, double t)
 * For easier analysis with Octave, comment out the n and t output
 */   
public static void put_snapshot(double[] mass, double[][] pos, double[][] vel, int n, double[] t) {
    System.out.println(n);  //number of particles     //comment out n and t if wanting to run the output .dat data in basic octave plotting code.
    System.out.println(t[0]); //current time
    for(int i = 0; i < n; i++) {
        System.out.print(mass[i]); //print mass first
        for(int k=0; k<NDIM;k++) {
            System.out.print(" ");
            System.out.print(pos[i][k]);    //print 3 pos coords
    }
        for(int k=0;k<NDIM;k++) {
            System.out.print(" ");
            System.out.print(vel[i][k]); //print 3 vel coords
        }
        System.out.println("");
    }
        

}

/*-----------------------------------------------------------------------------
* evolve -- integrates an N-body system, for a total duration dt_tot.
* Snapshots are sent to the standard output stream once every
* time interval dt_out. Diagnostics are sent to the standard
* error stream once every time interval dt_dia.
*
* note: the integration time step, shared by all particles at any given time,
* is variable. Before each integration step we use coll_time (short
* for collision time, an estimate of the time scale for any significant
* change in configuration to happen), multiplying it by dt_param (the
* accuracy parameter governing the size of dt in units of coll_time),
* to obtain the new time step size.
*
* Before moving any particles, we start with an initial diagnostics output
* and snapshot output if desired. In order to write the diagnostics, we
* first have to calculate the potential energy, with get_acc_jerk_pot_coll().
* That function also calculates accelerations, jerks, and an estimate for the
* collision time scale, all of which are needed before we can enter the main
* integration loop below.
* In the main loop, we take as many integration time steps as needed to
* reach the next output time, do the output required, and continue taking
* integration steps and invoking output this way until the final time is
* reached, which triggers a `break' to jump out of the infinite loop set up
* with `while(true)'.
* 
* t needs to be passed back to the main method!!
*-----------------------------------------------------------------------------
*/

public static void evolve(double[] mass, double[][] pos, double[][] vel,
                            int n, double[] t, double dt_param, double dt_dia,
                            double dt_out, double dt_tot, boolean init_out,
                            boolean x_flag) {
     double t_dia = t[0] + dt_dia;    //next time for diagnostic output
      double t_out = t[0] + dt_out;    //next time for snapshot output
      double t_end = t[0] + dt_tot;    //final time to finish integration
      
     System.err.println("Starting a Hermite integration for a " + n +
                            "-body system from time t = " + t[0] + " with time step control parameter dt_param = "
                            + dt_param + " until time " + t_end + 
                            " with diagnostic output interval dt_dia = " + dt_dia + 
                            " and snapshot output interval dt_out = " + dt_out + ".");
      
      double acc[][] = new double[n][NDIM];     //accelerations
      double jerk[][] = new double[n][NDIM];    // and jerks for particles
      double epot[] = new double[1];                              //potential energy
      double coll_time[] = new double[1];                         //collision time scale
      
      get_acc_jerk_pot_coll(mass,pos,vel,acc,jerk,n,epot,coll_time);
      
      int nsteps = 0;       //number of integration timesteps completed
      double einit[] = new double[1];         //initial total energy of system
      
      write_diagnostics(mass,pos,vel,acc,jerk,n,t,epot,nsteps,einit,true,x_flag);
      
      if(init_out) {
          put_snapshot(mass,pos,vel,n,t);
        }
        
      while(true) {
          while(t[0] < t_dia && t[0] < t_out && t[0] < t_end) {
              double dt = dt_param*coll_time[0];
              evolve_step(mass,pos,vel,acc,jerk,n,t,dt,epot,coll_time);
              nsteps++;
            }
          if(t[0] >= t_dia) {
              write_diagnostics(mass,pos,vel,acc,jerk,n,t,epot,nsteps,einit,false,x_flag);
              t_dia += dt_dia;
            }
          if(t[0] >= t_out) {
              put_snapshot(mass,pos,vel,n,t);
              t_out += dt_out;
            }
          if(t[0] >= t_end) {
              break;
            }
        }
}

/*-----------------------------------------------------------------------------
* evolve_step -- takes one integration step for an N-body system, using the
* Hermite algorithm.
*------------------------------------------------------------------------
*time t, epot and coll_time need to be passed by reference
*/
public static void evolve_step(double[] mass, double[][] pos, double[][] vel,
                                double[][] acc, double[][] jerk, int n, double[] t,
                                double dt, double[] epot, double[] coll_time) {
                
    //for new Hermite integrator we need to have r,v,a and j at time step i(old) and i+1 
    double old_pos[][] = new double[n][NDIM];
    double old_vel[][] = new double[n][NDIM];
    double old_acc[][] = new double[n][NDIM];
    double old_jerk[][] = new double[n][NDIM];
    
    for (int i = 0; i<n; i++) {
          for(int k = 0; k<NDIM;k++) {
              old_pos[i][k] = pos[i][k];    //assign previous value to old
              old_vel[i][k] = vel[i][k];
              old_acc[i][k] = acc[i][k];
              old_jerk[i][k] = jerk[i][k];
              
            }
        }
        predict_step(pos,vel,acc,jerk,n,dt);
        get_acc_jerk_pot_coll(mass,pos,vel,acc,jerk,n,epot,coll_time);
        correct_step(pos,vel,acc,jerk,old_pos,old_vel,old_acc,old_jerk,n,dt);
        
        t[0] += dt;
}

/*-----------------------------------------------------------------------------
* predict_step -- takes the first approximation of one Hermite integration
* step, advancing the positions and velocities through a
* Taylor series development up to the order of the jerks.
*-----------------------------------------------------------------------------
*/
public static void predict_step(double[][] pos, double[][] vel, double[][] acc,
                                double[][] jerk, int n, double dt) {
  for (int i = 0; i < n; i++) {
      for(int k=0;k<NDIM;k++) {
   pos[i][k] += vel[i][k]*dt + acc[i][k]*dt*dt/2 + jerk[i][k]*dt*dt*dt/6;
   vel[i][k] += acc[i][k]*dt + jerk[i][k]*dt*dt/2; //calculates the pos and velocities at the next step
              //the hermite integrator requires knowledge of now and 1 step forward
  }
}
}

/*-----------------------------------------------------------------------------
* correct_step -- takes one iteration to improve the new values of position
* and velocities, effectively by using a higher-order
* Taylor series constructed from the terms up to jerk at
* the beginning and the end of the time step.
*-----------------------------------------------------------------------------
*/
public static void correct_step(double[][] pos, double[][] vel, double[][] acc,
                                double[][] jerk, double[][] old_pos, double[][] old_vel,
                                double[][] old_acc, double[][] old_jerk, int n, double dt) {
     for(int i = 0; i < n; i++) {
     for(int k = 0; k < NDIM; k++) { //v and r are calculated this way round since v is needed in the r equation
         vel[i][k] = old_vel[i][k] + (old_acc[i][k] + acc[i][k])*dt/2 + (old_jerk[i][k] - jerk[i][k])*dt*dt/12;
         pos[i][k] = old_pos[i][k] + (old_vel[i][k] + vel[i][k])*dt/2 + (old_acc[i][k] - acc[i][k])*dt*dt/12;
        }
    }                               
 }
/*-----------------------------------------------------------------------------
* get_acc_jerk_pot_coll -- calculates accelerations and jerks, and as side
* effects also calculates potential energy and
* * the time scale coll_time for significant changes
* in local configurations to occur.
* __ __
* | --> --> |
* M M | r . v |
* --> j --> --> j | --> ji ji --> |
* a == -------- r ; j == -------- | v - 3 --------- r |
* ji |--> |3 ji ji |--> |3 | ji |--> |2 ji |
* | r | | r | | | r | |
* | ji | | ji | |__ | ji | __|
*
* note: it would be cleaner to calculate potential energy and collision time
* in a separate function. However, the current function is by far the
* most time consuming part of the whole program, with a double loop
* over all particles that is executed every time step. Splitting off
* some of the work to another function would significantly increase
* the total computer time (by an amount close to a factor two).
*
* We determine the values of all four quantities of interest by walking
* through the system in a double {i,j} loop. The first three, acceleration,
* jerk, and potential energy, are calculated by adding successive terms;
* the last, the estimate for the collision time, is found by determining the
* minimum value over all particle pairs and over the two choices of collision
* time, position/velocity and sqrt(position/acceleration), where position and
* velocity indicate their relative values between the two particles, while
* acceleration indicates their pairwise acceleration. At the start, the
* first three quantities are set to zero, to prepare for accumulation, while
* the last one is set to a very large number, to prepare for minimization.
* The integration loops only over half of the pairs, with j > i, since
* the contributions to the acceleration and jerk of particle j on particle i
* is the same as those of particle i on particle j, apart from a minus sign
* and a different mass factor.
*-----------------------------------------------------------------------------
*/
public static void get_acc_jerk_pot_coll(double mass[], double pos[][],
                                     double vel[][], double acc[][],double jerk[][],
                                     int n, double[] epot, double[] coll_time){ 
                                         
       // a method for calculating the acceleration and jerk of each body
       //set initial accelerations and jerk!
    for (int i = 0; i < n; i++) {
        for(int k = 0; k < NDIM; k++) {
            acc[i][k] = 0.0;          //initialise all accelerationto 0
            jerk[i][k] = 0.0;   //initialise the jerk
        }
  }
   //coll_time_sq[] is a single value array so that it is passed by reference                                  
   epot[0] = 0;    
   double coll_time_q = 1e300;  //collision time to 4th power
   double coll_est_q;           //collision time estimate to 4th power
  
  for(int i=0;i<n;i++) {
      for(int j = i+1; j < n; j++) {
          double rji[] = new double[NDIM]; //distance between body i and j
          double vji[] = new double[NDIM];
          for(int k = 0; k < NDIM; k++) {
              rji[k] = pos[j][k] - pos[i][k];    //x,y,z components of distance
              vji[k] = vel[j][k] - vel[i][k];
            }
        double r2 = 0.0;
        double v2 = 0.0;
        double rv_r2 = 0;   // (rij . vij) / |rij|^2
        for(int k = 0; k<NDIM;k++){
            r2 += rji[k]*rji[k];
            v2 += vji[k]*vji[k];
            rv_r2 += rji[k]*vji[k];
        }
        rv_r2 /= r2;
        double r = Math.sqrt(r2);
        double r3 = r2*r;
       
      // add the {i,j} contribution to the total pot energy
      
        epot[0] -= mass[i]*mass[j]/r;
        
      // add the {j,i} contribution to the {i,j} values of acc and jerk:
        
        double da[] = new double[NDIM]; //main terms in pairwise
        double dj[] = new double[NDIM]; // acceleration and jerk
        for (int k = 0; k<NDIM; k++) {
            da[k] = rji[k]/r3;
            dj[k] = (vji[k] - 3*rv_r2*rji[k])/r3;
        }
        for(int k = 0; k<NDIM;k++){
            acc[i][k] += mass[j]*da[k]; //acceleration of particle i due to j
            acc[j][k] -= mass[i]*da[k]; // acceleration of j due to i is just opposite 
            jerk[i][k] += mass[j]*dj[k];
            jerk[j][k] -= mass[i]*dj[k]; 
        }
      //first collision time estimate based on linear, unaccelerated motion  
        coll_est_q = (r2*r2)/(v2*v2); //to give t^4 units!
        if (coll_time_q > coll_est_q) {
            coll_time_q = coll_est_q;
        }
      //second collision time estimate based on free fall (acceleration):  
        double da2 = 0;
        for (int k=0; k<NDIM; k++) {
            da2 += da[k]*da[k]; //find the relative acceleration of both bodies
        }
        double mij = mass[i] + mass[j];
        da2 *= mij*mij;
        coll_est_q = r2/da2;    //units are t^4
        if(coll_time_q > coll_est_q) {
            coll_time_q = coll_est_q;
        }
    }
}          
    coll_time[0] = Math.sqrt(Math.sqrt(coll_time_q));                                    
 }

 /*-----------------------------------------------------------------------------
*                                           \\ o
* end of file: nbody_sh1.C                  /\\' O
*                                          /\ |
*=============================================================================
*/
}