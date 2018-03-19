import java.util.Random;
import java.util.Arrays;
import java.io.*;
/**
 * Creates a snapshot for input to the nbody_sh1 N-body simulation in the format:
 * N
 * t
 * m1 r1_X r1_y r1_z v_x v_y v_z
 * m2........
 * 
 * This code generates the initial conditions according to the Plummer model, given by Aarseth et al
 * 
 * v2 Includes code for quiet start - ensuring proper distribution of stars
 *      Also included is a new method to bring the centre of mass to 0 and the velocity
 *      of the centre of mass = 0
 * 
 * We will use standard units r_virial = 1 and a = 3PI/16 - without the scale factor a = 1.
 * This code is written from the Ruby version given in Vol.11 of the ACS N body lab tutorials.
 * A comparison of Numerical Methods for the Study of
 * Star Cluster Dynamics, by Sverre Aarseth, Michel Henon, and Roland Wielen,
 * which appeared ages ago, in 1974, in Astron. Astroph. 37, 183.
 * 
 * The quartiles method checks that the code gives us the expected Plummer model
 */
public class mkPlummer2
{

 public static final int NDIM = 3;

    public static void main(String[] args) throws IOException{
        System.setErr(new PrintStream("diagnosticIC.txt"));
        if(args.length == 1) {
        int n = Integer.parseInt(args[0]);     //number of particles to generate positions for
        double mass[] = new double[n];          //masses of all particles
        double pos[][] = new double[n][NDIM];   //position and velocity
        double vel[][] = new double[n][NDIM];   // arrays
        
        plummer(mass,pos,vel,n);         //generate random positions in a unit sphere
        adjust_CoM(pos,vel,mass,n);     //centres the cluster before outputting snapshot
        put_snapshot(mass,pos,vel,n,0); //output snapshot at time t=0
        //quartiles(pos);                 //TEST PLUMMER MODEL
    } else{
        System.out.println("Error, include the number of bodies in the command line");
    }
        
        
        
        
        
    }

/* plummer - uses a random number generator to add positions of n particles
 * inside a unit sphere
 * 
*/
 public static void plummer(double[] mass, double[][] pos, double[][] vel, int n) {
        double scale_factor = 16/(Math.PI*3);   //gets distances and velocities into the correct, standard units. Without this factor a = 1
        double cumulative_mass_min = 0.0;
        double cumulative_mass_max = 1.0/n;
        for(int i = 0; i < n; i++) {
            mass[i] = 1.0/n;            //n equal mass bodies
            double cumulative_mass = frand(cumulative_mass_min, cumulative_mass_max);       //random mass in layer, 
            cumulative_mass_min = cumulative_mass_max;
            cumulative_mass_max += 1.0/n;
            double r = 1.0/(Math.sqrt(Math.pow(cumulative_mass,(-2.0/3.0)) - 1.0));        //r(m) distribution from Aarseth et al
            double theta = Math.acos(frand(-1,1));    //uniformly distributed angle between 0 and PI
            double phi = frand(0,2*Math.PI);;   //angle between 0 and 2PI
            pos[i][0] = r*Math.sin(theta)*Math.cos(phi)/scale_factor;    //x coord in spherical coords
            pos[i][1] = r*Math.sin(theta)*Math.sin(phi)/scale_factor;    //y coord in spherical coords
            pos[i][2] = r*Math.cos(theta)/scale_factor;                  //z coord in spherical coords
            
            //assigning velocities according to Aarseth et al.
            //rejection technique
            double x = 0.0;
            double y = 0.1;
            while(y > Math.pow(x*x*(1.0-x*x),3.5)) {
                x = frand(0,1);
                y = frand(0,0.1);
            }
            double v = x*Math.sqrt(2.0)*Math.pow((1.0+r*r),-0.25);  //escape velocity of particle/max velocity possible
            theta = Math.acos(frand(-1,1));
            phi = frand(0,2*Math.PI);
            vel[i][0] = v*Math.sin(theta)*Math.cos(phi)*Math.sqrt(scale_factor);
            vel[i][1] = v*Math.sin(theta)*Math.sin(phi)*Math.sqrt(scale_factor);
            vel[i][2] = v*Math.cos(theta)*Math.sqrt(scale_factor);
        }
    }
    
    
    /*
 * put_snapshot(double[] mass, double [][] pos, double[][] vel, int n, double t)
 */   
public static void put_snapshot(double[] mass, double[][] pos, double[][] vel, int n, double t) {
    System.out.println(n);  //number of particles
    System.out.println(t); //current time
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
/*
 * frand generates a random number between a and b
 */
public static double frand(double a, double b) {
    double r = Math.random();
    r = a + r*(b-a);
    return r;
}

/*quartiles - gives the radius within which 1/4,1/2 and 3/4 of the mass is contained
 * 
 * Since each star has the same mass, if we order all the stars by radius
 * the first 1/4 has a 1/4 of the mass etc.
 */
public static void quartiles(double[][] pos) {
    int n = pos.length;
    double r2[] = new double[n];
    for(int i = 0; i < n; i++) {
        double rsq = pos[i][0]*pos[i][0] + pos[i][1]*pos[i][1] + pos[i][2]*pos[i][2];
        r2[i] = rsq;
    }
    //order all square radii and find the radius at 1/4,1/2 and 3/4
    Arrays.sort(r2);
    double r_1 = Math.sqrt(r2[n/4 - 1]);
    double r_2 = Math.sqrt(r2[n/2 - 1]);
    double r_3 = Math.sqrt(r2[3*n/4 - 1]);
    System.err.println("The 3 quartiles for r(M) are:");
    System.err.println("r(1/4) = " + r_1);
    System.err.println("r(1/2) = " + r_2);
    System.err.println("r(3/4) = " + r_3);
}

/*
 * adjust_CoM - takes in the position and velocity vectors and adjusts them
 * so that the centre of mass (CoM) is at the origin with 0 velocity
 */
public static void adjust_CoM(double[][] pos, double[][] vel, double[] mass, int n) {
    double[] vel_com = new double[NDIM];
    double[] pos_com = new double[NDIM];
    for(int i = 0; i < n; i++) {    //finds the centre of mass and CoM velocity
        for(int k = 0; k < NDIM; k++) {
            pos_com[k] += pos[i][k]*mass[i];
            vel_com[k] += vel[i][k]*mass[i];
        }
    }
    for(int i = 0; i < n; i++) {    //subtract the centre of mass and CoM velocity from
        for(int k = 0; k < NDIM; k++) { // each body in order to bring the CoM to origin
            pos[i][k] -= pos_com[k];
            vel[i][k] -= vel_com[k];
        }
    }
}
}
