# nbody
A Java implementation of the N-body integrator of Hut, P., Makino, J. &amp; McMillan, S.
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

Code to generate specific types of snapshot has been included in the repository:

writeSphere: generates a spherical distribution of particles
mkPlummer: generates a distribution of particles according to the Plummer model. 

The outputs of these java programs can then be input into nbody_sh1.
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

June 2016 David Reid (from v1 of the Art of Computational Physics by Piet Hut and Jun Makino.
