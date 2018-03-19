## Copyright (C) 2016 dpryd
## 
## v2: saves individual plots ready to encode them as a movie file for better animation.

## this script is intended to animate a sequence of plots from output data from
## nbody_sh1.java
## To make this program simpler the output of n and t from nbody_sh1 needs to be suppressed
## so that only the positions (and velocities) of the bodies are output to the .dat file
## REQUIRED INPUT: n the number of bodies in the simulation
##		
		
data = load ("reverse.dat");
n = 100;	%the number of bodies in the simulation
%dt_out = 0.1;
%dt_tot = 10;
%nout = dt_tot/dt_out;	% the number of data outputs given from the nbody simulation
x = data(:,2);
y = data(:,3);
z = data(:,4);
figure 1
im = 0;
for i=1:n:rows(data)	% gather data at 1 timestep for a plot
	xo = x(i:i+n-1);
	yo = y(i:i+n-1);
	zo = z(i:i+n-1);
	%figure(1,'visible','off')
	%scatter(xo,zo)			%change between 2D and 3D plot
	scatter3(xo,yo,zo);
	%axis([-1 1 -1 1], "manual", "equal")
	%axis ([-20 0 -5 5 -5 5], "manual", "equal");	%set fixed axis size
	xlabel("X");
	ylabel("Y");
	zlabel("Z");
	pause(0.0001);
	filename= strcat("reverse",sprintf('%04d',im),'.png');
    print(filename);
	im+=1;
endfor
	

	
