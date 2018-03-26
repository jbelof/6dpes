# 6dpes

For molecules whose symmetry is described by 3 translational, 2 rotational and 1 vibrational mode, this code produces the potential energy surface as a function of those internal variables.  The 6D PES resulting from this calculation can include an external influence, such as a nano surface, external electric field or mean field description of neighbors.  With the 6D PES in hand, additional calculations (such as coupled quantum translations/rotations/vibrations can be applied).

The main use of this code is to calculate vibrational spectroscopy, by direct diagonalization of the QM Hamiltonian matrix, for effects such as the vibrational Stark effect (vibrational line splitting).

For use of this 6dpes code, please cite I. Matanovic, J.L. Belof, B. Space, K. Sillar, J. Sauer, J. Eckert and B. Zlatko, "Hydrogen adsorbed in a metal organic framework-5: Coupled translation-rotation eigenstates from quantum five-dimensional calculations", J. Chem. Phys. 137:014701 (2012),  http://doi.org/10.1063/1.4730906

For use of the hydrogen potential used in an example, please cite J.L. Belof, A.C. Stern and B. Space, "An accurate and transferable intermolecular diatomic hydrogen potential for condensed phase simulation", J. Chem. Theory. Comput., 4:1332 (2008), http://doi.org/10.1021/ct800155q

For use of the nanomaterial potential (MOF-5) used in an example, please cite J.L. Belof, A.C. Abraham and B. Space, "A predictive model of hydrogen sorption for metal-organic materials", J. Phys. Chem. C., 113:9316 (2009), http://doi.org/10.1021/jp901988e


## Getting Started

After obtaining the source code, please consult the Makefile to set any specific compiler flags (defaults are gcc with std C library).

The code implements a potential surface that described long-range electrostatics, electronic repulsion and dispersion forces and many-body polarization.


## Installing

Compilation is simple and relies on only standard libraries:  

$ make  
gcc -c -O3 -DPOLAR_SELF -I. main.c  
gcc -c -O3 -DPOLAR_SELF -I. cleanup.c  
gcc -c -O3 -DPOLAR_SELF -I. input.c  
gcc -c -O3 -DPOLAR_SELF -I. pairs.c  
gcc -c -O3 -DPOLAR_SELF -I. pbc.c  
gcc -c -O3 -DPOLAR_SELF -I. surface.c  
gcc -c -O3 -DPOLAR_SELF -I. energy.c  
gcc -c -O3 -DPOLAR_SELF -I. polar.c  
gcc -O3 -DPOLAR_SELF *.o -o 6dpes  
$  


## Running the examples

Run the 6dpes binary without arguments to obtain the command line input.  In addition, examine some of the scripts included for examples of vibrational line splitting in hydrogen due to either nanosurface interface or an imposed external electric field.

$ ./6dpes   
usage: ./6dpes <PDB filename> [<b1_x> <b1_y> <b1_z> <b2_x> <b2_y> <b2_z> <b3_x> <b3_y> <b3_z>] [<xi> <yi> <zi> <ri> <xf> <yf> <zf> <rf>] [<dx> <dy> <dz> <dtheta> <dphi> <dr>] [ efield_x efield_y efield_z ]  
	[b1,...,b3] : the cartesian basis vectors (in A) of the unit cell, each b vector is a row of the basis matrix  
	[xi,...,zf,rf] : the initial and final xyz c.o.m. and r coordinates for the PES generation  
	[dx,...,dphi,dr] : the step size for c.o.m. coordinates, r, and for the spherical polar angles (theta angle of rotation around the Y axis, phi is the angle around the Z axis)  
	[efield_x,...efield_z] : apply an electric field to the system, in units of e/A^2  
	columnar output is [x,y,z,theta,phi,r,E,BOND,LJ,ES,POL)]  


MOF5 example:  
	$ ./6dpes MOF5+BSS.pdb 25.669 0.0 0.0 0.0 25.669 0.0 0.0 0.0 25.669 -10.0 -10.0 -10.0 0.2 10.0 10.0 10.0 2.0 0.001 0.001 0.001 0.1 0.1 0.01 0.0 0.0 0.0  
would map the PES at each point within a subcube of the unit cell spanned by [-10,-10,-10]x[10,10,10] for all angles theta=0-pi,phi=0-2*pi and bond distance 0.5 to 1.0 A in 0.1 A increments  


for a vibrational surface only:  
	$ ./6dpes MOF5+BSS.pdb 25.669 0.0 0.0 0.0 25.669 0.0 0.0 0.0 25.669 -10.0 -10.0 -10.0 0.2 -10.0 -10.0 -10.0 4.0 1000.0 1000.0 1000.0 1000.0 1000.0 0.001 0.0 0.0 0.0 | awk '{print $6 " " $7}'  


for the vibrational surface of an H2 in a 0.5 e/A^2 field transverse to it's axis:  


	$ ./6dpes BSSP.pdb 1000.0 0.0 0.0 0.0 1000.0 0.0 0.0 0.0 1000.0 -10.0 -10.0 -10.0 0.2 -10.0 -10.0 -10.0 4.0 1000.0 1000.0 1000.0 1000.0 1000.0 0.001 0.0 0.0 0.5 | awk '{print $6 " " $7}'  
or to evaluate the vib frequency (cm^-1) by second derivative (central difference) at the minimum:  


	$ ./6dpes BSSP.pdb 1000.0 0.0 0.0 0.0 1000.0 0.0 0.0 0.0 1000.0 -10.0 -10.0 -10.0 0.7419 -10.0 -10.0 -10.0 0.7421 1000.0 1000.0 1000.0 1000.0 1000.0 0.0001 0.00003 0.0 0.0 | awk '{print $7}' | column | awk '{}{printf("%.16f\n",sqrt(($3+$1-2.0*$2)/1.0e-8)*42.796);}{}'  


or to evaluate the vib frequency (cm^-1) by second derivative at global min in MOF5:  
	$ ./6dpes MOF5+H2.globalmin.pdb 25.669 0.0 0.0 0.0 25.669 0.0 0.0 0.0 25.669 -4.176 -4.216 8.437 0.7419 -4.176 -4.216 8.437 0.7421 1000.0 1000.0 1000.0 1000.0 1000.0 0.0001 0.0 0.0 0.0 | awk '{print $7}' | column | awk '{}{printf("%.16f\n",sqrt(($3+$1-2.0*$2)/1.0e-8)*42.796);}{}'  



## Authors

* **Jon Belof** [jbelof@github](https://github.com/jbelof)  

[google scholar](https://scholar.google.com/citations?user=gNrlNbwAAAAJ&hl=en)  
[research gate](https://www.researchgate.net/profile/Jon_Belof)  
[linkedin](http://www.linkedin.com/in/jbelof)  
[web profile](http://jbelof.academia.edu)  


## License

This project is licensed under the GNU General Public License v3, please see GPL_license.txt for details.


