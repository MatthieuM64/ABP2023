# ABP2023
Sedimention of active Brownian particles (ABPs) in a box

Codes used in the scientific publications (unpublished).

<b>Interacting ABPs</b></br>
C++ code on the sedimention of interacting active Brownian particles in a box.</br>
Exportations: dynamics with Delta t = 1 between two datafiles, steady state density, polarization and current.</br>
Compile: g++ interactingABPs.cpp -lgsl -lgslcblas -lm -O3 -s -o interactingABPs.out.</br>
Run: ./interactingABPs.out -parameter=value.</br>
List of parameters: Pe, alpha, F0, LX, LY, init, RAN, dt, teq, tmax, Npart (details as comments in the code).

<b>Ideal ABPs</b></br>
FreeFem++ code (https://freefem.org/).</br>
Exportations: steady state density and polarization.</br>
Run: FreeFem++ -v 0 -nw idealABPs.edp -parameter value.</br>
List of parameters: vp, vg, LX, LY, Nx, Ny, Nt (details as comments in the code).
