# ABP2023
Sedimention of active Brownian particles in a box

Codes used in the scientific publications (unpublished).

<b>Interacting ABPs</b>
C++ code on the sedimention of interacting active Brownian particles in a box.
Compile: g++ interactingABPs.cpp -lgsl -lgslcblas -lm -O3 -s -o interactingABPs.out.
Run: ./interactingABPs.out -parameter=value.
List of parameters: Pe, alpha, F0, LX, LY, init, RAN, dt, teq, tmax, Npart (details as comments in the code).

<b>Ideal ABPs</b>
FreeFem++ code (https://freefem.org/).
Run: FreeFem++ -v 0 -nw idealABPs.edp -parameter value</br>
List of parameters: vp, vg, LX, LY, Nx, Ny, Nt (details as comments in the code).
