# Sedimentation of active Brownian particles
<a href="https://dx.doi.org/10.5281/zenodo.8353561"><img src="https://zenodo.org/badge/399397571.svg" alt="DOI"></a>

Codes used in the scientific publication:</br>
M. Mangeat, S. Chakraborty, A. Wysocki, and H. Rieger, <i><a href='https://journals.aps.org/pre/abstract/10.1103/PhysRevE.109.014616'>Stationary particle currents in sedimenting active matter wetting a wall</a></i>, Phys. Rev. E 109, 014616 (2024). Preprint available on <a href='https://arxiv.org/abs/2309.09714'>arXiv</a>.

## Interacting ABPs

C++ code on the sedimentation of interacting active Brownian particles in a box.</br>
Exportations: dynamics with Delta t = 1 between two datafiles, steady-state density, polarization, and current.</br>
Compile: g++ interactingABPs.cpp -lgsl -lgslcblas -lm -O3 -s -o interactingABPs.out.</br>
Run: ./interactingABPs.out -parameter=value.</br>
List of parameters: Pe, alpha, F0, LX, LY, init, RAN, dt, teq, tmax, Npart (details as comments in the code).</br>
Generate the movie: python figure_interactingABPs_dynamics.py -parameter=value (parameters: Pe, alpha, F0, LX, LY, init, ran, NCPU, movie).</br>
Generate the steady-state figure: python figure_interactingABPs_steadystate.py -parameter=value (parameters: Pe, alpha, F0, LX, LY, init, Nran).

## Ideal ABPs

<a href='https://freefem.org/'>FreeFem++ code</a> on the sedimentation of non-interacting active Brownian particles in a box.</br>
Exportations: steady-state density and polarization.</br>
Run: FreeFem++ -v 0 -nw idealABPs.edp -parameter value.</br>
List of parameters: vp, vg, LX, LY, Nx, Ny, Nt (details as comments in the code).</br>
Generate the steady-state figure: python figure_idealABPs_steadystate.py -parameter=value (parameters: vp, vg, LX, LY).
