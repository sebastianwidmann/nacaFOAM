# nacaFOAM

CFD simulation pipeline using the `OpenFOAM` framework to generate a 2D database for machine learning purposes using
NACA symmetrical, 4- and 5-digit airfoils.

Machine learning database consists of subsonic simulations (Mach < 0.8) with varying angles of
attack `alpha = [-20, 20] deg`.

## Dependencies

* OpenFOAM v2112 or newer
* Python 3.10 or newer
* NumPy 1.22.4 or newer
* matplotlib 3.5.2 or newer

## Running
Before running the database, execute `source /lib/openfoam/openfoam2206/etc/bashrc` in the terminal