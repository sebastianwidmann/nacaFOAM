# nacaFOAM

CFD simulation pipeline using the `OpenFOAM` framework to generate a 2D 
database for machine learning purposes using NACA symmetrical, 4- and 
5-digit airfoils.

The pipeline can generate results of subsonic simulations for Mach numbers 
`M < 0.6` and any angle of attack. However, it is advised to stay within 
`alpha = [-20, 20] deg` as any lower or higher angle will only predict 
separation.

The detailed explanation about the simulation pipeline as well as the 
performed verification and validation can be found here (LINK TO PAPER).

To reference this algorithm, please cite the repository as:
```
@misc{widmann2023nacaFOAM,
  title={Automated Aerodynamic 2D Data Generation for Machine Learning},
  author={Widmann, Sebastian and Schlichter, Philipp and Reck, Michaela and Indinger, Thomas},
  year={2023},
  publisher={GitHub},
  howpublished={\url{https://github.com/sebastianwidmann/nacaFOAM}},
}
```

## 1. AirfoilMNIST: A Large-Scale Dataset based on Two-Dimensional RANS Simulations of Airfoils

The subsets of the *airfoilMNIST* datasets are subsequently linked here once the data has been published.
* airfoilMNIST-raw: Not uploaded yet
* airfoilMNIST: Not uploaded yet
* airfoilMNIST-incompressible: https://syncandshare.lrz.de/getlink/fiHGg7CzJ1tajtpDxwR65D/

## 2. Running the nacaFOAM pipeline

### 2.1 Dependencies

* OpenFOAM v2206 or newer
* Python packages as set in *requirements.txt*

### 2.2 Set OpenFOAM environment

* To execute the main script `nacaFOAM.py`, one has to enter the OpenFOAM 
  environment or source toward the bashrc script
    * **Ubuntu**: `source /lib/openfoam/openfoam2206/etc/bashrc`

### 2.3 Set Python environment on HPC systems
* Typically, Python is installed in the `usr/bin/python3.XX` directory on 
  local Linux machines. However, this is often not true for HPC systems.
  * By default, `nacaFOAM.py` will look for the executable of Python in the 
    following directory -> `usr/bin/python3`
  * If a virtual environment is used, the directory to the Python executable 
    must be specified in the `shebang line` at the beginning of `nacaFOAM.py`

### 2.4 Running parameters & core count settings

* To execute the main script, 3 parameters must be provided in the form of a 
  parameter list (see line 335 in `nacaFOAM.py`)
    * **Angles of attack**
    * **Mach numbers**
    * **Airfoil codes**
* The simulations will be saved into the `database` directory which will be 
  generated once the main script is executed.
  To minimise the storage capacity, it is advised to delete the raw data 
  folders and only keep the `VTK` files in the
  database folder.

**Optional** If `rhoSimpleFOAM` should run on a different core count than 4, 
respecify the following arguments

* `nacaFOAM.py`: the argument `nProc` in line 336,
* `nacaFOAM.py`: the argument `nr` in line 262,
* `nacaFOAM-template/system/decomposeParDict`: modify the argument 
  `numberOfSubdomains` in line 17

### 2.5 Run multiple simulations

* Specify range of angles of attack as `list` or `array` of **strings** in 
  line 325
    *     angles = np.arange(-5, 21, 1)
* Specify range of Mach numbers as `list` or `array` in line 326
    *     angles = np.arange(0.05, 0.65, 0.05)
* Specify range of NACA airfoils through function `generateNacaAirfoils()`
    *     nacas = generateNacaAirfoils()
        * In the function the different digits can be modified
            * NACA 4-series --> see lines 42-44
            * NACA 5-series --> see lines 53-56
* Execute `nacaFOAM.py` script

### 2.6 Run single simulation

* Easiest way to run a single simulation is by replacing lines 325-327 
  with the following lines
    *     angles = np.array([0])
    *     machs = np.array([0.15])
    *     nacas = ['0012']
* Execute `nacaFOAM.py` script

### 2.7 Postprocessing

* To generate the textfile with the lift, drag and moment coefficients for 
  all simulations, simply execute the `postProcessing.py` script after the 
  `nacaFOAM.py` script has been executed sucessfully. This script will generate
  a new `.csv` file within the `database` directory.



