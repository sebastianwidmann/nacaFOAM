#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------
# Created By  : Sebastian Widmann
# Institution : TU Munich, Department of Aerospace and Geodesy
# Created Date: June 18, 2022
# version ='1.0'
# ---------------------------------------------------------------------------
"""
Class implementation to generate the initial condition files in the zero
directory. Initial values generated for kinematic eddy viscosity nut,
initial pressure p, specific dissipation rate omega, initial temperature T,
turbulent kinetic energy k, turbulent thermal diffusivity alpha_t and initial
velocity u. Initial condition files written into the "/0.orig" directory.
"""
# ---------------------------------------------------------------------------

import numpy as np
import argparse

# Define global constants for air
R = 287.058  # [J*kg^-1*K^-1] Specific gas constant for dry air
gamma = 1.4  # [-] Heat capacity ratio

# Define ICAO standard atmosphere properties
p0 = 101325  # [Pa] Total pressure
T0 = 288.15  # [K] Total temperature


class generateInitialConditions(object):
    def __init__(self, mach):
        self.mach = mach
        self.p = None  # [Pa] Static pressure
        self.T = None  # [K] Static temperature
        self.u = None  # [m * s^-1] Freestream velocity in x-direction
        self.k = None  # [m^2 * s^-2] Turbulent kinetic energy
        self.omega = None  # [s^-1] Turbulence specific dissipation rate

        # Constants
        self.I = 0.05  # [%] Turbulence intensity
        self.C_mu = 0.09  # [-] Constant related to k-omega SST model
        self.L = 1  # [m] Reference length

        self.calculatePressure()
        self.calculateTemperature()
        self.calculateVelocity()
        self.calculateTurbulentKineticEnergy()
        self.calculateSpecificDissipationRate()

        self.writeToFile()

    def writeToFile(self):
        self.writeToFile_KinematicEddyViscosity()
        self.writeToFile_Pressure()
        self.writeToFile_SpecificDissipationRate()
        self.writeToFile_Temperature()
        self.writeToFile_TurbulentKineticEnergy()
        self.writeToFile_TurbulentThermalDiffusivity()
        self.writeToFile_Velocity()

    def calculatePressure(self):
        self.p = p0 / np.power(1 + 0.5 * (gamma - 1) * self.mach ** 2, - gamma / (gamma - 1))

    def calculateTemperature(self):
        self.T = T0 / (1 + 0.5 * (gamma - 1) * self.mach ** 2)

    def calculateVelocity(self):
        a = np.sqrt(gamma * R * self.T)  # [m*s^-1] Speed of sound
        self.u = self.mach * a  # [m*s^-1] Freestream velocity

    def calculateTurbulentKineticEnergy(self):
        self.k = 1.5 * np.power(self.I * np.absolute(self.u), 2)

        print(np.absolute(self.u))

    def calculateSpecificDissipationRate(self):
        self.omega = np.power(self.k, 0.5) / (np.power(self.C_mu, 0.25) * self.L)

    def writeToFile_KinematicEddyViscosity(self):
        f = open('0.orig/nut', 'w+')

        f.write('/*--------------------------------*- C++ -*----------------------------------*\\   \n')
        f.write('| =========                 |                                                 |    \n')
        f.write('| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |  \n')
        f.write('|  \\\\    /   O peration     | Version:  v2112                                 |  \n')
        f.write('|   \\\\  /    A nd           | Website:  www.openfoam.com                      |  \n')
        f.write('|    \\\\/     M anipulation  |                                                 |  \n')
        f.write('\\*---------------------------------------------------------------------------*/   \n')
        f.write('FoamFile                                                                           \n')
        f.write('{                                                                                  \n')
        f.write('    version     2.0;                                                               \n')
        f.write('    format      ascii;                                                             \n')
        f.write('    class       volScalarField;                                                    \n')
        f.write('    object      nut;                                                               \n')
        f.write('}                                                                                  \n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //    \n')
        f.write('                                                                                   \n')
        f.write('dimensions     [0 2 -1 0 0 0 0];                                                   \n')
        f.write('                                                                                   \n')
        f.write('internalField  uniform 0;                                                          \n')
        f.write('                                                                                   \n')
        f.write('boundaryField                                                                      \n')
        f.write('{                                                                                  \n')
        f.write('   freestream                                                                      \n')
        f.write('   {                                                                               \n')
        f.write('       type            calculated;                                                 \n')
        f.write('       value           uniform 0;                                                  \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   wall                                                                            \n')
        f.write('   {                                                                               \n')
        f.write('       type            nutkWallFunction;                                           \n')
        f.write('       value           uniform 0;                                                  \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   #includeEtc "caseDicts/setConstraintTypes"                                      \n')
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('// ************************************************************************* //    \n')
        f.close()

    def writeToFile_Pressure(self):
        f = open('0.orig/p', 'w+')

        f.write('/*--------------------------------*- C++ -*----------------------------------*\\   \n')
        f.write('| =========                 |                                                 |    \n')
        f.write('| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |  \n')
        f.write('|  \\\\    /   O peration     | Version:  v2112                                 |  \n')
        f.write('|   \\\\  /    A nd           | Website:  www.openfoam.com                      |  \n')
        f.write('|    \\\\/     M anipulation  |                                                 |  \n')
        f.write('\\*---------------------------------------------------------------------------*/   \n')
        f.write('FoamFile                                                                           \n')
        f.write('{                                                                                  \n')
        f.write('    version     2.0;                                                               \n')
        f.write('    format      ascii;                                                             \n')
        f.write('    class       volScalarField;                                                    \n')
        f.write('    object      p;                                                                 \n')
        f.write('}                                                                                  \n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //    \n')
        f.write('                                                                                   \n')
        f.write('pOut           {}; \n'.format(self.p))
        f.write('                                                                                   \n')
        f.write('dimensions     [1 -1 -2 0 0 0 0];                                                  \n')
        f.write('                                                                                   \n')
        f.write('internalField  uniform $pOut;                                                      \n')
        f.write('                                                                                   \n')
        f.write('boundaryField                                                                      \n')
        f.write('{                                                                                  \n')
        f.write('   freestream                                                                      \n')
        f.write('   {                                                                               \n')
        f.write('       type            freestreamPressure;                                         \n')
        f.write('       freestreamValue uniform $pOut;                                              \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   wall                                                                            \n')
        f.write('   {                                                                               \n')
        f.write('       type            zeroGradient;                                               \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   #includeEtc "caseDicts/setConstraintTypes"                                      \n')
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('// ************************************************************************* //    \n')
        f.close()

    def writeToFile_SpecificDissipationRate(self):
        f = open('0.orig/omega', 'w+')

        f.write('/*--------------------------------*- C++ -*----------------------------------*\\   \n')
        f.write('| =========                 |                                                 |    \n')
        f.write('| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |  \n')
        f.write('|  \\\\    /   O peration     | Version:  v2112                                 |  \n')
        f.write('|   \\\\  /    A nd           | Website:  www.openfoam.com                      |  \n')
        f.write('|    \\\\/     M anipulation  |                                                 |  \n')
        f.write('\\*---------------------------------------------------------------------------*/   \n')
        f.write('FoamFile                                                                           \n')
        f.write('{                                                                                  \n')
        f.write('    version     2.0;                                                               \n')
        f.write('    format      ascii;                                                             \n')
        f.write('    class       volScalarField;                                                    \n')
        f.write('    object      omega;                                                             \n')
        f.write('}                                                                                  \n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //    \n')
        f.write('                                                                                   \n')
        f.write('omegaInlet     {}; \n'.format(self.omega))
        f.write('                                                                                   \n')
        f.write('dimensions     [0 0 -1 0 0 0 0];                                                   \n')
        f.write('                                                                                   \n')
        f.write('internalField  uniform $omegaInlet;                                                \n')
        f.write('                                                                                   \n')
        f.write('boundaryField                                                                      \n')
        f.write('{                                                                                  \n')
        f.write('   freestream                                                                      \n')
        f.write('   {                                                                               \n')
        f.write('       type            inletOutlet;                                                \n')
        f.write('       inletValue      uniform $omegaInlet;                                        \n')
        f.write('       value           uniform $omegaInlet;                                        \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   wall                                                                            \n')
        f.write('   {                                                                               \n')
        f.write('       type            omegaWallFunction;                                          \n')
        f.write('       value           uniform $omegaInlet;                                        \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   #includeEtc "caseDicts/setConstraintTypes"                                      \n')
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('// ************************************************************************* //    \n')
        f.close()

    def writeToFile_Temperature(self):
        f = open('0.orig/T', 'w+')

        f.write('/*--------------------------------*- C++ -*----------------------------------*\\   \n')
        f.write('| =========                 |                                                 |    \n')
        f.write('| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |  \n')
        f.write('|  \\\\    /   O peration     | Version:  v2112                                 |  \n')
        f.write('|   \\\\  /    A nd           | Website:  www.openfoam.com                      |  \n')
        f.write('|    \\\\/     M anipulation  |                                                 |  \n')
        f.write('\\*---------------------------------------------------------------------------*/   \n')
        f.write('FoamFile                                                                           \n')
        f.write('{                                                                                  \n')
        f.write('    version     2.0;                                                               \n')
        f.write('    format      ascii;                                                             \n')
        f.write('    class       volScalarField;                                                    \n')
        f.write('    object      T;                                                                 \n')
        f.write('}                                                                                  \n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //    \n')
        f.write('                                                                                   \n')
        f.write('Tinlet         {}; \n'.format(self.T))
        f.write('                                                                                   \n')
        f.write('dimensions     [0 0 0 1 0 0 0];                                                    \n')
        f.write('                                                                                   \n')
        f.write('internalField  uniform $Tinlet;                                                    \n')
        f.write('                                                                                   \n')
        f.write('boundaryField                                                                      \n')
        f.write('{                                                                                  \n')
        f.write('   freestream                                                                      \n')
        f.write('   {                                                                               \n')
        f.write('       type            inletOutlet;                                                \n')
        f.write('       inletValue      uniform $Tinlet;                                            \n')
        f.write('       value           $inletValue;                                                \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   wall                                                                            \n')
        f.write('   {                                                                               \n')
        f.write('       type            zeroGradient;                                               \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   #includeEtc "caseDicts/setConstraintTypes"                                      \n')
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('// ************************************************************************* //    \n')
        f.close()

    def writeToFile_TurbulentKineticEnergy(self):
        f = open('0.orig/k', 'w+')

        f.write('/*--------------------------------*- C++ -*----------------------------------*\\   \n')
        f.write('| =========                 |                                                 |    \n')
        f.write('| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |  \n')
        f.write('|  \\\\    /   O peration     | Version:  v2112                                 |  \n')
        f.write('|   \\\\  /    A nd           | Website:  www.openfoam.com                      |  \n')
        f.write('|    \\\\/     M anipulation  |                                                 |  \n')
        f.write('\\*---------------------------------------------------------------------------*/   \n')
        f.write('FoamFile                                                                           \n')
        f.write('{                                                                                  \n')
        f.write('    version     2.0;                                                               \n')
        f.write('    format      ascii;                                                             \n')
        f.write('    class       volScalarField;                                                    \n')
        f.write('    object      k;                                                                 \n')
        f.write('}                                                                                  \n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //    \n')
        f.write('                                                                                   \n')
        f.write('kInlet         {}; \n'.format(self.k))
        f.write('                                                                                   \n')
        f.write('dimensions     [0 2 -2 0 0 0 0];                                                   \n')
        f.write('                                                                                   \n')
        f.write('internalField  uniform $kInlet;                                                    \n')
        f.write('                                                                                   \n')
        f.write('boundaryField                                                                      \n')
        f.write('{                                                                                  \n')
        f.write('   freestream                                                                      \n')
        f.write('   {                                                                               \n')
        f.write('       type            inletOutlet;                                                \n')
        f.write('       inletValue      uniform $kInlet;                                            \n')
        f.write('       value           uniform $kInlet;                                            \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   wall                                                                            \n')
        f.write('   {                                                                               \n')
        f.write('       type            kqRWallFunction;                                            \n')
        f.write('       value           uniform $kInlet;                                            \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   #includeEtc "caseDicts/setConstraintTypes"                                      \n')
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('// ************************************************************************* //    \n')
        f.close()

    def writeToFile_TurbulentThermalDiffusivity(self):
        f = open('0.orig/alphat', 'w+')

        f.write('/*--------------------------------*- C++ -*----------------------------------*\\   \n')
        f.write('| =========                 |                                                 |    \n')
        f.write('| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |  \n')
        f.write('|  \\\\    /   O peration     | Version:  v2112                                 |  \n')
        f.write('|   \\\\  /    A nd           | Website:  www.openfoam.com                      |  \n')
        f.write('|    \\\\/     M anipulation  |                                                 |  \n')
        f.write('\\*---------------------------------------------------------------------------*/   \n')
        f.write('FoamFile                                                                           \n')
        f.write('{                                                                                  \n')
        f.write('    version     2.0;                                                               \n')
        f.write('    format      ascii;                                                             \n')
        f.write('    class       volScalarField;                                                    \n')
        f.write('    object      alphat;                                                            \n')
        f.write('}                                                                                  \n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //    \n')
        f.write('                                                                                   \n')
        f.write('dimensions     [1 -1 -1 0 0 0 0];                                                  \n')
        f.write('                                                                                   \n')
        f.write('internalField  uniform 0;                                                          \n')
        f.write('                                                                                   \n')
        f.write('boundaryField                                                                      \n')
        f.write('{                                                                                  \n')
        f.write('   freestream                                                                      \n')
        f.write('   {                                                                               \n')
        f.write('       type            calculated;                                                 \n')
        f.write('       value           uniform 0;                                                  \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   wall                                                                            \n')
        f.write('   {                                                                               \n')
        f.write('       type            compressible::alphatWallFunction;                           \n')
        f.write('       value           uniform 0;                                                  \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   #includeEtc "caseDicts/setConstraintTypes"                                      \n')
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('// ************************************************************************* //    \n')
        f.close()

    def writeToFile_Velocity(self):
        f = open('0.orig/U', 'w+')

        f.write('/*--------------------------------*- C++ -*----------------------------------*\\   \n')
        f.write('| =========                 |                                                 |    \n')
        f.write('| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |  \n')
        f.write('|  \\\\    /   O peration     | Version:  v2112                                 |  \n')
        f.write('|   \\\\  /    A nd           | Website:  www.openfoam.com                      |  \n')
        f.write('|    \\\\/     M anipulation  |                                                 |  \n')
        f.write('\\*---------------------------------------------------------------------------*/   \n')
        f.write('FoamFile                                                                           \n')
        f.write('{                                                                                  \n')
        f.write('    version     2.0;                                                               \n')
        f.write('    format      ascii;                                                             \n')
        f.write('    class       volVectorField;                                                    \n')
        f.write('    object      U;                                                                 \n')
        f.write('}                                                                                  \n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //    \n')
        f.write('                                                                                   \n')
        f.write('Uinlet         ({} 0 0); \n'.format(self.u))
        f.write('                                                                                   \n')
        f.write('dimensions     [0 1 -1 0 0 0 0];                                                   \n')
        f.write('                                                                                   \n')
        f.write('internalField  uniform $Uinlet;                                                    \n')
        f.write('                                                                                   \n')
        f.write('boundaryField                                                                      \n')
        f.write('{                                                                                  \n')
        f.write('   freestream                                                                      \n')
        f.write('   {                                                                               \n')
        f.write('       type            freestreamVelocity;                                         \n')
        f.write('       freestreamValue uniform $Uinlet;                                            \n')
        f.write('       value           uniform $Uinlet;                                            \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   wall                                                                            \n')
        f.write('   {                                                                               \n')
        f.write('       type            noSlip;                                                     \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   #includeEtc "caseDicts/setConstraintTypes"                                      \n')
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('// ************************************************************************* //    \n')
        f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate Initial conditions files')
    parser.add_argument('mach', type=float, help='Freestream Mach number [-]')
    args = parser.parse_args()

    generateInitialConditions(args.mach)