#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------
# Created By  : Sebastian Widmann
# Institution : TU Munich, Department of Aerospace and Geodesy
# Created Date: September 21, 2022
# version ='1.0'
# ---------------------------------------------------------------------------
"""
Class implementation to generate controlDict. controlDict will be written
into the "/system" directory.
"""
# ---------------------------------------------------------------------------
import argparse
from flowProperties import *

class generateForceCoefficients(object):
    def __init__(self):
        self.p = calculateStaticPressure(args.mach)
        self.T = calculateStaticTemperature(args.mach)
        self.rho = calculateStaticDensity(self.p, self.T)
        self.u = args.mach * calculateSpeedofSound(self.T)

        self.writeToFile()

    def writeToFile(self):
        f = open('system/foForceCoeffs', 'w+')

        f.write('/*--------------------------------*- C++ -*----------------------------------*\\   \n')
        f.write('| =========                 |                                                 |    \n')
        f.write('| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |  \n')
        f.write('|  \\\\    /   O peration     | Version:  v2206                                 |  \n')
        f.write('|   \\\\  /    A nd           | Website:  www.openfoam.com                      |  \n')
        f.write('|    \\\\/     M anipulation  |                                                 |  \n')
        f.write('\\*---------------------------------------------------------------------------*/   \n')
        f.write('FoamFile                                                                           \n')
        f.write('{                                                                                  \n')
        f.write('    version     2.0;                                                               \n')
        f.write('    format      ascii;                                                             \n')
        f.write('    class       dictionary;                                                        \n')
        f.write('    object      forceCoeffs;                                                       \n')
        f.write('}                                                                                  \n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //    \n')
        f.write('                                                                                   \n')
        f.write('forceCoeffs                                                                        \n')
        f.write('{                                                                                  \n')
        f.write('   type            forceCoeffs;                                                    \n')
        f.write('   libs            ("libforces.so");                                               \n')
        f.write('   patches         (wall);                                                         \n')
        f.write('                                                                                   \n')
        f.write('   writeControl    timeStep;                                                       \n')
        f.write('   writeInterval   1;                                                              \n')
        f.write('   executeControl  timeStep;                                                       \n')
        f.write('   executeInterval 1;                                                              \n')
        f.write('   timeStart       0;                                                              \n')
        f.write('                                                                                   \n')
        f.write('   p               p;                                                              \n')
        f.write('   U               U;                                                              \n')
        f.write('   rho             rho;                                                            \n')
        f.write('   pRef            {}; \n'.format(self.p))
        f.write('   rhoInf          {}; \n'.format(self.rho))
        f.write('   writeFields     no;                                                             \n')
        f.write('   CofR            (0 0 0);                                                        \n')
        f.write('   liftDir         (0 1 0);                                                        \n')
        f.write('   dragDir         (1 0 0);                                                        \n')
        f.write('   pitchAxis       (0 0 1);                                                        \n')
        f.write('   magUInf         {}; \n'.format(self.u))
        f.write('   lRef            1;                                                              \n')
        f.write('   Aref            1;                                                              \n')
        f.write('                                                                                   \n')
        f.write('   binData                                                                         \n')
        f.write('   {                                                                               \n')
        f.write('       nBin        20;                                                             \n')
        f.write('       direction   (1 0 0);                                                        \n')
        f.write('       cumulative  yes;                                                            \n')
        f.write('   }                                                                               \n')
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('// ************************************************************************* //    \n')
        f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate functions sub-dictionary for controlDict. File is saved into "system/foForceCoeffs"')
    parser.add_argument('mach', type=float, help='Freestream Mach number [-]')
    args = parser.parse_args()

    generateForceCoefficients()
