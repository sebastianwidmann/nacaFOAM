# ----------------------------------------------------------------------------
# Created By  : Sebastian Widmann
# Institution : TU Munich, Department of Aerospace and Geodesy
# Created Date: September 27, 2022
# version ='1.0'
# ---------------------------------------------------------------------------
"""
Class implementation to generate pressure function within controlDict.
pressureCoeffs will be written into the "/system" directory.
"""
# ---------------------------------------------------------------------------
import argparse
from flowProperties import *

class generatePressureCoefficients(object):
    def __init__(self):
        self.p = calculateStaticPressure(args.mach)
        self.T = calculateStaticTemperature(args.mach)
        self.rho = calculateStaticDensity(self.p, self.T)
        self.u = args.mach * calculateSpeedofSound(self.T)

        self.writeToFile()

    def writeToFile(self):
        f = open('system/foPressureCoeffs', 'w+')

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
        f.write('    object      pressure;                                                          \n')
        f.write('}                                                                                  \n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //    \n')
        f.write('                                                                                   \n')
        f.write('force                                                                              \n')
        f.write('{                                                                                  \n')
        f.write('   type            pressure;                                                       \n')
        f.write('   libs            (fieldFunctionObjects);                                         \n')
        f.write('   mode            staticCoeff;                                                    \n')
        f.write('                                                                                   \n')
        f.write('   writeFields     no;                                                             \n')
        f.write('                                                                                   \n')
        f.write('   p               p;                                                              \n')
        f.write('   U               U;                                                              \n')
        f.write('   rho             rho;                                                            \n')
        f.write('   pRef            {}; \n'.format(p0))
        f.write('   rhoInf          {}; \n'.format(self.rho))
        f.write('   hydroStaticMode none;                                                           \n')
        f.write('   g               (0 -9.81 0);                                                    \n')
        f.write('   hRef            0;                                                              \n')
        f.write('   pInf            {}; \n'.format(self.p))
        f.write('   UInf            ({} 0 0); \n'.format(self.u))
        f.write('   lRef            1;                                                              \n')
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('// ************************************************************************* //    \n')
        f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate functions sub-dictionary for controlDict. File is saved into "system/pressureCoeffs"')
    parser.add_argument('mach', type=float, help='Freestream Mach number [-]')
    args = parser.parse_args()

    generatePressureCoefficients()
