# ----------------------------------------------------------------------------
# Created By  : Sebastian Widmann
# Institution : TU Munich, Department of Aerospace and Geodesy
# Created Date: September 12, 2022
# version ='1.0'
# ---------------------------------------------------------------------------
"""
Class implementation to generate fvOptions. During meshing process fvOptions should not exist
as displacement motion solver will try to enforce temperature limit and hence crash. fvOptions
will be written into the "/system" directory.
"""
# ---------------------------------------------------------------------------

import os

class generateFvOptionsDict(object):
    def __init__(self, caseDir):
        self.caseDir = caseDir
        self.writeToFile()

    def writeToFile(self):
        saveDir = os.path.join(self.caseDir, 'system/fvOptions')
        f = open(saveDir, 'w+')

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
        f.write('    class       dictionary;                                                        \n')
        f.write('    object      fvOptions;                                                         \n')
        f.write('}                                                                                  \n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //    \n')
        f.write('                                                                                   \n')
        f.write('limitT                                                                             \n')
        f.write('{                                                                                  \n')
        f.write('   type            limitTemperature;                                               \n')
        f.write('   min             273.15;                                                         \n')
        f.write('   max             500;                                                            \n')
        f.write('   selectionMode   all;                                                            \n')
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('// ************************************************************************* //    \n')
        f.close()
