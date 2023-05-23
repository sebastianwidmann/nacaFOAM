# ----------------------------------------------------------------------------
# Created By  : Sebastian Widmann
# Institution : TU Munich, Department of Aerospace and Geodesy
# Created Date: August 31, 2022
# version ='1.0'
# ---------------------------------------------------------------------------
"""
Class implementation to generate explicitly use surfaceExtractFeature by
creating the ".eMesh" file located in "constant/triSurface" directory
for improved surface snapping of snappyHexMesh.
"""
# ---------------------------------------------------------------------------

import os

class generateSurfaceFeatureExtractDict(object):
    def __init__(self, caseDir, airfoil, angle):
        self.caseDir = caseDir
        self.airfoil = airfoil
        self.angle = angle

        self.writeToFile()

    def writeToFile(self):

        saveDir = os.path.join(self.caseDir, 'system/surfaceFeatureExtractDict')
        f = open(saveDir, 'w+')

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
        f.write('    object      surfaceFeatureExtractDict;                                         \n')
        f.write('}                                                                                  \n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //    \n')
        f.write('                                                                                   \n')
        f.write('naca{}_{}.stl \n'.format(self.airfoil, self.angle))
        f.write('{                                                                                  \n')
        f.write('   extractionMethod        extractFromSurface;                                     \n')
        f.write('                                                                                   \n')
        f.write('   extractFromSurfaceCoeffs                                                        \n')
        f.write('   {                                                                               \n')
        f.write('       includedAngle       150;                                                    \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   subsetFeatures                                                                  \n')
        f.write('   {                                                                               \n')
        f.write('       nonManifoldEdges    yes;                                                    \n')
        f.write('       openEdges           yes;                                                    \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   writeObj                no;                                                     \n')
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('// ************************************************************************* //    \n')
        f.close()
