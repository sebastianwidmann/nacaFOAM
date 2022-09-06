#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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

import argparse
import numpy as np

class generateSurfaceFeatureExtract(object):
    def __init__(self, airfoil, angle):
        self.airfoil = airfoil
        self.alpha = angle

        self.writeToFile()

    def writeToFile(self):
        f = open('system/surfaceFeatureExtractDict', 'w+')

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
        f.write('    object      surfaceFeatureExtractDict;                                         \n')
        f.write('}                                                                                  \n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //    \n')
        f.write('                                                                                   \n')
        f.write('naca{}_{}deg.stl \n'.format(self.airfoil, self.alpha))
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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate snappyHexMeshDict file and save into "system/blockMeshFile"')
    parser.add_argument('airfoil', type=str, help='NACA airfoil digits')
    parser.add_argument('alpha', type=float, help='Angle of attack [deg]')
    args = parser.parse_args()

    # error handling for incorrect parsing of NACA airfoil or AoA
    if len(args.airfoil) < 4 or len(args.airfoil) > 5:
        parser.error('Please enter a 4- or 5-digit NACA airfoil.')
    if np.deg2rad(args.alpha) <= -np.pi / 2 or np.deg2rad(args.alpha) >= np.pi / 2:
        parser.error('Please enter an angle of attack between ({}, {}) degree.'.format(-90, 90))

    generateSurfaceFeatureExtract(args.airfoil, args.alpha)
