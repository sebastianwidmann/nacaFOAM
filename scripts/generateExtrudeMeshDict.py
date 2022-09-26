#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------
# Created By  : Sebastian Widmann
# Institution : TU Munich, Department of Aerospace and Geodesy
# Created Date: September 16, 2022
# version ='1.0'
# ---------------------------------------------------------------------------
"""
Class implementation to generate extrudeMeshDict. Used to convert 3D mesh into
2D mesh. extrudeMeshDict will be written into the "/system" directory.
"""
# ---------------------------------------------------------------------------

import argparse


class generateExtrudeMeshDict(object):
    def __init__(self):
        self.writeToFile()

    def writeToFile(self):
        f = open('system/extrudeMeshDict', 'w+')

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
        f.write('    object      extrudeMeshDict;                                                   \n')
        f.write('}                                                                                  \n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //    \n')
        f.write('                                                                                   \n')
        f.write('constructFrom      patch;                                                          \n')
        f.write('sourceCase         ".";                                                            \n')
        f.write('                                                                                   \n')
        if args.meshing == 0:
            f.write('sourcePatches      (symFront);                                                 \n')
            f.write('exposedPatchName   symBack;                                                    \n')
        elif args.meshing == 1:
            f.write('sourcePatches      (symBack);                                                  \n')
            f.write('exposedPatchName   symFront;                                                   \n')
        f.write('                                                                                   \n')
        f.write('flipNormals        false;                                                          \n')
        f.write('                                                                                   \n')
        f.write('extrudeModel       linearDirection;                                                \n')
        f.write('                                                                                   \n')
        f.write('nLayers            1;                                                              \n')
        f.write('                                                                                   \n')
        f.write('expansionRatio     1.0;                                                            \n')
        f.write('                                                                                   \n')
        f.write('linearDirectionCoeffs                                                              \n')
        f.write('{                                                                                  \n')
        if args.meshing == 0:
            f.write('   direction       (0 0 1);                                                    \n')
            f.write('   thickness       {}; \n'.format(args.thickness))
        elif args.meshing == 1:
            f.write('   direction       (0 0 -1);                                                   \n')
            f.write('   thickness       0.5;                                                        \n')
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('mergeFaces         false;                                                          \n')
        f.write('                                                                                   \n')
        f.write('mergeTol           1e-9;                                                           \n')
        f.write('                                                                                   \n')
        f.write('// ************************************************************************* //    \n')
        f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate extrudeMeshDict file.')
    parser.add_argument('-thickness', type=float, help='extrusion thickness')
    parser.add_argument('-meshing', dest='meshing', type=int,
                        help='integer corresponding to meshing subpart \n 0 = generate refinementRegions w/o geometry '
                             '1 = generate castellatedMesh plus snapping around airfoil')
    args = parser.parse_args()

    generateExtrudeMeshDict()
