#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------
# Created By  : Sebastian Widmann
# Institution : TU Munich, Department of Aerospace and Geodesy
# Created Date: June 16, 2022
# version ='1.0'
# ---------------------------------------------------------------------------
"""
Class implementation to generate 2D shells of NACA airfoils for blockMesh
process. Airfoil is specified by either 4- or 5-digit code together with the
given angle of attack. blockMeshDict written into the "/system" directory
"""
# ---------------------------------------------------------------------------

import argparse
import numpy as np
from math import ceil
import matplotlib.pyplot as plt

from flowProperties import *


class generateBlockMeshDict(object):
    def __init__(self, airfoil, angle, mach):
        self.foil = airfoil  # string used for creating Wavefront object files
        self.airfoil = [int(i) for i in airfoil]
        self.aoa = angle
        self.alpha = np.deg2rad(angle)
        self.mach = mach

        # Airfoil parameter
        self.nPoints = 400
        self.points = np.zeros((2 * self.nPoints, 2))
        self.chord = 1
        self.m = self.p = self.q = self.t = None
        self.xLead = self.xTrail = self.xUpper = self.xLower = None
        self.xLead = self.zTrail = self.zUpper = self.zLower = None
        self.spline47 = self.spline75 = self.spline48 = self.spline85 = None

        # Mesh parameter
        self.xMin = self.xMax = None
        self.zMin = self.zMax = None
        self.yMin = self.yMax = None
        self.xCellsB0 = self.xCellsB1 = self.xCellsB2 = None
        self.xCellsB3 = self.xCellsB4 = self.xCellsB5 = None
        self.xGradingB0 = self.xGradingB1 = self.xGradingB2 = None
        self.xGradingB3 = self.xGradingB4 = self.xGradingB5 = None
        self.zCells = self.zGrading = None

    def generateAirfoil(self):
        if len(self.airfoil) == 4:
            self.m = 0.01 * self.airfoil[0]
            self.p = 0.1 * self.airfoil[1]
            self.t = 0.01 * (10 * self.airfoil[2] + self.airfoil[3])
        elif len(self.airfoil) == 5:
            self.m = 0.15 * self.airfoil[0]
            self.p = 0.05 * self.airfoil[1]
            self.q = self.airfoil[2]
            self.t = 0.01 * (10 * self.airfoil[3] + self.airfoil[4])

        beta = np.linspace(0, np.pi, self.nPoints)
        x = 0.5 * (1 - np.cos(beta))

        # Coefficients to calculate half thickness z_t
        a0 = 0.2969
        a1 = -0.126
        a2 = -0.3516
        a3 = 0.2843
        a4 = -0.1036  # definition for closed trailing edge

        z_t = None
        z_c = None
        dz_c = None

        if len(self.airfoil) == 4:
            # Calculate half thickness (centerline to surface)
            z_t = 5 * self.t * (
                    a0 * np.sqrt(x) + a1 * x + a2 * np.power(x, 2) + a3 * np.power(x, 3) + a4 * np.power(x, 4))

            # Calculate mean camber line and gradient of mean camber line
            # if-condition for symmetrical, else-condition for cambered airfoils
            if self.p == 0:
                z_c = np.zeros(x.shape)
                dz_c = np.zeros(x.shape)
            else:
                z_c = np.where(x < self.p, self.m / np.power(self.p, 2) * (2 * self.p * x - np.power(x, 2)),
                               self.m / np.power(1 - self.p, 2) * (1 - 2 * self.p + 2 * self.p * x - np.power(x, 2)))
                dz_c = np.where(x < self.p, 2 * self.m / np.power(self.p, 2) * (self.p - x),
                                2 * self.m / np.power(1 - self.p, 2) * (self.p - x))

        elif len(self.airfoil) == 5:
            # Calculate half thickness (centerline to surface)
            z_t = 5 * self.t * (
                    a0 * np.sqrt(x) + a1 * x + a2 * np.power(x, 2) + a3 * np.power(x, 3) + a4 * np.power(x, 4))

            # Values for constants r, k1, k2/k1
            digits = np.array([10, 20, 30, 40, 50, 21, 31, 41, 51])
            r_values = np.array([0.0580, 0.1260, 0.2025, 0.2900, 0.3910, 0.1300, 0.2170, 0.3180, 0.4410])
            k1_values = np.array([361.400, 51.640, 15.957, 6.643, 3.230, 51.990, 15.793, 6.520, 3.191])
            k2k1_values = np.array([0, 0, 0, 0, 0, 0.000764, 0.00677, 0.0303, 0.1355])

            # Find array index for corresponding PQ-value
            index = np.where(digits == (10 * self.airfoil[1] + self.airfoil[2]))
            r = r_values[index]
            k1 = k1_values[index]
            k2k1 = k2k1_values[index]

            # Calculate mean camber line
            if self.q == 0:
                z_c = np.where(x < r, k1 / 6 * (np.power(x, 3) - 3 * r * np.power(x, 2) + np.power(r, 2) * x * (3 - r)),
                               k1 * np.power(r, 3) / 6 * (1 - x))
                dz_c = np.where(x < r, k1 / 6 * (3 * np.power(x, 2) - 6 * r * x + (3 - r) * np.power(r, 2)),
                                -k1 * np.power(r, 3) / 6)
            else:
                z_c = np.where(x < r, k1 / 6 * (
                        np.power(x - r, 3) - k2k1 * x * np.power(1 - r, 3) - x * np.power(r, 3) + np.power(r, 3)),
                               k1 / 6 * (k2k1 * np.power(x - r, 3) - k2k1 * x * np.power(1 - r, 3) - x * np.power(r,
                                                                                                                  3) + np.power(
                                   r, 3)))
                dz_c = np.where(x < r, k1 / 6 * (3 * np.power(x - r, 2) - k2k1 * np.power(1 - r, 3) - np.power(r, 3)),
                                k1 / 6 * (3 * k2k1 * np.power(x - r, 2) - k2k1 * np.power(1 - r, 3) - np.power(r, 3)))

        theta = np.arctan(dz_c)

        x_upper = x - z_t * np.sin(theta)
        x_lower = x + z_t * np.sin(theta)

        z_upper = z_c + z_t * np.cos(theta)
        z_lower = z_c - z_t * np.cos(theta)

        vertexX = np.concatenate((x_upper, x_lower), axis=0)
        vertexZ = np.concatenate((z_upper, z_lower), axis=0)

        self.points[:, 0], self.points[:, 1] = vertexX.transpose(), vertexZ.transpose()

        # # Translate aerodynamic centre (0.25c) to coordinate origin for rotation
        self.points[:, 0] -= 0.25 * self.chord
        rotationMatrix = np.array([[np.cos(self.alpha), np.sin(self.alpha)], [-np.sin(self.alpha), np.cos(self.alpha)]])

        self.points = np.einsum('ij, kj -> ki', rotationMatrix, self.points)

        # Access values of leading edge, trailing edge,
        # upper and lower points of max. thickness
        idx_min = np.argmin(self.points[:, 0])
        self.xLead = self.points[idx_min][0]
        self.zLead = self.points[idx_min][1]

        idx_max = np.argmax(self.points[:, 0])
        self.xTrail = self.points[idx_max][0]
        self.zTrail = self.points[idx_max][1]

        idz_max = np.argmax(z_t)
        self.xUpper = self.points[idz_max][0]
        self.zUpper = self.points[idz_max][1]
        self.xLower = self.points[idz_max + self.nPoints][0]
        self.zLower = self.points[idz_max + self.nPoints][1]

        # Create point subsets for splines within the different blocks
        self.spline47 = self.points[self.nPoints:idz_max + self.nPoints, ]
        self.spline75 = self.points[idz_max + self.nPoints:, ]
        self.spline48 = self.points[0:idz_max, ]
        self.spline85 = self.points[idz_max:self.nPoints, ]

    def setDomainSize(self, xMin=-20, xMax=40, zMin=-20, zMax=20):
        # Domain size parameter
        self.xMin = xMin
        self.xMax = xMax
        self.zMin = zMin
        self.zMax = zMax
        self.yMin = -0.1
        self.yMax = 0.1

    def calculateMinimumWallHeight(self):
        '''
            Calculate minimum cell height above airfoil to achieve wall y+ < 1
        '''
        yplus = 1

        p = calculateStaticPressure(self.mach)
        T = calculateStaticTemperature(self.mach)
        rho = calculateStaticDensity(p, T)
        u = self.mach * calculateSpeedofSound(T)
        Re = rho * u * self.chord / mu  # [-] Freestream Reynolds number
        cf = np.power(2 * np.log10(Re) - 0.65, -2.3)  # [-] Skin friction coefficient based on Schlichting
        Tau_w = 0.5 * cf * rho * u ** 2  # [Pa] Wall shear stress
        u_star = np.sqrt(Tau_w / rho)  # [m*s^-1] Friction velocity
        return yplus * mu / (rho * u_star)  # First layer thickness in z-direction

    def calculateNumberOfCells(self, cRatio, cHeight, Height):
        if cRatio == 1:
            return ceil(Height / cHeight)
        else:
            return ceil((np.log(1 - Height / cHeight * (1 - cRatio))) / np.log(cRatio))

    def calculateExpansionRatio(self, cRatio, nCells):
        return np.power(cRatio, nCells - 1)

    def calculateMaximumCellHeight(self, cHeight, ratio):
        return ratio * cHeight

    def setMeshSize(self):
        cRatio = 1.05  # [-] Cell-to-cell expansion ratio

        # Mesh parameter from airfoil walls to inlet in radial direction
        self.zCells = self.calculateNumberOfCells(cRatio, self.calculateMinimumWallHeight(), self.zMax)
        self.zGrading = self.calculateExpansionRatio(cRatio, self.zCells)

        # Mesh parameter from airfoil walls to inlet in azimuth direction for blocks 0, 3
        xACells = 60  # Total cells along inlet arc
        self.xCellsB3 = round(xACells * ((0.5 * np.pi - self.alpha) / np.pi))  # upstream cells in block 3
        self.xCellsB0 = xACells - self.xCellsB3  # upstream cells in block 1

        self.xGradingB0 = 5
        self.xGradingB3 = 5

        # Mesh parameter from airfoil walls to far field in x-direction for blocks 1, 4
        self.xCellsB1 = 50
        self.xCellsB4 = 30

        self.xGradingB1 = 1
        self.xGradingB4 = 1

        # Mesh parameter from trailing edge to outlet in x-direction for blocks 2, 6
        self.xCellsB2 = 100
        self.xCellsB5 = self.xCellsB2

        self.xGradingB2 = 10
        self.xGradingB5 = self.xGradingB2

        self.zCells = 80  # aerofoil to far field
        self.xUCells = 30  # upstream
        self.xMCells = 30  # middle
        self.xDCells = 40  # downstream

        self.zGrading = 40  # aerofoil to far field
        self.xUGrading = 5  # towards centre upstream
        self.leadGrading = 0.2  # towards leading edge
        self.xDGrading = 10  # downstream
        #
        # self.xUCells = self.calculateNumberOfCells(cRatio, 0.01,
        #                                            np.abs(self.xMin))  # Cells between inlet and leading edge
        # self.xUGrading = 5  # towards centre upstream
        #
        # self.xMCells = self.calculateNumberOfCells(1, 0.01, np.abs(
        #     self.xTrail - self.xUpper))  # Cells between max camber and trailing edge
        #
        # self.xDCells = self.calculateNumberOfCells(cRatio, 0.01, self.xMax)  # Cells from trailing edge to outlet
        # self.xDGrading = self.calculateExpansionRatio(cRatio, self.xDCells)  # downstream
        #
        # self.leadGrading = 5  # towards leading edge

    def pltShow(self):
        plt.scatter(self.points[:, 0], self.points[:, 1], s=1, color='black')
        plt.scatter(self.xLead, self.zLead, s=0.5, color='red')
        plt.scatter(self.xUpper, self.zUpper, s=0.5, color='red')
        plt.scatter(self.xTrail, self.zTrail, s=0.5, color='red')
        plt.scatter(self.xLower, self.zLower, s=0.5, color='red')
        plt.ylim(-0.5, 0.5)
        plt.show()

    def writeToFile(self):
        f = open('system/blockMeshDict', 'w+')

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
        f.write('    object      blockMeshDict;                                                     \n')
        f.write('}                                                                                  \n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //    \n')
        f.write('                                                                                   \n')
        f.write('domain                                                                             \n')
        f.write('{                                                                                  \n')
        f.write('   xMin        {}; \n'.format(self.xMin))
        f.write('   xMax        {}; \n'.format(self.xMax))
        f.write('   zMin        {}; \n'.format(self.zMin))
        f.write('   zMax        {}; \n'.format(self.zMax))
        f.write('                                                                                   \n')
        f.write('   // Number of cells                                                              \n')
        f.write('   xCellsB0    {}; // leading edge to xLower \n'.format(self.xCellsB0))
        f.write('   xCellsB1    {}; // xLower to trailing edge \n'.format(self.xCellsB1))
        f.write('   xCellsB2    {}; // trailing edge to outlet \n'.format(self.xCellsB2))
        f.write('   xCellsB3    {}; // leading edge to xUpper \n'.format(self.xCellsB3))
        f.write('   xCellsB4    {}; // xUpper to trailing edge \n'.format(self.xCellsB4))
        f.write('   xCellsB5    {}; // trailing edge to outlet \n'.format(self.xCellsB5))
        f.write('   zCells      {}; // aerofoil to far field \n'.format(self.zCells))
        f.write('                                                                                   \n')
        f.write('   // Mesh grading                                                                 \n')
        f.write('   xGradingB0  {}; // from 0 to 3 \n'.format(self.xGradingB0))
        f.write('   xGradingB1  {}; // from 0 to 1 \n'.format(self.xGradingB1))
        f.write('   xGradingB2  {}; // from 1 to 2 \n'.format(self.xGradingB2))
        f.write('   xGradingB3  {}; // from 9 to 3 \n'.format(self.xGradingB3))
        f.write('   xGradingB4  {}; // from 9 to 10 \n'.format(self.xGradingB4))
        f.write('   xGradingB5  {}; // from 10 to 11 \n'.format(self.xGradingB5))
        f.write('   zGrading    {}; // from airfoil to far field \n'.format(self.zGrading))
        f.write('   leadGrading {}; // from xUpper/xLower to LE \n'.format(self.leadGrading))
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('aerofoil                                                                           \n')
        f.write('{                                                                                  \n')
        f.write('   xLead       {}; \n'.format(self.xLead))
        f.write('   zLead       {}; \n'.format(self.zLead))
        f.write('   xTrail      {}; \n'.format(self.xTrail))
        f.write('   zTrail      {}; \n'.format(self.zTrail))
        f.write('   xUpper      {}; \n'.format(self.xUpper))
        f.write('   zUpper      {}; \n'.format(self.zUpper))
        f.write('   xLower      {}; \n'.format(self.xLower))
        f.write('   zLower      {}; \n'.format(self.zLower))
        f.write('   zInlet      {}; \n'.format(np.abs(self.xMin) * np.tan(self.alpha)))
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('geometry                                                                           \n')
        f.write('{                                                                                  \n')
        f.write('   // aerofoil = NACA{} \n'.format(self.foil))
        f.write('   // AoA      = {} deg \n'.format(self.aoa))
        f.write('                                                                                   \n')
        f.write('   cylinder                                                                        \n')
        f.write('   {                                                                               \n')
        f.write('       type    cylinder;                                                           \n')
        f.write('       point1  ($:aerofoil.xUpper -1e3 0);                                         \n')
        f.write('       point2  ($:aerofoil.xUpper  1e3 0);                                         \n')
        f.write('       radius  $:domain.zMax;                                                      \n')
        f.write('   }                                                                               \n')
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('vertices                                                                           \n')
        f.write('(                                                                                  \n')
        f.write('   project ($aerofoil.xLower   -0.1 $domain.zMin)      (cylinder)                  \n')
        f.write('           ($aerofoil.xTrail   -0.1 $domain.zMin)                                  \n')
        f.write('           ($domain.xMax       -0.1 $domain.zMin)                                  \n')
        f.write('                                                                                   \n')
        f.write('   project ($domain.xMin       -0.1 $aerofoil.zInlet)   (cylinder)                 \n')
        f.write('           ($aerofoil.xLead    -0.1 $aerofoil.zLead)                               \n')
        f.write('           ($aerofoil.xTrail   -0.1 $aerofoil.zTrail)                              \n')
        f.write('           ($domain.xMax       -0.1 $aerofoil.zTrail)                              \n')
        f.write('                                                                                   \n')
        f.write('           ($aerofoil.xLower   -0.1 $aerofoil.zLower)                              \n')
        f.write('           ($aerofoil.xUpper   -0.1 $aerofoil.zUpper)                              \n')
        f.write('                                                                                   \n')
        f.write('           ($aerofoil.xUpper   -0.1 $domain.zMax)                                  \n')
        f.write('           ($aerofoil.xTrail   -0.1 $domain.zMax)                                  \n')
        f.write('           ($domain.xMax       -0.1 $domain.zMax)                                  \n')
        f.write('                                                                                   \n')
        f.write('   project ($aerofoil.xLower   0.1 $domain.zMin)       (cylinder)                  \n')
        f.write('           ($aerofoil.xTrail   0.1 $domain.zMin)                                   \n')
        f.write('           ($domain.xMax       0.1 $domain.zMin)                                   \n')
        f.write('                                                                                   \n')
        f.write('   project ($domain.xMin       0.1 $aerofoil.zInlet)    (cylinder)                 \n')
        f.write('           ($aerofoil.xLead    0.1 $aerofoil.zLead)                                \n')
        f.write('           ($aerofoil.xTrail   0.1 $aerofoil.zTrail)                               \n')
        f.write('           ($domain.xMax       0.1 $aerofoil.zTrail)                               \n')
        f.write('                                                                                   \n')
        f.write('           ($aerofoil.xLower   0.1 $aerofoil.zLower)                               \n')
        f.write('           ($aerofoil.xUpper   0.1 $aerofoil.zUpper)                               \n')
        f.write('                                                                                   \n')
        f.write('           ($aerofoil.xUpper   0.1 $domain.zMax)                                   \n')
        f.write('           ($aerofoil.xTrail   0.1 $domain.zMax)                                   \n')
        f.write('           ($domain.xMax       0.1 $domain.zMax)                                   \n')
        f.write(');                                                                                 \n')
        f.write('                                                                                   \n')
        f.write('blocks                                                                             \n')
        f.write('(                                                                                  \n')
        f.write('   hex (7 4 16 19 0 3 15 12)                                                       \n')
        f.write('   ($:domain.xCellsB0 1 $:domain.zCells)                                           \n')
        f.write('   edgeGrading                                                                     \n')
        f.write('   (                                                                               \n')
        f.write('       $:domain.leadGrading $:domain.leadGrading $:domain.xGradingB0 $:domain.xGradingB0 \n')
        f.write('       1 1 1 1                                                                     \n')
        f.write('       $:domain.zGrading $:domain.zGrading $:domain.zGrading $:domain.zGrading     \n')
        f.write('   )                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   hex (5 7 19 17 1 0 12 13)                                                       \n')
        f.write('   ($:domain.xCellsB1 1 $:domain.zCells)                                           \n')
        f.write('   simpleGrading ($:domain.xGradingB1 1 $:domain.zGrading)                         \n')
        f.write('                                                                                   \n')
        f.write('   hex (17 18 6 5 13 14 2 1)                                                       \n')
        f.write('   ($:domain.xCellsB2 1 $:domain.zCells)                                           \n')
        f.write('   simpleGrading ($:domain.xGradingB2 1 $:domain.zGrading)                         \n')
        f.write('                                                                                   \n')
        f.write('   hex (20 16 4 8 21 15 3 9)                                                       \n')
        f.write('   ($:domain.xCellsB3 1 $:domain.zCells)                                           \n')
        f.write('   edgeGrading                                                                     \n')
        f.write('   (                                                                               \n')
        f.write('       $:domain.leadGrading $:domain.leadGrading $:domain.xGradingB3 $:domain.xGradingB3 \n')
        f.write('       1 1 1 1                                                                     \n')
        f.write('       $:domain.zGrading $:domain.zGrading $:domain.zGrading $:domain.zGrading     \n')
        f.write('   )                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   hex (17 20 8 5 22 21 9 10)                                                      \n')
        f.write('   ($:domain.xCellsB4 1 $:domain.zCells)                                           \n')
        f.write('   simpleGrading ($:domain.xGradingB4 1 $:domain.zGrading)                         \n')
        f.write('                                                                                   \n')
        f.write('   hex (5 6 18 17 10 11 23 22)                                                     \n')
        f.write('   ($:domain.xCellsB5 1 $:domain.zCells)                                           \n')
        f.write('   simpleGrading ($:domain.xGradingB5 1 $:domain.zGrading)                         \n')
        f.write(');                                                                                 \n')
        f.write('                                                                                   \n')
        f.write('edges                                                                              \n')
        f.write('(                                                                                  \n')
        f.write('   spline 4 7                                                                      \n')
        f.write('   (                                                                               \n')
        for row in self.spline47:
            f.write('        ({} {} {})\n'.format(row[0], self.yMin, row[1]))
        f.write('   )                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   spline 7 5                                                                      \n')
        f.write('   (                                                                               \n')
        for row in self.spline75:
            f.write('        ({} {} {})\n'.format(row[0], self.yMin, row[1]))
        f.write('   )                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   spline 4 8                                                                     \n')
        f.write('   (                                                                               \n')
        for row in self.spline48:
            f.write('        ({} {} {})\n'.format(row[0], self.yMin, row[1]))
        f.write('   )                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   spline 8 5                                                                      \n')
        f.write('   (                                                                               \n')
        for row in self.spline85:
            f.write('        ({} {} {})\n'.format(row[0], self.yMin, row[1]))
        f.write('   )                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   spline 16 19                                                                    \n')
        f.write('   (                                                                               \n')
        for row in self.spline47:
            f.write('        ({} {} {})\n'.format(row[0], self.yMax, row[1]))
        f.write('   )                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   spline 19 17                                                                    \n')
        f.write('   (                                                                               \n')
        for row in self.spline75:
            f.write('        ({} {} {})\n'.format(row[0], self.yMax, row[1]))
        f.write('   )                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   spline 16 20                                                                    \n')
        f.write('   (                                                                               \n')
        for row in self.spline48:
            f.write('        ({} {} {})\n'.format(row[0], self.yMax, row[1]))
        f.write('   )                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   spline 20 17                                                                    \n')
        f.write('   (                                                                               \n')
        for row in self.spline85:
            f.write('        ({} {} {})\n'.format(row[0], self.yMax, row[1]))
        f.write('   )                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   project 3 0 (cylinder)                                                          \n')
        f.write('   project 3 9 (cylinder)                                                          \n')
        f.write('   project 15 12 (cylinder)                                                        \n')
        f.write('   project 15 21 (cylinder)                                                        \n')
        f.write(');                                                                                 \n')
        f.write('                                                                                   \n')
        f.write('boundary                                                                           \n')
        f.write('(                                                                                  \n')
        f.write('   aerofoil                                                                        \n')
        f.write('   {                                                                               \n')
        f.write('       type wall;                                                                  \n')
        f.write('       faces                                                                       \n')
        f.write('       (                                                                           \n')
        f.write('           (4 7 19 16)                                                             \n')
        f.write('           (7 5 17 19)                                                             \n')
        f.write('           (5 8 20 17)                                                             \n')
        f.write('           (8 4 16 20)                                                             \n')
        f.write('       );                                                                          \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   inlet                                                                           \n')
        f.write('   {                                                                               \n')
        f.write('       type patch;                                                                 \n')
        f.write('       inGroups (freestream);                                                      \n')
        f.write('       faces                                                                       \n')
        f.write('       (                                                                           \n')
        f.write('           (3 0 12 15)                                                             \n')
        f.write('           (0 1 13 12)                                                             \n')
        f.write('           (1 2 14 13)                                                             \n')
        f.write('           (11 10 22 23)                                                           \n')
        f.write('           (10 9 21 22)                                                            \n')
        f.write('           (9 3 15 21)                                                             \n')
        f.write('       );                                                                          \n')
        f.write('                                                                                   \n')
        f.write('   }                                                                               \n')
        f.write('   outlet                                                                          \n')
        f.write('   {                                                                               \n')
        f.write('       type patch;                                                                 \n')
        f.write('       inGroups (freestream);                                                      \n')
        f.write('       faces                                                                       \n')
        f.write('       (                                                                           \n')
        f.write('           (2 6 18 14)                                                             \n')
        f.write('           (6 11 23 18)                                                            \n')
        f.write('       );                                                                          \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   back                                                                            \n')
        f.write('   {                                                                               \n')
        f.write('       type empty;                                                                 \n')
        f.write('       faces                                                                       \n')
        f.write('       (                                                                           \n')
        f.write('           (3 4 7 0)                                                               \n')
        f.write('           (7 5 1 0)                                                               \n')
        f.write('           (5 6 2 1)                                                               \n')
        f.write('           (3 9 8 4)                                                               \n')
        f.write('           (9 10 5 8)                                                              \n')
        f.write('           (10 11 6 5)                                                             \n')
        f.write('       );                                                                          \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   front                                                                           \n')
        f.write('   {                                                                               \n')
        f.write('       type empty;                                                                 \n')
        f.write('       faces                                                                       \n')
        f.write('       (                                                                           \n')
        f.write('           (15 16 19 12)                                                           \n')
        f.write('           (19 17 13 12)                                                           \n')
        f.write('           (17 18 14 13)                                                           \n')
        f.write('           (15 16 20 21)                                                           \n')
        f.write('           (20 17 22 21)                                                           \n')
        f.write('           (17 18 23 22)                                                           \n')
        f.write('       );                                                                          \n')
        f.write('   }                                                                               \n')
        f.write(');                                                                                 \n')
        f.write('                                                                                   \n')
        f.write('// ************************************************************************* //    \n')
        f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate blockMeshDict File')
    parser.add_argument('airfoil', type=str, help='NACA airfoil digits')
    parser.add_argument('aoa', type=float, help='Angle of attack [deg]')
    parser.add_argument('mach', type=float, help='Mach number [-]')
    args = parser.parse_args()

    # error handling for incorrect parsing of NACA airfoil or AoA
    if len(args.airfoil) < 4 or len(args.airfoil) > 5:
        parser.error('Please enter a 4- or 5-digit NACA airfoil.')
    if np.deg2rad(args.aoa) <= -np.pi / 2 or np.deg2rad(args.aoa) >= np.pi / 2:
        parser.error('Please enter an angle of attack between ({}, {}) degree.'.format(-90, 90))

    blockMesh = generateBlockMeshDict(args.airfoil, args.aoa, args.mach)

    blockMesh.generateAirfoil()
    blockMesh.setDomainSize()
    blockMesh.setMeshSize()

    blockMesh.writeToFile()
