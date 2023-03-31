#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------
# Created By  : Sebastian Widmann
# Institution : TU Munich, Department of Aerospace and Geodesy
# Created Date: June 16, 2022
# version ='1.0'
# ---------------------------------------------------------------------------
"""
Class implementation to generate 3D body of NACA airfoil for snappyHexMesh
process. Airfoil is specified by either 4- or 5-digit code together with the
given angle of attack. File has type ".stl" and will be written into
the "/constant/triSurface" directory.
"""
# ---------------------------------------------------------------------------

import numpy as np
import os
from stl import mesh

class generateGeometryFile(object):
    def __init__(self, caseDir, airfoil, angle):
        self.caseDir = str(caseDir)
        self.foil = airfoil  # string used for creating Wavefront object files
        self.airfoil = [int(i) for i in airfoil]
        self.angle = angle
        self.alpha = np.deg2rad(angle)

        # Airfoil parameter
        self.nPoints = 500
        self.chord = 1
        self.m = self.p = self.q = self.t = None

        self.vertices = None
        self.faces = None

        self.generateVertices()
        self.generateFaces()
        self.generateStandardTriangleLanguageFile()

    def generateVertices(self):
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

        y_t = None
        y_c = None
        dy_c = None

        if len(self.airfoil) == 4:
            # Calculate half thickness (centerline to surface)
            y_t = 5 * self.t * (
                    a0 * np.sqrt(x) + a1 * x + a2 * np.power(x, 2) + a3 * np.power(x, 3) + a4 * np.power(x, 4))

            # Calculate mean camber line and gradient of mean camber line
            # if-condition for symmetrical, else-condition for cambered airfoils
            if self.p == 0:
                y_c = np.zeros(x.shape)
                dy_c = np.zeros(x.shape)
            else:
                y_c = np.where(x < self.p, self.m / np.power(self.p, 2) * (2 * self.p * x - np.power(x, 2)),
                               self.m / np.power(1 - self.p, 2) * (1 - 2 * self.p + 2 * self.p * x - np.power(x, 2)))
                dy_c = np.where(x < self.p, 2 * self.m / np.power(self.p, 2) * (self.p - x),
                                2 * self.m / np.power(1 - self.p, 2) * (self.p - x))

        elif len(self.airfoil) == 5:
            # Calculate half thickness (centerline to surface)
            y_t = 5 * self.t * (
                    a0 * np.sqrt(x) + a1 * x + a2 * np.power(x, 2) + a3 * np.power(x, 3) + a4 * np.power(x, 4))

            # Values for constants r, k1, k2/k1 at Cl = 0.3
            digits = np.array([10, 20, 30, 40, 50, 21, 31, 41, 51])
            r_values = np.array([0.0580, 0.1260, 0.2025, 0.2900, 0.3910, 0.1300, 0.2170, 0.3180, 0.4410])
            k1_values = np.array([361.400, 51.640, 15.957, 6.643, 3.230, 51.990, 15.793, 6.520, 3.191])
            k2k1_values = np.array([0, 0, 0, 0, 0, 0.000764, 0.00677, 0.0303, 0.1355])

            # Find array index for corresponding PQ-value
            index = np.where(digits == (10 * self.airfoil[1] + self.airfoil[2]))
            r = r_values[index]
            k1 = k1_values[index]
            k2k1 = k2k1_values[index]

            # Calculate mean camber line and linearly scale camber and gradient with the designed coefficient of lift
            Cl = 0.3
            if self.q == 0:
                y_c = self.m / Cl * np.where(x < r, k1 / 6 * (np.power(x, 3) - 3 * r * np.power(x, 2) + np.power(r, 2) * x * (3 - r)),
                               k1 * np.power(r, 3) / 6 * (1 - x))
                dy_c = self.m / Cl * np.where(x < r, k1 / 6 * (3 * np.power(x, 2) - 6 * r * x + (3 - r) * np.power(r, 2)),
                                -k1 * np.power(r, 3) / 6)
            else:
                y_c = self.m / Cl * np.where(x < r, k1 / 6 * (
                        np.power(x - r, 3) - k2k1 * x * np.power(1 - r, 3) - x * np.power(r, 3) + np.power(r, 3)),
                               k1 / 6 * (k2k1 * np.power(x - r, 3) - k2k1 * x * np.power(1 - r, 3) - x * np.power(r,
                                                                                                                  3) + np.power(
                                   r, 3)))
                dy_c = self.m / Cl * np.where(x < r, k1 / 6 * (3 * np.power(x - r, 2) - k2k1 * np.power(1 - r, 3) - np.power(r, 3)),
                                k1 / 6 * (3 * k2k1 * np.power(x - r, 2) - k2k1 * np.power(1 - r, 3) - np.power(r, 3)))

        theta = np.arctan(dy_c)

        # Translate aerodynamic centre (0.25c) to coordinate origin for rotation
        xUpper = x - y_t * np.sin(theta) - 0.25 * self.chord
        xLower = x + y_t * np.sin(theta) - 0.25 * self.chord

        yUpper = y_c + y_t * np.cos(theta)
        yLower = y_c - y_t * np.cos(theta)

        rotationMatrix = np.array([[np.cos(self.alpha), np.sin(self.alpha)], [-np.sin(self.alpha), np.cos(self.alpha)]])

        upper = np.einsum('ij, kj -> ki', rotationMatrix, np.transpose((xUpper, yUpper)))[:-1]
        lower = np.einsum('ij, kj -> ki', rotationMatrix, np.transpose((xLower, yLower)))[1:]

        vertices = np.concatenate((np.flip(lower, 0), upper))

        self.vertices = np.zeros((2 * vertices.shape[0], 3))
        self.vertices[:vertices.shape[0], 2] = 1
        self.vertices[vertices.shape[0]:, 2] = -1

        self.vertices[:vertices.shape[0], 0:2] = vertices
        self.vertices[vertices.shape[0]:, 0:2] = vertices

    def generateFaces(self):
        faces = np.zeros((self.vertices.shape[0], 3))

        for i in range(int(0.5 * faces.shape[0]) - 1):
            faces[2 * i] = [i, i + 1, i + int(0.5 * faces.shape[0])]
            faces[2 * i + 1] = [i + 1, i + 1 + int(0.5 * faces.shape[0]), i + int(0.5 * faces.shape[0])]

        faces[-2] = [int(0.5 * faces.shape[0]) - 1, 0, faces.shape[0] - 1]
        faces[-1] = [0, int(0.5 * faces.shape[0]), faces.shape[0] - 1]

        facesSide = np.zeros((self.vertices.shape[0] - 4, 3))

        for j in range(1, int(0.5 * facesSide.shape[0]) - 4):
            facesSide[2 * j - 1] = [j + 1, j, int(0.5 * self.vertices.shape[0]) - j]
            facesSide[2 * j] = [j + 1, int(0.5 * self.vertices.shape[0]) - j, int(0.5 * self.vertices.shape[0]) - j - 1]

        facesSide[0] = [1, 0, int(0.5 * self.vertices.shape[0]) - 1]
        facesSide[int(0.5 * facesSide.shape[0]) - 1] = [int(0.25 * self.vertices.shape[0]),
                                                        int(0.25 * self.vertices.shape[0]) - 1,
                                                        int(0.25 * self.vertices.shape[0]) + 1]

        facesSide[int(0.5 * facesSide.shape[0]):] = facesSide[:int(0.5 * facesSide.shape[0])] + int(
            0.5 * self.vertices.shape[0])

        faces = np.concatenate((faces, facesSide), axis=0)

        self.faces = faces.astype(int)

    def generateStandardTriangleLanguageFile(self):
        naca = mesh.Mesh(np.zeros(self.faces.shape[0], dtype=mesh.Mesh.dtype))
        for i, f in enumerate(self.faces):
            for j in range(3):
                naca.vectors[i][j] = self.vertices[f[j], :]

        folderDir = os.path.join(self.caseDir, 'constant/triSurface')

        if not os.path.exists(folderDir):
            os.makedirs(folderDir)

        stlFileName = "naca" + str(self.foil) + "_" + str(self.angle) + ".stl"
        saveDir = os.path.join(folderDir, stlFileName)

        naca.save(saveDir)
