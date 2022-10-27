#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------
# Created By  : Sebastian Widmann
# Institution : TU Munich, Department of Aerospace and Geodesy
# Created Date: October 9, 2022
# version ='1.0'
# ---------------------------------------------------------------------------
"""

"""
# ---------------------------------------------------------------------------
import numpy as np
import subprocess
import os, psutil, shutil
import itertools
import multiprocessing

from datetime import datetime

from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
from PyFoam.Execution.UtilityRunner import UtilityRunner
from PyFoam.Execution.BasicRunner import BasicRunner
from PyFoam.Error import error
from PyFoam.Basics.Utilities import copytree
from PyFoam.Execution.ParallelExecution import LAMMachine

from scripts.generateStlFile import generateGeometryFile
from scripts.generateSurfaceFeatureExtractDict import generateSurfaceFeatureExtractDict
from scripts.generateBlockMeshDict import generateBlockMeshDict
from scripts.generateSnappyHexMeshDict import generateSnappyHexMeshDict
from scripts.generateFvOptions import generateFvOptionsDict
from scripts.generateForceCoefficientsDict import generateForceCoefficientsDict
from scripts.generateZeroDirectoryFiles import generateInitialConditionsFiles


def generateNacaAirfoils():
    '''
        Generate all required possibilities of NACA 4- and 5-series
        airfoils
    Returns
    -------
        list of combined combinations
    '''
    naca4, naca5 = [], []
    maxCamber = np.arange(0, 10, 1)
    position = np.arange(0, 10, 1)
    thickness = np.arange(5, 41, 5)

    for m in maxCamber:
        for p in position:
            for t in thickness:
                if (m != 0 and p == 0) or (m == 0 and p != 0):
                    continue
                naca4.append(str(m) + str(p) + str(t).zfill(2))

    maxCamber = np.arange(1, 7, 1)
    position = np.arange(2, 6, 1)
    reflex = np.arange(1, 2, 1)  # only include reflex airfoils
    thickness = np.arange(5, 31, 5)

    for m in maxCamber:
        for p in position:
            for q in reflex:
                for t in thickness:
                    if (m != 0 and p == 0) or (m == 0 and p != 0):
                        continue
                    naca5.append(str(m) + str(p) + str(q) + str(t).zfill(2))

    return naca4 + naca5


def printProgressBar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█', printEnd="\r"):
    '''
        Call in a loop to create terminal progress bar
    Parameters
    ----------
    iteration: current iteration [int]
    total: total iteration [int]
    prefix: prefix string [str]
    suffix: suffix string [str]
    decimals: positive number of decimals in percent complete [int]
    length: character length of bar [int]
    fill:bar fill character [str]
    printEnd: end character (e.g. "\r", "\r\n") [str]

    Returns
    -------

    '''
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end=printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()


def deleteTempMeshFiles(meshDir):
    '''
        Function will delete mesh related files which can cause
        snappyHexMesh to crash based on the special implementation
        how the mesh process is executed to reduce the runtime
    Parameters
    ----------
    meshDir: path to the polyMesh folder in each directory [str]

    Returns
    -------

    '''
    subprocess.run(["rm", "-f", os.path.join(meshDir, "cellLevel")])
    subprocess.run(["rm", "-f", os.path.join(meshDir, "level0Edge")])
    subprocess.run(["rm", "-f", os.path.join(meshDir, "pointLevel")])
    subprocess.run(["rm", "-f", os.path.join(meshDir, "surfaceIndex")])


def deletePyFoamTempFiles(newCase):
    '''
        Function will delete temporary files created by PyFoam which
        carry no relevant information of simulation to clear up storage
    Parameters
    ----------
    newCase: case name [str]

    Returns
    -------

    '''
    subprocess.run(['rm', '-f', 'PyFoamState.CurrentTime'])
    subprocess.run(['rm', '-f', 'PyFoamState.LastOutputSeen'])
    subprocess.run(['rm', '-f', 'PyFoamState.StartedAt'])
    subprocess.run(['rm', '-f', 'PyFoamState.TheState'])
    subprocess.run(['rm', '-f', '{}.foam'.format(newCase)])


def runCase(args):
    naca, angle, mach = args  # unpack arguments

    solver = "rhoSimpleFoam"
    baseCellSize = 0.2  # [m]

    cwd = os.getcwd()

    templateCase = SolutionDirectory("nacaFOAM-template", archive=None, paraviewLink=False)
    newCase = str(naca) + "_" + str(angle) + "_" + str(mach)
    case = templateCase.cloneCase(newCase)
    os.chdir(case.name)  # move to case directory to avoid OpenFOAM error of not finding controlDict

    print('Starting {}'.format(newCase))

    # ----- generate geometry -----
    generateGeometryFile(case.name, naca, angle)
    generateSurfaceFeatureExtractDict(case.name, naca, angle)

    surfaceFeatureExtract = BasicRunner(argv=['surfaceFeatureExtract', "-case", case.name], silent=True,
                                        logname="surfaceFeatureExtract")
    surfaceFeatureExtract.start()
    if surfaceFeatureExtract.runOK():
        subprocess.run(["rm", "-f", "surfaceFeatureExtract.logfile"])
    if not surfaceFeatureExtract.runOK():
        deletePyFoamTempFiles(newCase)
        os.chdir(cwd)
        return '{}_sFE'.format(newCase)

    meshDir = os.path.join(case.name, "constant/polyMesh")

    if baseCellSize != 0.2:
        # ----- blockMesh -----
        generateBlockMeshDict(case.name, baseCellSize)

        blockMesh = BasicRunner(argv=["blockMesh", "-case", case.name], silent=True, logname="blockMesh")
        blockMesh.start()
        if blockMesh.runOK():
            subprocess.run(["rm", "-f", "blockMesh.logfile"])
        if not blockMesh.runOK():
            deletePyFoamTempFiles(newCase)
            os.chdir(cwd)
            return '{}_bM'.format(newCase)

        # ----- snappyHexMesh -----
        generateSnappyHexMeshDict(case.name, naca, angle, mach, baseCellSize, 0)

        snappyHexMesh_0 = BasicRunner(argv=["snappyHexMesh", "-case", case.name, "-overwrite"], silent=True,
                                      logname="sHM_0")
        snappyHexMesh_0.start()
        if snappyHexMesh_0.runOK():
            subprocess.run(["rm", "-f", "sHM_0.logfile"])
        if not snappyHexMesh_0.runOK():
            deletePyFoamTempFiles(newCase)
            os.chdir(cwd)
            return '{}_sHM0'.format(newCase)

        deleteTempMeshFiles(meshDir)

    generateSnappyHexMeshDict(case.name, naca, angle, mach, baseCellSize, 1)

    snappyHexMesh_1 = BasicRunner(argv=["snappyHexMesh", "-case", case.name, "-overwrite"], silent=True,
                                  logname="sHM_1")
    snappyHexMesh_1.start()
    if snappyHexMesh_1.runOK():
        subprocess.run(["rm", "-f", "sHM_1.logfile"])
    if not snappyHexMesh_1.runOK():
        deletePyFoamTempFiles(newCase)
        os.chdir(cwd)
        return '{}_sHM1'.format(newCase)

    deleteTempMeshFiles(meshDir)

    generateSnappyHexMeshDict(case.name, naca, angle, mach, baseCellSize, 2)
    fvOptionsDir = os.path.join(case.name, "system/fvOptions")
    subprocess.run(["rm", "-f", fvOptionsDir])

    snappyHexMesh_2 = BasicRunner(argv=["snappyHexMesh", "-case", case.name, "-overwrite"], silent=True,
                                  logname="sHM_2")
    snappyHexMesh_2.start()
    if snappyHexMesh_2.runOK():
        subprocess.run(["rm", "-f", "sHM_2.logfile"])
    if not snappyHexMesh_2.runOK():
        deletePyFoamTempFiles(newCase)
        os.chdir(cwd)
        return '{}_sHM2'.format(newCase)

    # ----- extrudeMesh -----
    extrudeMesh = BasicRunner(argv=["extrudeMesh", "-case", case.name], silent=True, logname="extrudeMesh")
    extrudeMesh.start()
    if extrudeMesh.runOK():
        subprocess.run(["rm", "-f", "extrudeMesh.logfile"])
    if not extrudeMesh.runOK():
        deletePyFoamTempFiles(newCase)
        os.chdir(cwd)
        return '{}_eM'.format(newCase)

    # ----- createPatch -----
    createPatch = BasicRunner(argv=["createPatch", "-case", case.name, "-overwrite"], silent=True,
                              logname="createPatch")
    createPatch.start()
    if createPatch.runOK():
        subprocess.run(["rm", "-f", "createPatch.logfile"])
    if not createPatch.runOK():
        deletePyFoamTempFiles(newCase)
        os.chdir(cwd)
        return '{}_cP'.format(newCase)

    # ----- renumberMesh -----
    renumberMesh = BasicRunner(argv=["renumberMesh", "-case", case.name, "-overwrite"], silent=True,
                               logname="renumberMesh")
    renumberMesh.start()
    if renumberMesh.runOK():
        subprocess.run(["rm", "-f", "renumberMesh.logfile"])
    if not renumberMesh.runOK():
        deletePyFoamTempFiles(newCase)
        os.chdir(cwd)
        return '{}_rM'.format(newCase)

    # ----- checkMesh -----
    checkMesh = BasicRunner(argv=["checkMesh", "-case", case.name, "-latestTime"], silent=True, logname="checkMesh")
    checkMesh.start()
    if not checkMesh.runOK():
        deletePyFoamTempFiles(newCase)
        os.chdir(cwd)
        return '{}_cM'.format(newCase)

    # ----- generate various dicts -----
    generateFvOptionsDict(case.name)
    generateForceCoefficientsDict(case.name, mach)

    # ----- restore 0 directory from 0.orig -----
    generateInitialConditionsFiles(case.name, mach)
    copytree(os.path.join(case.name, "0.org"), os.path.join(case.name, "0"))

    # ----- decomposePar -----
    decomposePar = UtilityRunner(argv=["decomposePar"], silent=True, logname="decomposePar")
    decomposePar.start()
    if decomposePar.runOK():
        subprocess.run(["rm", "-f", "decomposePar.logfile"])
    if not decomposePar.runOK():
        deletePyFoamTempFiles(newCase)
        os.chdir(cwd)
        return '{}_dP'.format(newCase)

    # ----- rhoSimpleFoam -----
    machine = LAMMachine(nr=4)
    rhoSimpleFoam = BasicRunner(argv=[solver, "-case", case.name], silent=True, lam=machine, logname="rhoSimpleFoam")
    rhoSimpleFoam.start()
    if not rhoSimpleFoam.runOK():
        deletePyFoamTempFiles(newCase)
        os.chdir(cwd)
        return '{}_solver'.format(newCase)

    # ----- reconstructPar -----
    reconstructPar = UtilityRunner(argv=["reconstructPar"], silent=True, logname="reconstructPar")
    reconstructPar.start()
    if reconstructPar.runOK():
        subprocess.run(["rm", "-f", "reconstructPar.logfile"])
    if not reconstructPar.runOK():
        deletePyFoamTempFiles(newCase)
        os.chdir(cwd)
        return '{}_rP'.format(newCase)

    # ----- foamToVTK -----
    foamToVTK = BasicRunner(argv=["foamToVTK", "-latestTime", "-no-boundary", "-fields",
                                  "'(TMean UMean alphatMean kMean nutMean omegaMean pMean rhoMean)'", "-overwrite",
                                  "-name", newCase], silent=True, logname="foamToVTK")
    foamToVTK.start()
    if foamToVTK.runOK():
        subprocess.run(["rm", "-f", "foamToVTK.logfile"])
    if not foamToVTK.runOK():
        deletePyFoamTempFiles(newCase)
        os.chdir(cwd)
        return '{}_vtk'.format(newCase)

    # ----- data handling -----
    databaseDir = os.path.join(cwd, "database")
    sourceDir = os.path.join(case.name, newCase)
    targetDir = os.path.join(databaseDir, newCase)

    try:
        shutil.copy(os.path.join(case.name, "postProcessing/forceCoeffs/0/coefficient.dat"),
                    os.path.join(sourceDir, "{}_forceCoeffs.dat".format(newCase)))
    except FileNotFoundError:
        pass

    copytree(sourceDir, targetDir)  # copy VTK folder into database

    try:
        deletePyFoamTempFiles(newCase)
        subprocess.run(["rm", "-rf", "decomposePar.analyzed"])
        subprocess.run(["rm", "-rf", "reconstructPar.analyzed"])
    except FileNotFoundError:
        pass

    os.chdir(cwd)


def main():
    databaseDir = os.path.join(os.getcwd(), "database")
    if not os.path.exists(databaseDir):
        os.makedirs(databaseDir)

    angles = np.arange(-5, 21, 1)
    machs = np.arange(0.05, 0.65, 0.05)
    nacas = generateNacaAirfoils()
    paramlist = list(itertools.product(nacas, angles, np.round(machs, 9)))

    '''
        OpenFOAM doesn't support hyperthreading
        -> multiprocessing.cpu_count() returns logical core count
        -> psutil.cpu_count(args) will return physical core count
    '''
    nCores = psutil.cpu_count(logical=False) - 2
    nProc = 4
    nProcs = int(nCores / nProc)

    startTime = datetime.now()

    with multiprocessing.Pool(processes=nProcs) as pool:
        with open('errors.txt', 'w') as f:
            for msg in pool.imap_unordered(runCase, paramlist):
                if msg is not None:
                    print(msg, file=f)

    print(datetime.now() - startTime)


if __name__ == '__main__':
    main()
