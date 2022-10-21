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
import os
import pandas as pd
import numpy as np

cwd = os.getcwd()

databaseDir = os.path.join(cwd, "database")

databaseList = os.listdir(databaseDir)
if 'nacaFOAM_Cl_Cd_Cm.csv' in databaseList:
    databaseList.remove('nacaFOAM_Cl_Cd_Cm.csv')

df = pd.DataFrame(columns=['naca', 'angle', 'mach', 'Cl', 'Cd', 'Cm'])

i = 0
for case in databaseList:
    caseDir = os.path.join(databaseDir, case)
    caseName = case.split("_")
    for file in os.listdir(caseDir):
        if file.lower().endswith('.dat'):
            fileDir = os.path.join(caseDir, file)
            temp = pd.read_csv(fileDir, sep='\t', header=12, usecols=[1, 4, 7])
            temp = temp.iloc[-100:]
            temp = round(temp.mean(), 6)

            df.loc[i] = {'naca': caseName[0], 'angle': caseName[1], 'mach': caseName[2], 'Cl': temp[1], 'Cd': temp[0],
                         'Cm': temp[2]}

            i += 1

saveDir = os.path.join(databaseDir, "nacaFOAM_Cl_Cd_Cm.csv")
df.to_csv(saveDir, sep='\t')

print('Post processing completed! Lift-drag data saved into: {}'.format(saveDir))
