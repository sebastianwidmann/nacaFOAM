# ----------------------------------------------------------------------------
# Created By  : Sebastian Widmann
# Institution : TU Munich, Department of Aerospace and Geodesy
# Created Date: October 9, 2022
# version ='1.0'
# ---------------------------------------------------------------------------
import os
import sys
import pandas as pd
import argparse
from tqdm import tqdm


def forceCoeffs_pp(src_dir: str):
    if os.path.isfile(os.path.join(
            os.path.dirname(src_dir), 'airfoilMNIST_forceCoeffs.csv')):
        proceed = input('forceCoeffs data file exists. Do you want to '
                        'overwrite? ')
        if proceed.upper() == 'Y' or proceed.upper() == "YES":
            pass
        elif proceed.upper() == "N" or proceed.upper() == "NO":
            sys.exit('Program aborted! File will not be overwritten.')

    dataset_list = os.listdir(src_dir)

    df = pd.DataFrame(
        columns=['naca', 'angle', 'mach', 'Cl', 'Cd', 'Cm'])

    i = 0
    for sim in tqdm(dataset_list):
        sim_config = sim[:-4]  # remove file ending
        naca, angle, mach = sim_config.split("_", 2)[:3]

        file_dir = os.path.join(src_dir, sim)

        data = pd.read_csv(file_dir, sep='\t', header=12,
                           usecols=[1, 4, 7])
        data = round(data.iloc[-100:].mean(), 6)

        df.loc[i] = {'naca': naca,
                     'angle': angle,
                     'mach': mach,
                     'Cl': data[1],
                     'Cd': data[0],
                     'Cm': data[2],
                     }

        i += 1

    save_dir = os.path.join(os.path.dirname(src_dir),
                            'airfoilMNIST_forceCoeffs.csv')

    df.to_csv(save_dir, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='airfoilMNIST Postprocessing',
        description='Generate .csv-file with mean force and moment coeffs'
    )

    parser.add_argument(
        'src_dir',
        type=str,
        help='Specify load directory.'
    )

    args = parser.parse_args()

    forceCoeffs_pp(args.src_dir)
