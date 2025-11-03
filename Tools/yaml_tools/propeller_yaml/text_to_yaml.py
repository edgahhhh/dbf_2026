""" text to yaml """
import os
import yaml
import pandas as pd
# import numpy as np

# pylint: disable=invalid-name

def main(fileIn, fileOut, D, p, N):
    """
    Creating a dictionary and dumping it in some yaml file

    Inputs:
    -------
    fileIn : string
        Location of .txt file to read from 
    fileOut : string
        path to file location relative to where this is run form
    D : float
        diameter in m
    p : float
        pitch in m
    N : float
        RPM that test data was run at

    """

    df = pd.read_csv(fileIn, sep=r'\s+')
    performance_data = df.to_dict(orient='list')
    
    geometry_data = {
        'diameter': D,
        'pitch': p,
        'units': 'm'
    }
    run_data = [
        {'rpm': N}
    ]
    all_data = {
        'geometry': geometry_data,
        'performance_data': performance_data,
        'run_data': run_data
    }

    with open(fileOut, 'w') as file:
        yaml.dump(all_data, file, sort_keys=False)


# Future work could include making these easier to type out as its prone to mistakes

trunkPath = os.getcwd()
dataInPath = 'tools/yaml_automation/data/prop'
dataOutPath = 'simulation'

filesIn = [  'apcff_4.2x4_0617rd_6043.txt',
             'apcff_4.2x4_0618rd_8043.txt',
             'apcff_4.2x4_0619rd_8067.txt',
             'apcff_4.2x4_0620rd_10042.txt',
             'apcff_4.2x4_0621rd_10071.txt']
filesOut = ['4.2x4_6043.yaml',
            '4.2x4_8043.yaml',
            '4.2x4_8067.yaml',
            '4.2x4_10042.yaml',
            '4.2x4_10071.yaml']

filesInPath = [os.path.join(trunkPath, dataInPath, filename) for filename in filesIn]
filesOutPath = [os.path.join(trunkPath, dataOutPath, filename) for filename in filesOut]

print(f' \n\n {filesInPath} \n\n')
diameter = 0.107
pitch = 0.102
rpm = [6043,
       8043,
       8067,
       10042,
       10071]

for a, b, c in zip(filesInPath, filesOutPath, rpm):
    main(a, b, diameter, pitch, c)
