""" drag polar csv to yaml """
import os
import yaml
import pandas as pd
import numpy as np

# pylint: disable=invalid-name
# pylint: disable=unspecified-encoding
# pylint: disable=too-many-locals
# pylint: disable=line-too-long


# df=pd.read_csv('data/rio_5-15_0m_s.csv', skipinitialspace=True)
# d = df.to_dict( orient='list')
# print(np.asarray(d['alpha']))

def csv_to_yaml(file_in):
    """ create yaml from csv file """
    name=os.path.basename(file_in)
    file_out=os.path.join(yaml_dir, name+'.yaml')

    df=pd.read_csv(file_in, skipinitialspace=True)
    data_dict=df.to_dict(orient='list')

    data={
        'alpha': data_dict['alpha'],
        'cL': data_dict['CL'],
        'cD': data_dict['CD'],
        'cDi': data_dict['CDi']
    }

    with open(file_out, 'w') as path:
        yaml.safe_dump(data, path, sort_keys=False)

data_dir = 'data'
yaml_dir = 'yaml'

for file in os.listdir(data_dir):
    # parse through each file in data directory
    _, ext = os.path.splitext(file)
    if ext != '.csv':
        # Ignore hidden files and non-txt files
        continue

    path=os.path.join(data_dir, file)
    print(f'processing file: {path}')
    csv_to_yaml(file_in=path)
