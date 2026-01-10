""" text to yaml 
    run from terminal at proper path,  cd props then python3 this script """
import os
import yaml
import pandas as pd
import numpy as np

# pylint: disable=invalid-name
# pylint: disable=unspecified-encoding
# pylint: disable=too-many-locals

def text_to_yaml(file_in):
    """ Create yaml from propeller .txt file
    @param file_in [string]
    """
    file_name = os.path.basename(file_in)
    base, _ = os.path.splitext(file_name)

    # Parse base to get the brand, diameter, pitch, and rpm
    # brand_diameterxpitch_text_rpm.txt
    parts = base.split('_')
    if len(parts)!=4:
        raise ValueError(f"Unexpected file name format: {file_name}")
    brand = parts[0]
    diameter_pitch = parts[1]
    rpm = parts[3]
    diam_pitch_parts = diameter_pitch.split('x')
    diam = diam_pitch_parts[0]
    pitch = diam_pitch_parts[1]

    file_out = os.path.join(yaml_dir, brand + '_' + diam + 'x' + pitch + '_' + rpm + '.yaml')

    df = pd.read_csv(file_in, sep=r'\s+')
    orig_performance_data = df.to_dict(orient='list')

    # Convert from imperical to metric
    orig_J = np.asarray(orig_performance_data['J'])
    orig_CT = np.asarray(orig_performance_data['CT'])
    orig_CP = np.asarray(orig_performance_data['CP'])
    orig_eta = np.asarray(orig_performance_data['eta'])

    # Convert to use n in rad/sec (diabled)
    new_J = orig_J #* (2*np.pi)
    new_CT = orig_CT #* (2*np.pi)**(2)
    new_CP = orig_CP #* (2*np.pi)**(3)
    new_eta = orig_eta

    performance_data = {
        'J': new_J.tolist(),
        'CT': new_CT.tolist(),
        'CP': new_CP.tolist(),
        'eta': new_eta.tolist()
    }

    geometry_data = {
        'diameter': float(diam)/39.37,
        'pitch': float(pitch)/39.37,
        'units': 'm, rev, s'
    }
    run_data = [
        {'rpm': float(rpm)}
    ]
    all_data = {
        'brand': brand,
        'geometry': geometry_data,
        'performance_data': performance_data,
        'run_data': run_data
    }

    with open(file_out, 'w') as path:
        yaml.dump(all_data, path, sort_keys=False)

data_dir = 'data'
yaml_dir = 'yaml'

for file in os.listdir(data_dir):
    # parse through each file in data directory
    _, ext = os.path.splitext(file)
    if ext != '.txt':
        # Ignore hidden files and non-txt files
        continue
    path_in = os.path.join(data_dir, file)
    text_to_yaml(file_in=path_in)
