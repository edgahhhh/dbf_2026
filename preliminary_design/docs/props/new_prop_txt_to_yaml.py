""" new text to yaml 
    run from terminal at proper path,  cd props then python3 this script """
import os
import yaml
import pandas as pd
import numpy as np

# pylint: disable=invalid-name
# pylint: disable=unspecified-encoding
# pylint: disable=too-many-locals
# pylint: diable=line-too-long

def dissect_name(file_path):
    """ Dissect the file name into the proper parts 
    Returns:

        brand: string
        rpm: string
        diameter: string
        pitch: string
    """
    name = os.path.basename(file_path)
    base, _ = os.path.splitext(name)

    parts = base.split('_')
    if len(parts)!=4:
        raise ValueError(f"Unexpected file name format: {name}")
    
    diam_pitch_parts = parts[1].split('x')
    return parts[0], parts[3], diam_pitch_parts[0], diam_pitch_parts[1]

def text_to_yaml(file_one, file_two):
    """ Create yaml from both propellers .txt files
    Args:
        file_one: string
        file_two: string
    """
    name_one = os.path.basename(file_one)
    if file_two is None:
        name_two = None
    else:
        name_two = os.path.basename(file_two)

    prop, rpm, diam, pit = dissect_name(name_one)

    if float(diam) >= 100:
        # Divide by 10 since some decimals get missed
        diam=str(float(diam)/10)

    file_out = os.path.join(yaml_dir, prop + '_' + diam + 'x' + pit + '_' + rpm + '.yaml')

    df_one = pd.read_csv(file_one, sep=r'\s+')
    orig_performance_data_one = df_one.to_dict(orient='list')
    if name_two is not None:
        df_two = pd.read_csv(file_two, sep=r'\s+')
        orig_performance_data_two = df_two.to_dict(orient='list')
        orig_J = np.concatenate(
            (np.asarray(orig_performance_data_one['J']),np.asarray(orig_performance_data_two['J'])),axis=0)
        orig_CT = np.concatenate(
            (np.asarray(orig_performance_data_one['CT']),np.asarray(orig_performance_data_two['CT'])),axis=0)
        orig_CP = np.concatenate(
            (np.asarray(orig_performance_data_one['CP']),np.asarray(orig_performance_data_two['CP'])),axis=0)
        orig_eta = np.concatenate(
            (np.asarray(orig_performance_data_one['eta']),np.asarray(orig_performance_data_two['eta'])),axis=0)

    else:
        orig_J = np.asarray(orig_performance_data_one['J'])
        orig_CT = np.asarray(orig_performance_data_one['CT'])
        orig_CP = np.asarray(orig_performance_data_one['CP'])
        orig_eta = np.asarray(orig_performance_data_one['eta'])

    # Convert to use n in rad/sec (disabled)
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

data_dir = 'new_data'
yaml_dir = 'new_yaml'

brand_check = []
diameter_check = []
pitch_check = []
processed = set()
for file in os.listdir(data_dir):
    # parse through each file in data directory
    _, ext = os.path.splitext(file)
    if ext != '.txt':
        # Ignore hidden files and non-txt files
        continue
    brand, rpm, diameter, pitch = dissect_name(file)
    combo = (brand.lower(), diameter.lower(), pitch.lower())

    if combo in processed:
        continue

    # if (
    #     brand in brand_check
    #     and diameter in diameter_check
    #     and pitch in pitch_check
    # ):
    #     print(brand, diameter, pitch)
    #     continue

    for other_file in os.listdir(data_dir):
        # Look for the other file that goes with current one
        _, other_ext = os.path.splitext(other_file)
        if other_ext != '.txt':
            # Ignore hidden files and non-txt files
            continue
        other_brand, other_rpm, other_diameter, other_pitch = dissect_name(other_file)

        # print(brand, other_brand)
        if (
            other_brand.lower() == brand.lower()
            and other_diameter.lower() == diameter.lower()
            and other_pitch.lower() == pitch.lower()
            and rpm.lower() != other_rpm.lower()
        ):
            file_pair = other_file
            break

    # Check to see if pair was found or not
    # Checking in loop was being weird
    if (
        other_brand.lower() != brand.lower()
        or other_diameter.lower() != diameter.lower()
        or other_pitch.lower() != pitch.lower()
    ):
        file_pair = None

    path_one = os.path.join(data_dir, file)
    if file_pair is not None:
        path_two = os.path.join(data_dir, file_pair)
    else:
        path_two = None
    print(f'processing files: {[path_one, path_two]}')
    text_to_yaml(
        file_one=path_one,
        file_two=path_two
        )
    # Add to set of complete propellers
    processed.add(combo)
