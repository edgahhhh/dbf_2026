""" Propeller Analysis """
# pylint: disable=invalid-name
# pylint: disable=too-many-public-methods
# pylint: disable=unused-import
import os
from tqdm import tqdm
import numpy as np
import yaml
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import Axes3D
from aircraft.propeller import Propeller
from aircraft.simple_plane import Aircraft
from propeller_analysis import PropellerAnalysis
plt.rcParams['axes.grid'] = True

