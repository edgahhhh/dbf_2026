""" Motor, Battery, and Propeller Analysis """
# pylint: disable=invalid-name
# pylint: disable=too-many-public-methods
# pylint: disable=unused-import
import os
from tqdm import tqdm
import numpy as np
# import yaml
import matplotlib.pyplot as plt
# from scipy.interpolate import interp1d
# from mpl_toolkits.mplot3d import Axes3D
# from aircraft.propeller import Propeller
from aircraft.simple_plane import Aircraft
from propeller_analysis import PropellerAnalysis
plt.rcParams['axes.grid'] = True

class Battery():
    """ Battery Model """
    def __init__(self, series, parallel, esr):
        """ initialize battery

        Not accounting for capacitance or anything dynamic

        Args:
            series: number of cells in series
            parallel: number of cells in parallel
            esr: equivalent series resistance of battery s
        """
        self.label = f'{series}S{parallel}P'
        self.info = set(series, parallel, esr)

        self.V_nominal = series * 3.75  # Volts
        self.esr = esr

        self.power_loss = 0 # W

    def calc_power_loss(self, power_out):
        """ Find power loss due to esr at some condition """
        current_out = power_out / self.V_nominal
        return current_out**(2)/self.esr

class Motor():
    """ Simple Motor """
    def __init__(self, back_emf, zero_load_current, terminal_resistance):
        """ Initialize motor 
        Args:
            back_emf: motor constant [rpm/Volt]
            zero_load_current: motor constant [Amperes]
            terminal_resistance: motor constant [Ohm]
        """
        self.kv = back_emf
        self.i_no = zero_load_current
        self.r_no = terminal_resistance
        self.eff = 0.9

def motor_selection(battery_model:Battery, prop_rpm ):
    """ Find motor constants for battery to match propeller """
    # return motor


# Find what rpm is needed to get some thrust at some airspeed
# Compare all the shaft powers for each propeller
# Shaft power needed directly correlates to energy consumption
yaml_dir = 'preliminary_design/docs/props/new_yaml'

conditions = {
    'altitude asl':397,
    'cruise altitude agl':10
    }
aircraft_design_parameters = {
    'aspect ratio' : 3.644246,
    'wing span':1.524,
    'CLmax':0.85,
    'CDmin_no_fuselage':0.04,
    'fineness ratio':6,
    'm2 bank angle':0.96607,
    'm3 bank angle': 0.5163
    }
aircraft_design_parameters['conditions'] = conditions
mission_parameters = {
    'cargo':1,
    'ducks':3,
    'mission two cruise speed': 16.587256,
    'banner length':8.3,
    'mission three cruise speed':12.510729,
    'banner aspect ratio':5
    }
course_parameters = {
    'ground run': 20,
    'time limit': 300,
    'length of straights': 304.8
    }
initial_guess = {
    'motor mass': 0.1,
    'm2 battery mass': 0.8,
    'm3 battery mass':0.5
    }

count = 0
for file in tqdm(os.listdir(yaml_dir), desc='Processing propeller files...'):
    path = os.path.join(yaml_dir, file)
    prop = PropellerAnalysis(path)
    # Solve for best banner
    l, v_cruise = prop.size_banner(
        lower_bound=2,
        upper_bound=20,
        iterations=20
        )
    # Create aircraft model for banner
    # Banner is sized now we have to see 
    if count == 0:
        count += 1
        break
