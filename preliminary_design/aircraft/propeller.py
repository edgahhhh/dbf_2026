""" propeller """

import yaml
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# pylint: disable=unspecified-encoding
# pylint: disable=invalid-name

# Integrate into a revamped aircraft model at some point
class Propeller():
    """ Propeller model"""
    def __init__(self, propeller_info):
        with open(propeller_info, 'r') as file:
            self.prop_info = yaml.safe_load(file)
        self.brand = self.prop_info['brand']    # string
        self.diameter = self.prop_info['geometry']['diameter']  # m
        self.pitch = self.prop_info['geometry']['pitch']    # m
        self.j_data = self.prop_info['performance_data']['J']   # dimensionless (1/rev)
        self.ct_data = self.prop_info['performance_data']['CT'] # dimensionless
        self.cp_data = self.prop_info['performance_data']['CP'] # dimensionless
        self.eta_data = self.prop_info['performance_data']['eta']   # dimensionless

        self.rho = self.density_from_alt(400)   # density at 400m altitude

        # Interpolators for data
        self.interp_ct = interp1d(self.j_data, self.ct_data, fill_value='extrapolate')
        self.interp_cp = interp1d(self.j_data, self.cp_data, fill_value='extrapolate')
        self.interp_eta = interp1d(self.j_data, self.eta_data, fill_value='extrapolate')

    def density_from_alt(self, height_asl):
        """ calculate density from altitude 
        @param height_asl: height above sea level [m]
        """
        theta = 1 + (-0.000022558) * height_asl
        rho = 1.225 * (theta**4.2561)
        return rho

    def get_advance_ratio(self, V, n):
        """ Calculate advance ratio
        @param V: velocity [m/s]
        @param n: rotation rate [rad/s]
        function converts rad/sec to rev/sec, as the propeller data uses rev/sec
        """
        J = V / n / self.diameter * 2*np.pi
        return J

    def calculate_thrust(self, V, n):
        """ Calculate thrust
        @param V: velocity [m/s]
        @param n: rotation rate [rad/s]
        """
        J = self.get_advance_ratio(V, n)
        ct = max(min(self.ct_data), min(self.interp_ct(J), max(self.ct_data)))
        D = self.diameter
        T = ct * self.rho * (n/2/np.pi)**2 * D**4
        return T

    def calculate_power(self, V, n):
        """ Calculate power
        @param V: velocity [m/s]
        @param n: rotation rate [rad/s]
        """
        J = self.get_advance_ratio(V, n)
        cp = max(min(self.cp_data), min(self.interp_cp(J), max(self.cp_data)))
        D = self.diameter
        P = cp * self.rho * (n/2/np.pi)**3 * D**5
        return P

    def calculate_efficiency(self, V, n):
        """ Calculate efficiency
        @param V: velocity [m/s]
        @param n: rotation rate [rad/s]
        """
        J = self.get_advance_ratio(V, n)
        eta = max(min(self.eta_data), min(self.interp_eta(J), max(self.eta_data)))
        return eta

    def plot_original_data(self):
        """ Plot the original propeller data """
        plt.figure()
        plt.subplot(3,1,1)
        plt.plot(self.j_data, self.ct_data, 'o-')
        plt.xlabel('J')
        plt.ylabel('CT')

        plt.subplot(3,1,2)
        plt.plot(self.j_data, self.cp_data, 'o-')
        plt.xlabel('J')
        plt.ylabel('CP')

        plt.subplot(3,1,3)
        plt.plot(self.j_data, self.eta_data, 'o-')
        plt.xlabel('J')
        plt.ylabel('eta')

        plt.tight_layout()
        plt.show()
