""" propeller """
import warnings
import yaml
import numpy as np
from typing import Optional, Tuple, NamedTuple
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# pylint: disable=unspecified-encoding
# pylint: disable=invalid-name

class PropellerAnalysis(NamedTuple):
    """ Propeller analysis class """
    thrust: float
    power: float
    eta: float
    omega: float

class Propeller():
    """ Propeller class to model propellers using a yaml data file"""
    def __init__(self, propeller_info, altitude=400.0):
        """ load in propeller data from yaml file 
            @param propeller_info: path to propeller yaml file
            @param altitude: altitude above sea level [m]
        """
        with open(propeller_info, 'r') as file:
            self._prop_info = yaml.safe_load(file)

        # TODO: protect some of these variables that are sensitive
        self.brand = self._prop_info['brand']    # string
        self.diameter = self._prop_info['geometry']['diameter']  # m
        self.pitch = self._prop_info['geometry']['pitch']    # m
        self._j_data = self._prop_info['performance_data']['J']   # dimensionless (1/rev)
        self._ct_data = self._prop_info['performance_data']['CT'] # dimensionless
        self._cp_data = self._prop_info['performance_data']['CP'] # dimensionless
        self._eta_data = self._prop_info['performance_data']['eta']   # dimensionless

        self._rho = self.density_from_alt(altitude)   # density at 400m altitude

        # Interpolators for data
        self._interp_ct = interp1d(self._j_data, self._ct_data, fill_value='extrapolate')
        self._interp_cp = interp1d(self._j_data, self._cp_data, fill_value='extrapolate')
        self._interp_eta = interp1d(self._j_data, self._eta_data, fill_value='extrapolate')

        self._name = f"{self.brand}_{self.diameter*39.37:.0f}x{self.pitch*39.37:.0f}"

    @property
    def name(self):
        """ get name of propeller """
        return self._name

    @staticmethod
    def density_from_alt(height_asl):
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
            @return J: advance ratio [dimensionless, ref. revs]
        """
        J = V / (n/2/np.pi * self.diameter)
        return J

    def calculate_thrust(self, V, n):
        """ Calculate thrust with saturation
            @param V: velocity [m/s]
            @param n: rotation rate [rad/s]
            @return T: thrust force [N]
        """
        J = self.get_advance_ratio(V, n)
        ct = max(min(self._ct_data), min(self._interp_ct(J), max(self._ct_data)))
        D = self.diameter
        T = ct * self._rho * (n/2/np.pi)**(2) * D**(4)
        return T

    def calculate_power(self, V, n):
        """ Calculate power with saturation
            @param V: velocity [m/s]
            @param n: rotation rate [rad/s]
            @return P: power [W]
        """
        J = self.get_advance_ratio(V, n)
        cp = max(min(self._cp_data), min(self._interp_cp(J), max(self._cp_data)))
        D = self.diameter
        P = cp * self._rho * (n/2/np.pi)**(3) * D**(5)
        return P

    def calculate_efficiency(self, V, n):
        """ Calculate efficiency with saturation
            @param V: velocity [m/s]
            @param n: rotation rate [rad/s]
            @return eta: efficiency [dimensionless]
        """
        J = self.get_advance_ratio(V, n)
        eta = max(min(self._eta_data), min(self._interp_eta(J), max(self._eta_data)))
        return eta

    def plot_original_data(self):
        """ Plot the original propeller data """
        plt.figure()
        plt.subplot(3,1,1)
        plt.plot(self._j_data, self._ct_data, 'o-')
        plt.xlabel('J')
        plt.ylabel('CT')

        plt.subplot(3,1,2)
        plt.plot(self._j_data, self._cp_data, 'o-')
        plt.xlabel('J')
        plt.ylabel('CP')

        plt.subplot(3,1,3)
        plt.plot(self._j_data, self._eta_data, 'o-')
        plt.xlabel('J')
        plt.ylabel('eta')

        plt.title(f'Original data for {[self.brand, self.diameter, self.pitch]}')
        plt.tight_layout()
        plt.show()

    def find_propeller_speed(self, airspeed:float, drag:float):
        """ find propeller speed to meet thrust requirement at airspeed 
            @param airspeed: velocity [m/s]
            @param drag: required thrust [N]
            @return omega: propeller speed [rad/sec] 
        """
        a = 1
        b = 100000
        fa = self.calculate_thrust(airspeed,n=a) - drag
        for i in range(100):
            c = (a+b)/2 # rad/sec
            f = self.calculate_thrust(airspeed,n=c) - drag
            if abs(f) <= 0.01:
                # print(f'object: {self.label}, converged with find_rpm() after {i} iterations')
                return c
            elif fa * f < 0:
                b = c
            else:
                a = c
                fa = f
        warnings.warn('solution did not converge', UserWarning, stacklevel=2)
        return None

    def analysis(
                self, airspeed:float, drag:float
                 ) -> Optional[PropellerAnalysis]:
        """perform analysis for propeller

        Arguments:
            airspeed [m/s]
            drag [N]
        Returns:
            thrust [N]
            power [W]
            eta [1]
            omega [rad/sec]
        """
        omega = self.find_propeller_speed(airspeed, drag=drag)  # rad/sec
        if omega is None:
            return None
        thrust=self.calculate_thrust(V=airspeed, n=omega)
        power=self.calculate_power(V=airspeed, n=omega)
        eta = self.calculate_efficiency(V=airspeed, n=omega)
        return PropellerAnalysis(thrust, power, eta, omega)
    
# TODO: change doc strings from @param and @return to standard Arguments: and Returns:
