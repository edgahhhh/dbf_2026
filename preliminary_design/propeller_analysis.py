""" Propeller Analysis """
# pylint: disable=invalid-name
# pylint: disable=too-many-public-methods
# pylint: disable=unused-import
import numpy as np
import matplotlib.pyplot as plt
from aircraft.propeller import Propeller
from aircraft.simple_plane import Aircraft
plt.rcParams['axes.grid'] = True

class PropellerAnalysis():
    """ object for finding adequate props """
    def __init__(self, propeller_info,banner_path='preliminary_design/docs/banner_data.yaml'):
        """ Initialize object for each propeller """
        self.propeller = Propeller(propeller_info)
        self.label = f"{self.propeller.brand}_{self.propeller.diameter*39.37:.0f}x{self.propeller.pitch*39.37:.0f}"

        # Initialize an aircraft object for placeholders and to use its public functions
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
        self.dummy_plane = Aircraft(
            aircraft_design_parameters,
            mission_parameters,
            course_parameters,
            initial_guess,
            banner_data=banner_path)
        self.dummy_plane.size_aircraft_all_missions()

    def find_propeller_speed(self, airspeed:float, drag:float):
        """ find propeller speed to meet thrust requirement at airspeed 
            @param airspeed: velocity [m/s]
            @param drag: required thrust [N]
            @return omega: propeller speed [rad/sec] 
        """
        a = 1
        b = 100000
        fa = self.propeller.calculate_thrust(airspeed,n=a) - drag
        for i in range(100):
            c = (a+b)/2 # rad/sec
            f = self.propeller.calculate_thrust(airspeed,n=c) - drag
            if abs(f) <= 0.01:
                # print(f'object: {self.label}, converged with find_rpm() after {i} iterations')
                return c
            elif fa * f < 0:
                b = c
            else:
                a = c
                fa = f
        print(f'object: {self.label} did not converge with find_rpm() after {i} iterations...')
        return 1    # debugging

    def analysis(self, airspeed:float, drag:float):
        """ perform analysis for propeller 
            @param airspeed [m/s]
            @param drag [N]

            @return thrust [N]
            @return power [W]
            @return eta [1]
            @return omega [rad/sec]
        """
        omega = self.find_propeller_speed(airspeed, drag=drag)  # rad/sec
        thrust = self.propeller.calculate_thrust(V=airspeed, n=omega)
        power = self.propeller.calculate_power(V=airspeed, n=omega)
        eta = self.propeller.calculate_efficiency(V=airspeed, n=omega)
        return thrust, power, eta, omega

    def energy_from_banner(self, banner):
        """ Calculate energy for a banner at given airspeed
        function uses the aircraft object solved at __init__()
            @param banner: length of banner [m]
            
            @return energy_total: total energy needed for given condition [Whr]
            @return cruise_speed: cruising speed [m/s]
            @return extras: set of extra useful information
                [thrust_cruise,
                power_cruise,
                eta_cruise,
                omega_cruise,
                thrust_turn,
                power_turn,
                eta_turn,
                omega_turn,
                turn_drag,
                tof_speed]
            """
        energy_tof = 0.55264     # Whr
        energy_climb = 0.3668    # Whr
        t_tof = 3.4         # s
        t_climb = 1.385     # s
        mass_no_banner = self.dummy_plane.mass_m3_gross-self.dummy_plane.mass_banner
        area_banner = banner**(2)/5 # m^2
        mass_banner = area_banner * 0.0359  # kg, linear trend w/ area
        mass_total = mass_no_banner + mass_banner
        cruise_airspeed = 1.2*np.sqrt(mass_total*9.81 / (
            1/2*self.dummy_plane.rho*self.dummy_plane.S*self.dummy_plane.CLmax))
        tof_speed= 1.1*np.sqrt(mass_total*9.81 / (
            1/2*self.dummy_plane.rho*self.dummy_plane.S*self.dummy_plane.CLmax))
        
        # find cl at cruise current config
        cl_cruise = mass_total*9.81 / (
            (1/2*self.dummy_plane.rho*cruise_airspeed**(2)*self.dummy_plane.S))

        # find drag of banner at cruise and turn
        temp_cd = self.dummy_plane.banner_cd_interp_fun(banner)
        cd_banner = max( min(self.dummy_plane.banner_data_cd), 
                        min(temp_cd, max(self.dummy_plane.banner_data_cd)))*(
                            area_banner/self.dummy_plane.S)

        cd_total_cruise = self.dummy_plane.CD_m3_cruise - self.dummy_plane.CD_m3_banner_cruise  + cd_banner + self.dummy_plane.k*cl_cruise**(2)
        cd_total_turn = self.dummy_plane.CD_m3_cruise - self.dummy_plane.CD_m3_banner_cruise  + cd_banner + self.dummy_plane.k*(self.dummy_plane.n_m3*cl_cruise)**(2)

        cruise_drag = cd_total_cruise * 1/2 * self.dummy_plane.rho * cruise_airspeed**(2) * self.dummy_plane.S   # N
        turn_drag = cd_total_turn * 1/2 * self.dummy_plane.rho * cruise_airspeed**(2) * self.dummy_plane.S    # N

        # find shaft power for cruising and turning
        thrust_cruise, power_cruise, eta_cruise, omega_cruise = self.analysis(
            airspeed=cruise_airspeed,
            drag=cruise_drag
            )
        thrust_turn, power_turn, eta_turn, omega_turn = self.analysis(
            airspeed=cruise_airspeed,
            drag=turn_drag
            )
        if omega_cruise == 1 or omega_turn == 1:
            return 1,1,1
        power_cruise_batt = power_cruise * 0.85
        power_turn_batt = power_turn * 0.85
        _, _, energy_cruise = self.dummy_plane.calculate_cruise_performance(
            level_power_batt=power_cruise_batt,
            turn_power_batt=power_turn_batt,
            mass=0,
            cruise_speed=cruise_airspeed,
            time_tof=t_tof,
            time_climb=t_climb,
            phi=self.dummy_plane.phi_m3
        )
        energy_total = (energy_cruise + energy_tof + energy_climb)/0.75
        extras=[thrust_cruise,power_cruise,eta_cruise,omega_cruise,
                thrust_turn,power_turn,eta_turn,omega_turn,
                turn_drag, tof_speed]
        return energy_total, cruise_airspeed, extras

    def size_banner(self, lower_bound=1, upper_bound=20, iterations=20):
        """ Estimate the size of banner capable of being towed 
        Description:
            banner length was previously explicitly defined, now we can try to solve
            for a banner size that is best for some given airspeed
        Args:
        
            lower_bound: lower banner length bracket [m]

            upper_bound: upper banner length bracket [m]

        Returns:
            banner_length: length of banner [m]

            cruise_speed: cruising airspeed [m/s]

            extras: set of other useful data
                [thrust_cruise,
                power_cruise,
                eta_cruise,
                omega_cruise,
                thrust_turn,
                power_turn,
                eta_turn,
                omega_turn,
                turn_drag,
                tof_speed]
        """
        target_energy = 99.9
        a = lower_bound
        b = upper_bound

        fa, _, _ = self.energy_from_banner(a)
        fa = fa - target_energy
        for i in range(iterations):
            c = (a+b)/2
            f, v_cruise, extras = self.energy_from_banner(c)
            f = f - target_energy

            if abs(f) <= 0.1:
                # print(f'object: {self.label}, converged with size_banner() after {i} iterations')
                return c, v_cruise, extras
            elif f == -98:
                break
            
            if fa * f < 0:
                b = c
            else:
                a = c
                fa = f
        print(f'object: {self.label}, did not converge with size_banner() after {i} iterations...')
        return 1, 1, 1
