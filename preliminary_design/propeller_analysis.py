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
plt.rcParams['axes.grid'] = True

class PropellerAnalysis():
    """ object for finding adequate props """
    def __init__(self, propeller_info):
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
            banner_data='preliminary_design/docs/banner_data.yaml')
        self.dummy_plane.size_aircraft_all_missions()

    def find_rpm(self, airspeed:float, drag:float):
        """ find rpm to meet thrust requirement at airspeed 
        brute force solver, could use a more omtimized version for convergence
            @param airspeed: velocity [m/s]
            @param drag: required thrust [N]
        """
        rpm_array = np.linspace(100, 60000, 10000)
        n_array = rpm_array * 2*np.pi/60  # rad/s

        threshold = 0.1  # N tolerance
        for n in n_array:
            # Step through speeds to find required frequency
            thrust = self.propeller.calculate_thrust(V=airspeed, n=n)
            if np.isclose(thrust, drag, atol=threshold):
                return n*60/2/np.pi
        print(f'object: {self.label} did not converge with find_rpm()...')
        return 1    # debug value
        
    def analysis(self, airspeed:float, drag:float):
        """ perform analysis for propeller 
            @param airspeed [m/s]
            @param drag [N]

            @return thrust [N]
            @return power [W]
            @return eta [1]
            @return rpm [rev/sec]
        """
        rpm = self.find_rpm(airspeed, drag=drag)
        n = rpm * 2*np.pi/60  # rad/s
        thrust = self.propeller.calculate_thrust(airspeed, n)
        power = self.propeller.calculate_power(airspeed, n)
        # eta = self.propeller.calculate_efficiency(airspeed, n)
        eta = airspeed*thrust / power
        return thrust, power, eta, rpm

    def energy_from_banner(self, banner):
        """ Calculate energy for a banner at given airspeed
        function uses the aircraft object solved at __init__()
            @param banner: length of banner [m]
            @param cruise_airspeed [m/s]
            
            @return energy_total: total energy needed for given condition [Whr]
            """
        energy_tof = 0.55264     # Whr
        energy_climb = 0.3668    # Whr
        t_tof = 3.4         # s
        t_climb = 1.385  
        mass_no_banner = self.dummy_plane.mass_m3_gross - self.dummy_plane.mass_banner
        area_banner = banner**(2)/5 # m^2
        mass_banner = area_banner * 0.0359  # kg, linear trend w/ area
        mass_total = mass_no_banner + mass_banner
        cruise_airspeed = 1.2*np.sqrt(mass_total*9.81/(1/2*self.dummy_plane.rho*self.dummy_plane.S*self.dummy_plane.CLmax))


        # find cl at cruise current config
        cl_cruise = mass_total*9.81 / (1/2*self.dummy_plane.rho*cruise_airspeed**(2)*self.dummy_plane.S)

        # find drag of banner at cruise and turn
        # cd_total = 0.081644+cd_banner   # from cd_cruise-cd_banner for conceptual_sizing results
        temp_cd = self.dummy_plane.banner_cd_interp_fun(banner)
        cd_banner = max( min(self.dummy_plane.banner_data_cd), min(temp_cd, max(self.dummy_plane.banner_data_cd)))*area_banner/self.dummy_plane.S

        cd_total_cruise = self.dummy_plane.CD_m3_cruise - self.dummy_plane.CD_m3_banner_cruise  + cd_banner + self.dummy_plane.k*cl_cruise**(2)
        cd_total_turn = self.dummy_plane.CD_m3_cruise - self.dummy_plane.CD_m3_banner_cruise  + cd_banner + self.dummy_plane.k*(self.dummy_plane.n_m3*cl_cruise)**(2)

        cruise_drag = cd_total_cruise * 1/2 * self.dummy_plane.rho * cruise_airspeed**(2) * self.dummy_plane.S   # N
        turn_drag = cd_total_turn * 1/2 * self.dummy_plane.rho * cruise_airspeed**(2) * self.dummy_plane.S    # N

        # find shaft power for cruising and turning
        thrust_cruise, power_cruise, eta_cruise, rpm_cruise = self.analysis(
            airspeed=cruise_airspeed,
            drag=cruise_drag
            )
        thrust_turn, power_turn, eta_turn, rpm_turn = self.analysis(
            airspeed=cruise_airspeed,
            drag=turn_drag
            )
        if rpm_cruise == 1 or rpm_turn == 1:
            return 1
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
        return energy_total

#
    def size_banner(self, cruise_airspeed, lower_bound, upper_bound, iterations):
        """ Estimate the size of banner capable of being towed 
            description:
            banner length was previously explicitly defined, now we can try to solve
            for a banner size that is best for some given airspeed

            @param cruise_airspeed: cruising airspeed [m/s]

            @param lower_bound: lower banner length bracket [m]

            @param upper_bound: upper banner length bracket [m]

            @return banner_length: length of banner [m]
        """
        target_energy = 99
        a = lower_bound
        b = upper_bound

        fa = self.energy_from_banner(a) - target_energy
        for i in range(iterations):
            c = (a+b)/2
            f = self.energy_from_banner(c) - target_energy

            if abs(f) <= 0.1:
                print(f'object: {self.label}, converged with size_banner() after {i} iterations')
                return c
            elif f == -98:
                break
            
            if fa * f < 0:
                b = c
            else:
                a = c
                fa = f
        print(f'object: {self.label}, did not converge with size_banner()...')
        return 1

# Find what rpm is needed to get some thrust at some airspeed
# Compare all the shaft powers for each propeller
# Shaft power needed directly correlates to energy consumption
yaml_dir = 'preliminary_design/docs/props/yaml'
speed = 12.51  # m/s
drag = 23.445  # N

# Arrays for analysis
props = []
labels = []
diameters = []
pitches = []
thrusts = []
powers = []
effs = []
rpms = []
banner_lengths = []

count1 = 0
count2 = 0

for file in tqdm(os.listdir(yaml_dir), desc='Processing propeller files...'):
    path = os.path.join(yaml_dir, file)
    prop = PropellerAnalysis(path)
    Fn, P, eff, N = prop.analysis(speed, drag)
    if N ==1 :
        count1 += 1
        continue
    l = prop.size_banner(
        cruise_airspeed=speed,
        lower_bound=2,
        upper_bound=20,
        iterations=20
        )
    if l == 1:
        count2 += 1
        continue

    props.append(prop)
    labels.append(prop.label)
    diameters.append(prop.propeller.diameter)
    pitches.append(prop.propeller.pitch)
    thrusts.append(Fn)
    powers.append(P)
    effs.append(eff)
    rpms.append(N)
    banner_lengths.append(l)

print('---------------------\n',
        '  Analysis Complete  \n',
        f'Propeller analysis didnt converge {count1} times \n',
        f'banner analysis didnt converge {count2} times above that')

plt.figure()
plt.bar(labels, powers)
plt.xticks(rotation=45, ha='right')
plt.ylabel('Shaft Power (W)')
plt.title(f'Propeller Shaft Power at V={speed} m/s, T={drag} N')
plt.tight_layout()

plt.figure()
plt.bar(labels, effs)
plt.xticks(rotation=45, ha='right')
plt.ylabel('Propeller Efficiency')
plt.title(f'Propeller Efficiency at V={speed} m/s, T={drag} N')
plt.tight_layout()

plt.figure()
plt.bar(labels, banner_lengths)
plt.xticks(rotation=45, ha='right')
plt.ylabel('Banner (m)')
plt.title('Largest Banner Length for Each Propeller at 99 Whr')
plt.tight_layout()

plt.show()
