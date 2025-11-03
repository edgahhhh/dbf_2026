""" Propeller Analysis """
# pylint: disable=invalid-name
# pylint: disable=too-many-public-methods
# pylint: disable=unused-import
import os
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from propeller import Propeller
plt.rcParams['axes.grid'] = True

class PropellerAnalysis():
    """ object for finding adequate props """
    def __init__(self, propeller_info):
        """ Initialize object for each propeller """
        self.propeller = Propeller(propeller_info)
        self.label = f"{self.propeller.brand}_{self.propeller.diameter*39.37:.0f}x{self.propeller.pitch*39.37:.0f}"

    def find_rpm(self, airspeed:float, drag:float):
        """ find rpm to meet thrust requirement at airspeed 
        @param airspeed: velocity [m/s]
        @param drag: required thrust [N]
        """
        rpm_array = np.linspace(100, 50000, 3145)
        n_array = rpm_array * 2*np.pi/60  # rad/s

        threshold = 0.1  # N tolerance
        for n in n_array:
            # Step through speeds to find required frequency
            thrust = self.propeller.calculate_thrust(V=airspeed, n=n)
            if np.isclose(thrust, drag, atol=threshold):
                return n*60/2/np.pi
        print(f'object: {self.label}, did not converge, continuing...')
        return 1

    def analysis(self, airspeed:float, drag:float):
        """ perform analysis for propeller """
        rpm = self.find_rpm(airspeed, drag)
        n = rpm * 2*np.pi/60  # rad/s
        thrust = self.propeller.calculate_thrust(airspeed, n)
        power = self.propeller.calculate_power(airspeed, n)
        # eta = self.propeller.calculate_efficiency(airspeed, n)
        eta = airspeed*thrust / power
        return thrust, power, eta, rpm

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

for file in tqdm(os.listdir(yaml_dir), desc="Processing YAML files"):
    path = os.path.join(yaml_dir, file)
    prop = PropellerAnalysis(path)
    Fn, P, eff, N = prop.analysis(speed, drag)
    if N == 1:
        # error flag
        continue
    props.append(prop)
    labels.append(prop.label)
    diameters.append(prop.propeller.diameter)
    pitches.append(prop.propeller.pitch)
    thrusts.append(Fn)
    powers.append(P)
    effs.append(eff)
    rpms.append(N)

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
ax1 = plt.gcf().add_subplot(111, projection='3d')
d = np.asarray(diameters) *39.37
p = np.asarray(pitches) *39.37
sc = ax1.scatter(d, p, rpms, c=effs, marker='o')
ax1.set_xlabel('Diameter (in)')
ax1.set_ylabel('Pitch (in)')
ax1.set_zlabel('RPM')
ax1.set_title(f'Propeller RPM at V={speed} m/s, T={drag} N')
cb = plt.colorbar(sc, ax=ax1, shrink=0.6, location='left')
plt.tight_layout()

plt.figure()
ax2 = plt.gcf().add_subplot(111, projection='3d')
sc2 = ax2.scatter(d, p, powers, c=effs, marker='o')
ax2.set_xlabel('Diameter (in)')
ax2.set_ylabel('Pitch (in)')
ax2.set_zlabel('Shaft Power (W)')
ax2.set_title(f'Propeller Shaft Power at V={speed} m/s, T={drag} N')
cb2 = plt.colorbar(sc2, ax=ax2, shrink=0.6, location='left')
plt.tight_layout()
plt.show()
