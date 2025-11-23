""" Banner Analysis """
# pylint: disable=invalid-name
# pylint: disable=too-many-public-methods
# pylint: disable=unused-import
import os
from tqdm import tqdm
from propeller_analysis import PropellerAnalysis
import matplotlib.pyplot as plt

# Find what rpm is needed to get some thrust at some airspeed
# Compare all the shaft powers for each propeller
# Shaft power needed directly correlates to energy consumption
yaml_dir = 'preliminary_design/docs/props/new_yaml'
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
cruise_speed = []

count1 = 0
count2 = 0
sanity_count = 0
for file in tqdm(os.listdir(yaml_dir), desc='Processing propeller files...'):
    path = os.path.join(yaml_dir, file)
    prop = PropellerAnalysis(path) 
    if sanity_count == 0:
        print('======================= \n',
              '  sanity check here \n\n',)
        prop.dummy_plane.size_aircraft_all_missions(show_summary=False)
        print()
        sanity_count += 1

    # prop.propeller.plot_original_data()
    Fn, P, eff, N = prop.analysis(speed, drag)
    if N ==1 :
        count1 += 1
        continue
    l, V = prop.size_banner(
        lower_bound=2,
        upper_bound=20,
        iterations=20
        )
    if l == 1 :
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
    cruise_speed.append(V)

# print('---------------------\n',
#         '  Analysis Complete  \n',
#         f'Propeller analysis didnt converge {count1} times \n',
#         f'banner analysis didnt converge {count2} times above that')

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

# print(len(props))
# print(len(l))

plt.figure()
plt.bar(labels, banner_lengths)
plt.xticks(rotation=45, ha='right')
plt.ylabel('Banner Length (m)')
plt.title('Largest Banner Length for Each Propeller at 99 Whr')
plt.tight_layout()

plt.figure()
plt.bar(labels, cruise_speed)
plt.xticks(rotation=45, ha='right')
plt.ylabel('Cruise Speed (m/s)')
plt.title('Cruise Speed for Each Propeller at 99 Whr w/ Largest Banner (Vcruise=1.2Vstall)')
plt.tight_layout()

plt.show()
