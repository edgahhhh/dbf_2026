""" Banner Analysis """
# pylint: disable=invalid-name
# pylint: disable=too-many-public-methods
# pylint: disable=unused-import
import os
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from propeller_analysis import PropellerAnalysis

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
new_effs_cruise=[]
new_powers_cruise=[]
new_effs_turn=[]
new_powers_turn=[]

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
    l, V, data = prop.size_banner(
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
    rpms.append(N)#flag
    banner_lengths.append(l)
    cruise_speed.append(V)
    new_effs_cruise.append(data[2])
    new_powers_cruise.append(data[1])
    new_effs_turn.append(data[6])
    new_powers_turn.append(data[5])

# print('---------------------\n',
#         '  Analysis Complete  \n',
#         f'Propeller analysis didnt converge {count1} times \n',
#         f'banner analysis didnt converge {count2} times above that')

# Original propeller analysis stuff
f1=plt.figure()
f1.suptitle('Original Propeller Analysis')
a11 = f1.add_subplot(2,1,1)
a12 = f1.add_subplot(2,1,2)

a11.bar(labels, powers)
plt.setp(a11.get_xticklabels(), rotation=45, ha='right')
a11.set_ylabel('Shaft Power (W)')
a11.set_title(f'Propeller Shaft Power at V={speed} m/s, T={drag} N')

a12.bar(labels, effs)
plt.setp(a12.get_xticklabels(), rotation=45, ha='right')
a12.set_ylabel('Propeller Efficiency')
a12.set_title(f'Propeller Efficiency at V={speed} m/s, T={drag} N')

f1.tight_layout()

# Banner analysis stuff
f2=plt.figure()
f2.suptitle('Banner Propeller Analysis')
a21=f2.add_subplot(2,2,1)
a22=f2.add_subplot(2,2,2)
a23=f2.add_subplot(2,2,3)
a24=f2.add_subplot(2,2,4)

a21.bar(labels, banner_lengths)
plt.setp(a21.get_xticklabels(), rotation=45, ha='right')
a21.set_ylabel('Banner Length (m)')
a21.set_title('Largest Banner Length for Each Propeller at 99 Whr')

a22.bar(labels, cruise_speed)
plt.setp(a22.get_xticklabels(), rotation=45, ha='right')
a22.set_ylabel('Cruise Speed (m/s)')
a22.set_title('Cruise Speed for Each Propeller at 99 Whr w/ Largest Banner (Vcruise=1.2Vstall)')

width=0.35
x=np.arange(len(labels))
a23.bar(x-width/2, new_effs_cruise, width=width,label='cruise')
a23.bar(x+width/2, new_effs_turn, width=width,label='turn')
plt.setp(a23.get_xticklabels(), rotation=45, ha='right')
a23.set_ylabel('efficiencies')
a23.set_title('Efficiency for each propeller towing its largest banner')
a23.legend()
a23.set_xticks(x)
a23.set_xticklabels(labels=labels)

a24.bar(x-width/2, new_powers_cruise, width=width, label='cruise')
a24.bar(x+width/2, new_powers_turn, width=width, label='turn')
plt.setp(a24.get_xticklabels(), rotation=45, ha='right')
a24.set_ylabel('Shaft power (W)')
a24.set_title('Shaft powers for propellers with largest banners')
a24.legend()
a24.set_xticks(x)
a24.set_xticklabels(labels=labels)

f2.tight_layout()

# Revisiting propeller stuff
f3=plt.figure()


plt.show()
