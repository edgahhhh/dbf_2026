""" """
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from py_planes.systems_plane import Aircraft

# Goal is to use AVL to plot a few things:
# 1. Change in neutral point with change in tail volumes
# 2. Change in 

 




aircraft_info={
    'wing_info':{   # wing info
        'area':0.633984, # reference area
        'mac':0.416, # mean aerodynamic chord
        'span':1.524,    # wingspan
        'aspect_ratio':3.663462,   # aspect ratio
    },
    'tail_info':{   # tail info (surfaces)
        'h_area':0.09, # horizontal tail area
        'h_mac':0.1647,  # horizontal tail mean aerodynamic chord
        'h_span':0.5647, # horizontal tail span
        'h_AR':3.32,   # horizontal tail aspect ratio
        'v_area':0.0348, # vertical tail area
        'v_mac':0.1543,  # vertical tail mean aerodynamic chord
        'v_span':0.2254, # vertical tail span
        'v_AR':1.46,   # vertical tail aspect ratio
    },
    'fuselage_info':{   # fuselage info
        'length':1.266,    # fuselage length
        'diameter':0.1016,  # fuselage diameter
    },
    'empennage_info':{  # empennage info, square boom type
        'tail_arm':1.5,  # distance from wing AC to h-tail AC
        'boom_length':1.088, # distance from aft fus to h-tail AC
        'boom_width':0.0254,  # width of tail boom
    },
    'polar_info':'docs/polar/yaml/rio_5-15_0m_s.yaml',
    'empty_mass': 2.286,  # empty mass of the aircraft (no payload + no batts)
}


mission_two_info={
    'propulsion_info': 'docs/props/new_yaml/apce_21x13_2043.yaml',
    'gross_mass': 2.741346,
    'cruise_speed': 16.6,
    'ducks': 3,
    'pucks': 1,
}

mission_three_info={
    'propulsion_info': 'docs/props/new_yaml/apce_21x13_2043.yaml',
    'gross_mass':3.70024,
    'post_gross_mass':3.18146, 
    'cruise_speed': 12.0,
    'banner_info': 'docs/isaac_banner_data_2.yaml',
    'banner_length': 9.5,
}

aircraft=Aircraft(
    aircraft_info=aircraft_info,
    mission2_info=mission_two_info,
    mission3_info=mission_three_info,
    altitude=400
)

print(aircraft._mission_two_performance())
print(aircraft._mission_three_performance())
aircraft.dynamics.plot_polar()
aircraft.drag_buildup()
# plt.show()
