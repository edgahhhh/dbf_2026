""" tail sizing """
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from py_planes.avl import Aircraft as AvlAircraft
from py_planes.avl import Avl, get_coef

# [X] want to say here are the tail volumes range
# [X] here is the fixed tail arm
# [X] for each value in range, create a new plane
# [X] /runs/generated_vht_VALUE.avl
# [X] run the steps in avl from the template for each plane
# [X] save the output in that one file
# [X] for each file, open, get coefficient, and add to a larger dictionary

# ht aspect ratio = 1.44608
def size_horizontal_tail(min_vht, max_vht, n_vht, ar_ht, tail_arm, wing, mass):
    """ size horizontal tail based off tail volume and tail arm """
    vht_arr=np.linspace(min_vht, max_vht, n_vht)
    tail={ # remove horizontal tail stuff
        'l_tail': tail_arm,
        'h_foil': 'naca_0012.dat',
        'h_vertical_offset': 0.1127,
        'elev_chord': 0.08,
        'v_span': 0.2254,
        'v_chord': 0.1543,
        'v_foil': 'naca_0012.dat',
        'v_vertical_offset': None,  # if the tail is offset (skipping for now)
        'rudd_chord': 0.08,
        'incidence': -2.0 # i think idk
    }
    xnps=[]
    for vht in vht_arr:
        # iterate through each tail volume in range
        # get tails dimensions
        sht=vht*wing['chord']*wing['span']/tail_arm
        bht=np.sqrt(ar_ht*sht)
        cht = bht / ar_ht
        tail['h_span']=bht
        tail['h_chord']=cht
        # create a plane from dimensions and an avl file
        plane=AvlAircraft(wing, tail, mass)
        avl=Avl(f'gen_vht_{vht:.4f}.avl',
                plane,
                'Avl',
                'runs')
        avl.create_avl_file()
        # run the steps from the template associated
        # TODO: change from tail_sizing to horizontal_tail_sizing
        # remove file that may be there already
        file_out_og=f'Avl/runs/tail_sizing/{avl.name}_st.txt'
        if os.path.exists(file_out_og):
            # remove file
            os.remove(file_out_og)
            print(f'removing file: {file_out_og}')
        
        avl.run_avl('Avl/automation/steps_template_tail_sizing.txt')
        # rename the generically named file from the template to something
        # a little more descriptive of our current use cases
        file_out_name=f'Avl/runs/tail_sizing/{avl.name}_{vht:.4f}_st.txt'
        os.rename(
            file_out_og,
            file_out_name)
        # open file, and get its TextIO object to get coefficient from
        with open(file_out_name, 'r') as stab_file:
            file_out = stab_file
            xnp=get_coef(file_out, 'Xnp')
        # outs.append(out)
        xnps.append(xnp)
        # xnps.append(out['Xnp'])
    
    f=plt.figure()
    a=f.add_subplot(111)
    a.plot(vht_arr, xnps, marker='o', linestyle='-')
    a.set_xlabel('horizontal tail volumes')
    a.set_ylabel('neutral point locations (aft leading edge), m')
    a.grid()


wing={
    'span': 1.524,
    'chord': 0.416,
    'airfoil': 'clark_y.dat',
    'aileron_span': 0.381,   # span of EACH aileron
    'aileron_chord': 0.104,
    'incidence': 2.0
}

mass={  # cog of our plane, well have multiple cases, so we can make case yamls
    'x_ref': 0.1543,
    'y_ref': 0.0,
    'z_ref': 0.0,
}

size_horizontal_tail(0.1, 2.0, 25, 1.44608, 1.5, wing, mass)
plt.show()



# Goal is to use AVL to plot a few things:
# 1. Change in neutral point with change in tail volumes
# 2. Change in 

 




# aircraft_info={
#     'wing_info':{   # wing info
#         'area':0.633984, # reference area
#         'mac':0.416, # mean aerodynamic chord
#         'span':1.524,    # wingspan
#         'aspect_ratio':3.663462,   # aspect ratio
#     },
#     'tail_info':{   # tail info (surfaces)
#         'h_area':0.09, # horizontal tail area
#         'h_mac':0.1647,  # horizontal tail mean aerodynamic chord
#         'h_span':0.5647, # horizontal tail span
#         'h_AR':3.32,   # horizontal tail aspect ratio
#         'v_area':0.0348, # vertical tail area
#         'v_mac':0.1543,  # vertical tail mean aerodynamic chord
#         'v_span':0.2254, # vertical tail span
#         'v_AR':1.46,   # vertical tail aspect ratio
#     },
#     'fuselage_info':{   # fuselage info
#         'length':1.266,    # fuselage length
#         'diameter':0.1016,  # fuselage diameter
#     },
#     'empennage_info':{  # empennage info, square boom type
#         'tail_arm':1.5,  # distance from wing AC to h-tail AC
#         'boom_length':1.088, # distance from aft fus to h-tail AC
#         'boom_width':0.0254,  # width of tail boom
#     },
#     'polar_info':'docs/polar/yaml/rio_5-15_0m_s.yaml',
#     'empty_mass': 2.286,  # empty mass of the aircraft (no payload + no batts)
# }


# mission_two_info={
#     'propulsion_info': 'docs/props/new_yaml/apce_21x13_2043.yaml',
#     'gross_mass': 2.741346,
#     'cruise_speed': 16.6,
#     'ducks': 3,
#     'pucks': 1,
# }

# mission_three_info={
#     'propulsion_info': 'docs/props/new_yaml/apce_21x13_2043.yaml',
#     'gross_mass':3.70024,
#     'post_gross_mass':3.18146, 
#     'cruise_speed': 12.0,
#     'banner_info': 'docs/isaac_banner_data_2.yaml',
#     'banner_length': 9.5,
# }

# aircraft=Aircraft(
#     aircraft_info=aircraft_info,
#     mission2_info=mission_two_info,
#     mission3_info=mission_three_info,
#     altitude=400
# )

# print(aircraft._mission_two_performance())
# print(aircraft._mission_three_performance())
# aircraft.dynamics.plot_polar()
# aircraft.drag_buildup()
# # plt.show()
