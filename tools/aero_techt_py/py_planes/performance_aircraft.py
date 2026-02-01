""" systems plane 
    preliminary design aircraft object """
# import numpy as np
# import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from py_planes.aircraft import Aircraft
# pylint: disable=redefined-outer-name

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
#     'polar_info':'path/to/drag_polar.yaml',
#     'empty_mass': 2.286,  # empty mass of the aircraft (no payload + no batts)
# }
# mission_two_info={
#     'propulsion_info': 'docs/props/new_yaml/apce_21x13_2043.yaml',
#     'gross_mass': 2.741346,
#     'cruise_speed': 0.0,
#     'ducks': 3,
#     'pucks': 1,
# }
# mission_three_info={
#     'propulsion_info': 'docs/props/new_yaml/apce_21x13_2043.yaml',
#     'gross_mass':3.70024,
#     'post_gross_mass':3.18146, 
#     'cruise_speed': 0.0,
#     'banner_info': 'docs/isaac_banner_data_2.yaml',
#     'banner_length': 1.9,
# }

class PerformanceAircraft(Aircraft):
    """ Aircraft class for detailed design use """
    def __init__(self,
                aircraft_info,
                mission2_info,
                mission3_info,
                altitude=400.0):
        """ initialize aircraft model with basic info
            @param aircraft_info: dict with basic aircraft info
            @param mission2_info: dict with mission 2 info
            @param mission3_info: dict with mission 3 info
            @param altitude: altitude above sea level [m]"""
        # Initialize with inherited __init__ for cleanliness
        super().__init__(
            aircraft_info=aircraft_info,
            mission2_info=mission2_info,
            mission3_info=mission3_info,
            altitude=altitude)

    def _mission_two_performance(self):
        """ calculate mission two performance """
        # The overall goal I THINK is to find the optimal cruise speed for mission two
        # Also want to find stuff like lap times 
        self._m2_takeoff=self.take_off_performance(
            mass=self.m2_gross_mass)
        self._m2_climb=self.climb_performance(
            mass=self.m2_gross_mass )
        self._m2_E_laps, self._m2_t_laps, _ =  self.lap_performance(
            laps=5,
            mass=self.m2_gross_mass,
            v_cruise=self.m2_cruise_speed,
            bank_turn=0.96607,
        )
        energy_total = (self._m2_E_laps+self._m2_takeoff[5]+self._m2_climb[5])/0.75
        time_total = self._m2_t_laps + self._m2_takeoff[6] + self._m2_climb[6]
        return energy_total, time_total

    def _mission_three_performance(self):
        """ calculate mission three performance """
        # find the minimum speed for mission three, this will be the cruise speed
        self.m3_cruise_speed=self.minimum_speed(
            mass=self.m3_gross_mass
        )
        # print(self.m3_cruise_speed)

        # calculate the take off performance
        # This is just a single tuple for now
        self._m3_takeoff=self.take_off_performance(
            mass=self.m3_gross_mass, banner=True, stowed=True)
        
        # calculate the climb performance
        climb_perf, climb_power, climb_energy, prop_perf = self.climb_performance(
            self.m3_gross_mass, self.m3_prop, banner=True, stowed=True)
        self._m3_climb_performance = climb_perf
        self._m3_climb_power = climb_power
        self._m3_climb_energy = climb_energy
        self._m3_climb_prop_perf = prop_perf

        self._m3_laps_perf, self._m3_turn_prop_perf = self.lap_performance(
            laps=4,
            mass=self.m3_post_gross_mass,
            v_cruise=self.m3_cruise_speed,
            bank_turn=0.5163,
            banner=True,
            propeller=self.m3_prop
        )
        self._m3_E_laps, self._m3_t_laps, cl_cruise =  self.lap_performance(
            laps=4,
            mass=self.m3_post_gross_mass,
            v_cruise=self.m3_cruise_speed,
            bank_turn=0.5163,
            banner=True,
        )

        alpha_m3 = self.dynamics._aoa_from_cl_interp(cl_cruise)
        energy_total =(self._m3_E_laps+self._m3_takeoff[5]+self._m3_climb[5])/0.8
        time_total = self._m3_t_laps + self._m3_takeoff[6] + self._m3_climb[6]

        # get the maximum drag for this mission
        # assume its from banner towing just for fun
        # or we can return
        # drag: [t/o, climb, cruise, turn, land]
        return energy_total, time_total, alpha_m3

    def drag_buildup(self):
        """ drag buildup 
        graphs for each config and what are the cd's, bar graphs? """
        # wing/tail, fuselage, landing gear, payloads
        # m2 stuff
        m2_weight = self.m2_gross_mass*9.81
        cl_m2, cd_m2_total = self.dynamics.get_drag_coefficient_legacy(
            airspeed=self.m2_cruise_speed,
            weight=m2_weight,
            load_factor=1
        )
        cd_m2_fuselage=self.dynamics.fuselage_drag(airspeed=self.m2_cruise_speed)

        cd_m2_wing_tail = self.dynamics._cd_from_cl_interp(cl_m2)

        m3_weight = self.m3_gross_mass*9.81
        m3_post_weight = self.m3_post_gross_mass*9.81

        cl_m3_stowed, cd_m3_total_stowed=self.dynamics.get_drag_coefficient_legacy(
            airspeed=self.m3_v_lof,
            weight=m3_weight,
            load_factor=1.0,
            banner=True,
            stowed=True
        )

        _,cd_m3_total_tow = self.dynamics.get_drag_coefficient_legacy(
            airspeed=self.m3_cruise_speed,
            weight=m3_post_weight,
            load_factor=1.0,
            banner=True,
            stowed=False
        )
        cd_m3_wing_tail = self.dynamics._cd_from_cl_interp(cl_m3_stowed)

        cd_m3_fuselage_stowed = self.dynamics.fuselage_drag(
            airspeed=self.m3_v_lof
        )

        cd_m3_fuselage_tow = self.dynamics.fuselage_drag(
            airspeed=self.m3_cruise_speed
        )

        cd_banner_stowed = self.dynamics._banner.get_stowed_drag_coefficient(
            ref_area=self.ref_area,
        )
        cd_banner_tow = self.dynamics._banner.get_drag_coefficient(
            airspeed=self.m3_cruise_speed,
            ref_area=self.ref_area
        )

        cd_gear = self.dynamics.landing_gear_drag()


        f=plt.figure()
        a1=f.add_subplot(1,3,1)
        m2_cats=['Total', 'WingTail', 'Fuselage', 'Gear']
        m2_drags=[cd_m2_total, cd_m2_wing_tail, cd_m2_fuselage, cd_gear]

        a1.set_title('Mission two')
        a1.bar(m2_cats, m2_drags)
        a1.set_ylabel('Drag Coefficient (CD)')
        a1.grid()

        a2=f.add_subplot(1,3,2)
        m3_cats=['Total', 'Wingtail', 'Fuselage', 'Gear', 'Banner']
        m3_drags_stowed=[cd_m3_total_stowed, cd_m3_wing_tail, cd_m3_fuselage_stowed, cd_gear, cd_banner_stowed]
        a2.bar(m3_cats, m3_drags_stowed)
        a2.set_title('Mission three stowed')
        a2.grid()

        a3=f.add_subplot(1,3,3)
        m3_drags_tow=[cd_m3_total_tow, cd_m3_wing_tail, cd_m3_fuselage_tow, cd_gear, cd_banner_tow]
        a3.bar(m3_cats, m3_drags_tow)
        a3.set_title('Mission three towing')
        a3.grid()

        f.suptitle('Drag Buildup')

        # create another figure that shows drag numbers for each mission
        self._

    def calculate_performance(self,):
        """ this should calculate performance of the aircraft in all missions"""
        return None


def plt_3d(x,y,z,title,xlabel,ylabel,zlabel,
           c=None, clabel=None):
    """ plot 3d points """
    plt.figure()
    ax=plt.gcf().add_subplot(111,projection='3d')
    ax.set_title(title)
    if c is None:
        sc=ax.scatter(x,y,z)
    else:
        sc=ax.scatter(x,y,z,c=c,cmap='magma')
        cb=plt.gcf().colorbar(sc, ax=ax, pad=0.1)
        if clabel is None:
            # set to zlabel if not given
            cb.set_label(zlabel)
        else:
            cb.set_label(clabel)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    return sc
