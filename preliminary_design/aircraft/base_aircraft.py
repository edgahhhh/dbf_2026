""" base aircraft
    module to store some behind the scenes stuff for use in modeling
    aircraft. So any other models can inherit this.
 """
import numpy as np
import yaml
from scipy.interpolate import interp1d
from types import SimpleNamespace   # because im too lazy to make more objects
from propeller import Propeller
#pylint: disable=redefined-outer-name
#pylint: disable=line-too-long

class Banner:
    """ Banner class for use in analyzing the banner in flight """
    def __init__(self, banner_info, banner_length, aspect_ratio=5.0):
        """ load in banner data from yaml file
            @param banner_info: path to banner yaml file
            @param banner_length: length of banner [m]
            @param aspect_ratio: aspect ratio of banner
        """
        with open(banner_info, 'r') as file:
            self._banner_info = yaml.safe_load(file)
        self._reynolds_data=self._banner_info.get('reynolds_data', None)
        self._cd_data=self._banner_info.get('c_d_data', None)
        self._cd_interp=interp1d(
            self._reynolds_data, self._cd_data, fill_value='extrapolate'
        )
        self.length=banner_length
        self.aspect_ratio=aspect_ratio
        self.area=self.length**2 / self.aspect_ratio

    def get_drag_coefficient(self, reynolds_number, ref_area):
        """ get saturated drag coefficient from reynolds number 
            @param reynolds_number: reynolds number of banner in flight
            @param ref_area: reference area for coefficient [m^2]
            @return c_d: drag coefficient [dimensionless]
        """
        cd_non_ref=max(
            min(self._cd_data), min(
                self._cd_interp(reynolds_number), max(self._cd_data)))
        cd_ref = cd_non_ref * (self.area / ref_area)
        return cd_ref

aircraft_info={
    'wing_info':{   # wing info
        'area':0.0, # reference area
        'mac':0.0, # mean aerodynamic chord
        'span':0.0,    # wingspan
        'aspect_ratio':0.0,   # aspect ratio
    },
    'tail_info':{   # tail info (surfaces)
        'h_area':0.0, # horizontal tail area
        'h_mac':0.0,  # horizontal tail mean aerodynamic chord
        'h_span':0.0, # horizontal tail span
        'h_AR':0.0,   # horizontal tail aspect ratio
        'v_area':0.0, # vertical tail area
        'v_mac':0.0,  # vertical tail mean aerodynamic chord
        'v_span':0.0, # vertical tail span
        'v_AR':0.0,   # vertical tail aspect ratio
    },
    'fuselage_info':{   # fuselage info
        'length':0.0,    # fuselage length
        'diameter':0.0,  # fuselage diameter
    },
    'empennage_info':{  # empennage info, square boom type
        'tail_arm':0.0,  # distance from wing AC to h-tail AC
        'boom_length':0.0, # distance from aft fus to h-tail AC
        'boom_width':0.0,  # width of tail boom
    },
    'drag_info':'path/to/drag_polar.yaml',
    'empty_mass': 0.0,  # empty mass of the aircraft (no payload + no batts)
}

# TODO: clean up SimpleNameSpace with objects and loops and shorten lines
class BaseAircraft:
    """ Base Aircraft class """
    def __init__(self,
                aircraft_info:dict,
                mission2_info:dict,
                mission3_info:dict,
                altitude=0.0):
        """ initialize base aircraft model.
        Refer to dictionaries in source file for help in syntax.
            @param aircraft_info: dict with basic aircraft info
            @param mission2_info: dict with mission 2 info
            @param mission3_info: dict with mission 3 info
            @param altitude: altitude above sea level [m]
        """
        # Load aircraft info
            # some of these don't need to be objects so they can just
            # be simple namespaces for readability

        self._aircraft_info = aircraft_info
        self.wing=SimpleNamespace()
        self.tail=SimpleNamespace()
        self.fuselage=SimpleNamespace()
        self.empennage=SimpleNamespace()
        self.m2=SimpleNamespace()
        self.m3=SimpleNamespace()

        self.wing.area=self._aircraft_info.get('wing_info',{}).get('area', 0.0)
        self.wing.mac = self._aircraft_info.get('wing_info',{}).get('mac', 0.0)
        self.wing.span = self._aircraft_info.get('wing_info',{}).get('span', 0.0)
        self.wing.aspect_ratio = self._aircraft_info.get('wing_info',{}).get('aspect_ratio', 0.0) 

        self.tail.h.area=self._aircraft_info.get('tail_info',{}).get('h_area', 0.0)
        self.tail.h.mac=self._aircraft_info.get('tail_info',{}).get('h_mac', 0.0)
        self.tail.h.span=self._aircraft_info.get('tail_info',{}).get('h_span', 0.0)
        self.tail.h.aspect_ratio=self._aircraft_info.get('tail_info',{}).get('h_AR', 0.0)

        self.tail.v.area=self._aircraft_info.get('tail_info',{}).get('v_area', 0.0)
        self.tail.v.mac=self._aircraft_info.get('tail_info',{}).get('v_mac', 0.0)
        self.tail.v.span=self._aircraft_info.get('tail_info',{}).get('v_span', 0.0)
        self.tail.v.aspect_ratio=self._aircraft_info.get('tail_info',{}).get('v_AR', 0.0)

        self.fuselage.length=self._aircraft_info.get('fuselage_info',{}).get('length', 0.0)
        self.fuselage.diameter=self._aircraft_info.get('fuselage_info',{}).get('diameter', 0.0)
        self.fuselage.area=self._fuse_wet_area()
        self.fuselage.fineness=self.fuselage.length/self.fuselage.diameter if self.fuselage.diameter>0 else 0.0

        self.empennage.tail_arm=self._aircraft_info.get('empennage_info',{}).get('tail_arm', 0.0)
        self.empennage.boom_length=self._aircraft_info.get('empennage_info',{}).get('boom_length', 0.0)
        self.empennage.boom_width=self._aircraft_info.get('empennage_info',{}).get('boom_width', 0.0)
        self.ref_area=self.wing.area

        self._drag_info = self._aircraft_info.get(
            'drag_info', 'path/to/drag_polar.yaml')
        self.empty_mass = self._aircraft_info.get('empty_mass', 0.0)

        # Load mission 2 info
        self._m2_info = mission2_info
        self._m2_prop_info= self._m2_info.get('propulsion_info', None)
        self.m2_prop=Propeller(
            propeller_info=self._m2_prop_info, altitude=altitude)
        self.m2_gross_mass = self._m2_info.get('gross_mass', 0.0)
        self.m2_cruise_speed = self._m2_info.get('cruise_speed', 0.0)
        self.m2_ducks = self._m2_info.get('ducks', 0)
        self.m2_pucks = self._m2_info.get('pucks', 0)

        # Load mission 3 info
        self._m3_info = mission3_info
        self._m3_prop_info= self._m3_info.get('propulsion_info', None)
        self.m3_prop=Propeller(
            propeller_info=self._m3_prop_info, altitude=altitude)
        self._m3_banner_info = self._m3_info.get(
            'banner_info', 'path/to/banner.yaml')
        self.banner=Banner(
            banner_info=self._m3_banner_info,
            banner_length=self._m3_info.get('banner_length',0.0),
            aspect_ratio=5.0
        )

        self.m3_gross_mass = self._m3_info.get('gross_mass', 0.0)
        self.m3_cruise_speed = self._m3_info.get('cruise_speed', 0.0)
        self.m3_banner_length = self._m3_info.get('banner_length', 0.0)

        # Calculate density from altitude
        self.rho = self.density_from_alt(altitude)  # sea level density by default
        self._altitude=altitude

    def density_from_alt(self, height_asl):
        """ calculate density from altitude """
        theta = 1 + (-0.000022558) * height_asl
        rho = 1.225 * (theta**4.2561)
        return rho

    def _temp_from_alt(self, height_asl):
        """ calculate temperature from altitude """
        theta = 1 + (-0.000022558) * height_asl
        temp = 288.15 * theta
        return temp

    def reynolds_number(self, airspeed, length_wet):
        """ calculate the reynolds number in air from free-stream """
        # mu = 18*10**(-6) # Pa*s
        temp=self._temp_from_alt(height_asl=self._altitude)
        mu=1.458*10**(-6)*(temp**(1.5)/(110.4+temp))
        re = self.rho * airspeed * length_wet / mu
        return re

    def _fuse_wet_area(self):
        """ calculate the wetted area of the fuselage """
        s_wet=3/4*np.pi*self.fuselage.diameter*self.fuselage.length
        return s_wet

    def fuselage_drag(self, airspeed):
        """ calculate the fuselage drag coefficient """
        # find the skin friction coefficient for the fuselage
        re_fuselage = self.reynolds_number(
            airspeed=airspeed,
            length_wet=self.fuselage.length)
        # TODO: account for transition and cutoff
        # estimate cd as laminar
        cf = 1.328/np.sqrt(re_fuselage)
        # cd=1/S(cf*ff*Swet*IF)
        cd=1/self.ref_area*(
            cf*(1+1/self.fuselage.fineness)+0.11*(1/self.fuselage.fineness)**(2))
        return cd
