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
    def __init__(self,
                banner_info,
                banner_length,
                banner_density=913.44,
                banner_thickness=0.0000254,
                aspect_ratio=5.0,
                ):
        """ load in banner data from yaml file
            @param banner_info: path to banner yaml file
            @param banner_length: length of banner [m]
            @param banner_density [kg/m^3]
            @param banner_density [m]
            @param aspect_ratio: aspect ratio of banner
        """
        with open(banner_info, 'r') as file:
            self._banner_info = yaml.safe_load(file)
        self._reynolds_data=self._banner_info.get('reynolds_data', None)
        self._cd_data=self._banner_info.get('c_d_data', None)
        self._cd_interp=interp1d(
            self._reynolds_data, self._cd_data, fill_value='extrapolate'
        )
        self.thickness=banner_thickness
        self.length=banner_length
        self.aspect_ratio=aspect_ratio
        self.height=banner_length/aspect_ratio

        # banner calculations
        self._banner_rho=banner_density #kg/m^3
        self._volume=self.thickness*self.length**(2)/self.aspect_ratio
        self.mass=self._banner_rho*self._volume #mass in kg
        self.area=self.length**2 / self.aspect_ratio    # towing config
        self.area_frontal=self._volume/self.height  # stowed config

        # Other stuff needed
        self._rho=density_from_alt(height_asl=400.0)    # air density kg/m^3

    def get_drag_coefficient(self, airspeed, ref_area):
        """ get saturated drag coefficient from reynolds number 
            @param reynolds_number: reynolds number of banner in flight
            @param ref_area: reference area for coefficient [m^2]
            @return c_d: drag coefficient [dimensionless]
        """
        re = reynolds_number(
            rho=self._rho,
            airspeed=airspeed,
            length_wet=self.length
        )
        cd_non_ref=max(
            min(self._cd_data), min(
                self._cd_interp(re), max(self._cd_data)))
        cd_ref = cd_non_ref * (self.area / ref_area)
        return cd_ref


    def get_stowed_drag_coefficient(self, ref_area):
        """" get drag coefficient for banner while its stowed """
        cds=1.98    # flat frontal area
        cd=cds*self.area_frontal/ref_area
        return cd

# figure out how we can split up mission coefficients into the dynamics object???
# idk tbh, maybe only leave the polar in here, as well as some equations that may bloat
# other stuff....
# Create an interpolator for polar vs aoa, cd_fus vs Re and cd_banner vs Re?
class Dynamics:
    """ Dynamics class for storing info from given drag polar 
        CL
        ^               
        |           __=__
        |         [       ]
        |      [            ]
        |    [              ]
        |  [                +
        | [
        |+
        |------------------------> alpha"""
    def __init__(self,
                ref_area,
                aspect_ratio,
                polar_info,
                fus_length,
                fus_diameter,
                banner_info,
                banner_length,
                banner_aspect_ratio=5.00,
                crud_factor=0.1
                ):
        """ Load in drag polar and estimate all other drag 
        Assume that only wing and tail drag is given in the polar
        without a crud factor """
        self._ref_area=ref_area
        self._ar=aspect_ratio
        self._k = 1/np.pi/self._ar/0.9
        self._rho=density_from_alt(height_asl=400.0)
        self._crud=crud_factor
        # Load in polar data and make interpolators
        with open(polar_info, 'r') as file:
            self._polar_info=yaml.safe_load(file)
        self._cl_data=self._polar_info.get('cL')
        self._cd_data=self._polar_info.get('cD')
        self._aoa_data=self._polar_info.get('aoa')
        self._cl_interp=interp1d(self._aoa_data,self._cl_data)
        self._cd_interp=interp1d(self._aoa_data,self._cd_data)
        self._cd_from_cl_interp=interp1d(self._cl_data,self._cd_data)

        # Load in fuselage info
        self._fus_length=fus_length
        self._fus_diameter=fus_diameter
        self._fus_fineness=fus_length/fus_diameter
        self._fus_wet_area=3/4*np.pi*fus_diameter*fus_length

        # Create banner object here? for now yes
        self._banner=Banner(
            banner_info=banner_info,
            banner_length=banner_length,
            aspect_ratio=banner_aspect_ratio
        )

    def get_drag_coefficient(self,airspeed, aoa_plane, load_factor, stowed=False):
        """ Get all drag coefficients and apply a crud factor """
        # get banner cd
        if stowed is False:
            cd_banner=self._banner.get_drag_coefficient(
                airspeed=airspeed,
                ref_area=self._ref_area
            )
        else:
            cd_banner=self._banner.get_stowed_drag_coefficient(
                ref_area=self._ref_area)
        # get fuselage cd
        cd_fus=self.fuselage_drag(airspeed=airspeed)
        # get cd and cl of wing+tail
        cl_wing_tail, cd_wing_tail=self.wing_tail_dynamics(aoa=aoa_plane)

        cdi=self.induced_drag(
            cl=cl_wing_tail,
            load_factor=load_factor,
        )
        cd=cd_banner+cd_fus+cd_wing_tail+cdi
        cd=cd/self._crud
        return cd

    def get_drag_coefficients_legacy(self,airspeed,weight,load_factor,stowed=False):
        """ Get all drag coefficients and apply crud factor 
        legacy: not using aoa (which is more of a control use case)
        """
        # get banner cd
        if stowed is False:
            cd_banner=self._banner.get_drag_coefficient(
                airspeed=airspeed,
                ref_area=self._ref_area
            )
        else:
            cd_banner=self._banner.get_stowed_drag_coefficient(
                ref_area=self._ref_area
            )
        # get fuselage cd
        cd_fus=self.fuselage_drag(airspeed=airspeed)
        # get cd of wing+tail
        cl_wing_tail=load_factor*weight*2/self._rho/airspeed**(2)/self._ref_area
        cd_wing_tail=max(
            min(self._cd_data), min(
                self._cd_from_cl_interp(cl_wing_tail), max(self._cd_data)
            )
        )
        cdi=self.induced_drag(cl=cl_wing_tail, load_factor=load_factor)
        cd=cd_banner+cd_fus+cd_wing_tail+cdi
        cd=cd/self._crud
        return cd


    def wing_tail_dynamics(self,aoa):
        """ get saturated cL and cD from wing/tail polar 
            @param aoa: angle of attack of plane (deg)
            @return cl: cL of wing/tail
            @return cd: cD of wing/tail
        """
        # cd_non_ref=max(
        #     min(self._cd_data), min(
        #         self._cd_interp(re), max(self._cd_data)))
        cl=max(
            min(self._cl_data), min(
                self._cl_interp(aoa), max(self._cl_data)
            )
        )
        cd=max(
            min(self._cd_data), min(
                self._cd_interp(aoa), max(self._cd_data)
            )
        )
        return cl, cd

    def induced_drag(self,cl,load_factor):
        """ calculate the induced drag for some loading factor """
        # TODO: check if XFLR5 does induced drag already for the polar
        cdi=(load_factor*cl)**(2)*self._k
        return cdi

    def landing_gear_drag(self):
        """ calculate the cD of landing gear """


    def fuselage_drag(self,airspeed):
        """ calculate the fuselage drag coefficient by estimating
        it as a laminar flow 
            @param airspeed: true airspeed [m/s]
        """
        # find the skin friction coefficient for the fuselage
        re_fuselage = reynolds_number(
            rho=self._rho,
            airspeed=airspeed,
            length_wet=self._fus_length)
        # TODO: account for transition and cutoff
        # estimate cd as laminar
        cf = 1.328/np.sqrt(re_fuselage)
        # cd=1/S(cf*ff*Swet*IF) assume no interference
        cd=self._fus_wet_area/self.ref_area*(
            cf*(1+self.fuselage.fineness**(-1.5))+0.11*self.fuselage.fineness**(-2))
        return cd

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
    'polar_info':'path/to/drag_polar.yaml',
    'empty_mass': 0.0,  # empty mass of the aircraft (no payload + no batts)
}

# TODO: clean up SimpleNameSpace with objects and loops and shorten lines
class BaseAircraft:
    """ Base Aircraft class """
    def __init__(self,
                aircraft_info:dict,
                mission2_info:dict,
                mission3_info:dict,
                altitude=400.0):
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
        
        # Create objects for the dynamics of the aircraft
        self._polar_info = self._aircraft_info.get(
            'polar_info', 'path/to/drag_polar.yaml')
        
        self.dynamics=Dynamics(
            polar_info=self._polar_info,
            banner_info=self._m3_banner_info,
            banner_length=self._m3_info.get('banner_length', 0.0),
            banner_aspect_ratio=5.0
        )

        self.empty_mass = self._aircraft_info.get('empty_mass', 0.0)


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
    

    def takeoff_speed(self, mass):
        return None
        


    def take_off_performance(self):
        """ estimate the performance at takeoff """
        return None



# Local functions for all objects?
def density_from_alt(height_asl):
    """ calculate density from altitude """
    theta = 1 + (-0.000022558) * height_asl
    rho = 1.225 * (theta**4.2561)
    return rho

def temp_from_alt(height_asl):
    """ calculate temperature from altitude """
    theta = 1 + (-0.000022558) * height_asl
    temp = 288.15 * theta
    return temp

def reynolds_number(rho, airspeed, length_wet, height_asl=400.0):
    """ calculate the reynolds number in air from free-stream """
    # mu = 18*10**(-6) # Pa*s
    temp=temp_from_alt(height_asl=height_asl)
    mu=1.458*10**(-6)*a(temp**(1.5)/(110.4+temp))
    re = rho * airspeed * length_wet / mu
    return re