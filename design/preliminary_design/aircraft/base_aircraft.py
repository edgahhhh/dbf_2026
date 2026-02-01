""" base aircraft
    module to store some behind the scenes stuff for use in modeling
    aircraft. So any other models can inherit this.
 """
import numpy as np
import matplotlib.pyplot as plt
import yaml
from scipy.interpolate import interp1d
from types import SimpleNamespace   # because im too lazy to make more objects
from aircraft.propeller import Propeller
# from propeller_analysis import PropellerAnalysis
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
        self._reynolds_data=self._banner_info.get('Re', None)
        self._cd_data=self._banner_info.get('CD', None)
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
        |------------------------> alpha
    """
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
        self._cdi_data=self._polar_info.get('cDi')  # for correcting for loading factors
        self._aoa_data=self._polar_info.get('alpha')
        self._cl_interp=interp1d(self._aoa_data,self._cl_data, fill_value='extrapolate')
        self._cd_interp=interp1d(self._aoa_data,self._cd_data, fill_value='extrapolate')
        self._cd_from_cl_interp=interp1d(self._cl_data,self._cd_data, fill_value='extrapolate')
        self._cdi_from_cl_interp=interp1d(self._cl_data, self._cdi_data, fill_value='extrapolate')
        self._aoa_from_cl_interp=interp1d(self._cl_data, self._aoa_data, fill_value='extrapolate')
        self.cl_max=max(self._cl_data)

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

    def calculate_dynamics(self, airspeed, cl, cd):
        """ calculate lift and drag (eventually more) """
        q_bar=1/2*self._rho*airspeed**(2)*self._ref_area
        lift=q_bar*cl
        drag=q_bar*cd
        return lift, drag

    def plot_polar(self):
        """ plot polar for fun!! """
        f=plt.figure()
        a1=f.add_subplot(1,2,1)
        a1.plot(self._aoa_data, self._cl_data)
        a1.set_title('cL data')
        a1.set_xlabel('alpha')
        a1.set_ylabel('cL')
        
        a2=f.add_subplot(1,2,2)
        a2.plot(self._aoa_data, self._cd_data, label='cD')
        a2.plot(self._aoa_data, self._cdi_data, label='cDi', color='green')
        a2.set_title('cD data')
        a2.set_xlabel('alpha')
        a2.set_ylabel('cD')

        f.suptitle('Polar data')
        f.legend()
        a1.grid()
        a2.grid()

    def get_coefficients(self,airspeed, aoa_plane, load_factor,banner=False,stowed=False):
        """ Get coefficients and apply a crud factor to drag 
            @ return cd, cl"""
        # get banner cd
        if stowed is False and banner is True:
            cd_banner=self._banner.get_drag_coefficient(
                airspeed=airspeed,
                ref_area=self._ref_area
            )
        elif stowed is True and banner is True:
            cd_banner=self._banner.get_stowed_drag_coefficient(
                ref_area=self._ref_area)
        else:
            cd_banner=0
        # get fuselage cd
        cd_fus=self.fuselage_drag(airspeed=airspeed)
        # get cd and cl of wing+tail
        cl_wing_tail, cd_wing_tail=self.wing_tail_dynamics(aoa=aoa_plane)
        # Get landing gear coefficient
        cd_gear=self.landing_gear_drag()
        # Correct induced drag for loading factor
        cdi_old, cdi_new=self.induced_drag(cl=cl_wing_tail, load_factor=load_factor)
        # get total drag and apply the crud
        cd=cd_banner+cd_fus+cd_wing_tail+cdi_new-cdi_old+cd_gear
        cd=cd*(1+self._crud)
        return cl_wing_tail, cd

    def get_drag_coefficient_legacy(self,airspeed,weight,load_factor,banner=False,stowed=False):
        """ Get all drag coefficients and apply crud factor 
        legacy: not using aoa (which is more of a control use case)
        """
        # get banner cd
        if stowed is False and banner is True:
            cd_banner=self._banner.get_drag_coefficient(
                airspeed=airspeed,
                ref_area=self._ref_area
            )
        elif stowed is True and banner is True:
            cd_banner=self._banner.get_stowed_drag_coefficient(
                ref_area=self._ref_area)
        else:
            cd_banner=0
        # get fuselage cd
        cd_fus=self.fuselage_drag(airspeed=airspeed)
        # get cd of wing+tail
        cl_wing_tail=load_factor*weight*2/self._rho/airspeed**(2)/self._ref_area
        cd_wing_tail=max(
            min(self._cd_data), min(
                self._cd_from_cl_interp(cl_wing_tail), max(self._cd_data)
            )
        )
        # Get landing gear coefficient
        cd_gear=self.landing_gear_drag()
        # Correct induced drag for loading factor
        cdi_old, cdi_new=self.induced_drag(cl=cl_wing_tail, load_factor=load_factor)
        # get total drag and apply the crud
        cd=cd_banner+cd_fus+cd_wing_tail+cdi_new-cdi_old+cd_gear
        cd=cd*(1+self._crud)
        return cl_wing_tail, cd

    def wing_tail_dynamics(self,aoa):
        """ get saturated cL and cD from wing/tail polar 
            @param aoa: angle of attack of plane (deg)
            @return cl: cL of wing/tail
            @return cd: cD of wing/tail
        """
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

    def induced_drag(self,cl,load_factor=1.0):
        """ correct induced drag for load factor
            @param cl: lift coefficient
            @param load_factor
            @return cdi_raw: cDi from interpolator
            @return cdi_new: corrected cDi"""
        cdi_raw=max(
            min(self._cdi_data), min(
                self._cdi_from_cl_interp(cl), max(self._cdi_data)
            )
        )
        cdi_new=cdi_raw*load_factor**(2)
        return cdi_raw, cdi_new

    def landing_gear_drag(self,width=0.01,height=0.25):
        """ calculate the cD of landing gear
            @param width: width of tire in m
            @param height: height of tire in m
        """
        cd_fixed=0.5  # for fixed landing gear, Gud. type I
        cd_gear=cd_fixed*width*height/self._ref_area
        return cd_gear*2

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
        cd=self._fus_wet_area/self._ref_area*(
            cf*(1+self._fus_fineness**(-1.5))+0.11*self._fus_fineness**(-2))
        return cd

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
        # Calculate density from altitude
        self.rho = density_from_alt(height_asl=altitude)  # sea level density by default
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

        self.tail.h=SimpleNamespace()
        self.tail.h.area=self._aircraft_info.get('tail_info',{}).get('h_area', 0.0)
        self.tail.h.mac=self._aircraft_info.get('tail_info',{}).get('h_mac', 0.0)
        self.tail.h.span=self._aircraft_info.get('tail_info',{}).get('h_span', 0.0)
        self.tail.h.aspect_ratio=self._aircraft_info.get('tail_info',{}).get('h_AR', 0.0)

        self.tail.v=SimpleNamespace()
        self.tail.v.area=self._aircraft_info.get('tail_info',{}).get('v_area', 0.0)
        self.tail.v.mac=self._aircraft_info.get('tail_info',{}).get('v_mac', 0.0)
        self.tail.v.span=self._aircraft_info.get('tail_info',{}).get('v_span', 0.0)
        self.tail.v.aspect_ratio=self._aircraft_info.get('tail_info',{}).get('v_AR', 0.0)

        self.fuselage.length=self._aircraft_info.get('fuselage_info',{}).get('length', 0.0)
        self.fuselage.diameter=self._aircraft_info.get('fuselage_info',{}).get('diameter', 0.0)
        self.fuselage.area=3/4*np.pi*self.fuselage.length*self.fuselage.diameter
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
        self.m2_gross_mass = self._m2_info.get('gross_mass', 0.0)*1.1
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
        self.m3_gross_mass = self._m3_info.get('gross_mass', 0.0)*1.1
        self.m3_post_gross_mass=self._m3_info.get('post_gross_mass', 0.0)*1.1
        self.m3_cruise_speed = self._m3_info.get('cruise_speed', 0.0)
        self.m3_banner_length = self._m3_info.get('banner_length', 0.0)
            # cruise speed depends on configuration
    


        # Create objects for the dynamics of the aircraft
        self._polar_info = self._aircraft_info.get(
            'polar_info', 'path/to/drag_polar.yaml')

        self.dynamics=Dynamics(
            polar_info=self._polar_info,
            ref_area=self.ref_area,
            aspect_ratio=self.wing.aspect_ratio,
            fus_length=self.fuselage.length,
            fus_diameter=self.fuselage.diameter,
            banner_info=self._m3_banner_info,
            banner_length=self._m3_info.get('banner_length', 0.0),
            banner_aspect_ratio=5.0
        )

        self.m2_v_stall=self.stall_speed(mass=self.m2_gross_mass)
        self.m2_v_lof=self.takeoff_speed(mass=self.m2_gross_mass)
        self.m2_v_v=self.climb_speed(mass=self.m2_gross_mass)
        self.empty_mass = self._aircraft_info.get('empty_mass', 0.0)

        self.m3_v_stall=self.stall_speed(mass=self.m3_gross_mass)
        self.m3_v_lof=self.takeoff_speed(mass=self.m3_gross_mass)
        self.m3_v_v=self.climb_speed(mass=self.m3_gross_mass)

        # self.banner=Banner(
        #     banner_info=self._m3_banner_info,
        #     banner_length=self._m3_info.get('banner_length',0.0),
        #     aspect_ratio=5.0
        # )

        
        self._altitude=altitude # save altitude

    def stall_speed(self, mass):
        """ calculate the stalling speed for the aircraft """
        v_stall=np.sqrt(2*mass*9.81/self.rho/self.ref_area/self.dynamics.cl_max)
        return v_stall

    def takeoff_speed(self, mass):
        """ calculate takeoff speed for some given mass"""
        v_lof=1.556*np.sqrt(mass*9.81/self.rho/self.ref_area/self.dynamics.cl_max)
        return v_lof
    
    def climb_speed(self, mass):
        """ calculate climb speed for some mass """
        v_v=1.2*np.sqrt(2*mass*9.81/self.rho/self.ref_area/self.dynamics.cl_max)
        return v_v

    def minimum_speed(self, mass):
        """ find the minimum cruise speed for config 
        assume v_min~1.3*v_stall"""
        v_min=1.3*np.sqrt(2*mass*9.81/self.rho/self.ref_area/self.dynamics.cl_max)
        return v_min

    def ground_run(self):
        """ find ground run distance for configuration """
        # TODO: this depends on thrust, so is propeller dependant
        return None

    def optimal_speed(self):
        """ find the optimal cruise speed for config """
        # TODO: define what optimal may be going forward
        return None
    
    @staticmethod
    def _calculate_time_to_tof(v_lof, sg=17.5):
        """ Calculate time to take off """
        acceleration =v_lof**2 / 2 / sg # m/s^2
        return (v_lof/acceleration)  # seconds
    
    def take_off_performance(self,mass,sg=17.5,banner=False,stowed=False):
        """ estimate the performance at takeoff 
        return: {cl,cd,L,D,Pbatt,Ebatt,time}
        """
        # assume roll at 0 degrees
        # Find how much power is needed for lift off
        # TODO: add propeller support with propeller analysis functions
        mu=0.05
        v=self.takeoff_speed(mass=mass)/np.sqrt(2)
        cl, cd=self.dynamics.get_coefficients(
            airspeed=v,aoa_plane=0,load_factor=1,banner=banner,stowed=stowed
        )
        lift, drag=self.dynamics.calculate_dynamics(
            airspeed=v,cl=cl,cd=cd
        )
        weight=mass*9.81
        thrust=v**(2)*weight/2/9.81/sg + mu*(weight-lift) + drag
        power_shaft=thrust*v/0.6
        power_batt=power_shaft/0.9/0.9
        t_tof=self._calculate_time_to_tof(
            v_lof=self.takeoff_speed(mass=mass), sg=sg
        )
        energy_batt=power_batt*t_tof/3600
        return [cl,cd,lift,drag,power_batt,energy_batt,t_tof]
    
    @staticmethod
    def calculate_time_to_climb(climb_rate, cruise_alt=20):
        """ calculate time to finish climb sequence """
        return (cruise_alt / climb_rate)
    
    def climb_performance(self,mass,banner=False,stowed=False):
        # return: {cl,cd,L,D,Pbatt,Ebatt,time}
        weight=mass*9.81
        v_climb=self.climb_speed(mass=mass)
        cl=weight/(1/2*self.rho*v_climb**(2)*self.ref_area)
        _,cd=self.dynamics.get_drag_coefficient_legacy(
            airspeed=v_climb,weight=weight,load_factor=1,banner=banner,stowed=stowed
        )
        lift, drag=self.dynamics.calculate_dynamics(
            airspeed=v_climb,cl=cl,cd=cd
        )
        # excess_power=total_avail_power-drag*v_climb
        # v_v=excess_power/weight
        v_v=3   # m/s, shortcut and easy assumption
        power_batt=0.0
        t_climb=5
        energy_climb=1
        return [cl,cd,lift,drag,power_batt,energy_climb,t_climb]
    
    def lap_performance(self,
                        laps,
                        mass,
                        v_cruise,
                        bank_turn, 
                        banner=False):
        """ estimate the performance for a lap
            @return energy_total
            @return t_total
        """
        n=1/np.cos(bank_turn)
        weight=mass*9.81

        cl_cruise, cd_cruise = self.dynamics.get_drag_coefficient_legacy(
            airspeed=v_cruise, weight=weight, load_factor=1,banner=banner)
       
        print(cl_cruise)
       
        _, drag_cruise = self.dynamics.calculate_dynamics(
            airspeed=v_cruise, cl=cl_cruise, cd=cd_cruise)
        power_cruise_aero=drag_cruise*v_cruise
        power_cruise_batt=power_cruise_aero/0.6/0.9/0.9

        cl_turn, cd_turn = self.dynamics.get_drag_coefficient_legacy(
            airspeed=v_cruise,weight=weight,load_factor=n,banner=banner
        )
        _,drag_turn=self.dynamics.calculate_dynamics(
            airspeed=v_cruise,cl=cl_turn,cd=cd_turn
        )
        power_turn_aero=drag_turn*v_cruise
        power_turn_batt=power_turn_aero/0.6/0.9/0.9

        # Find how much time needed for total cruising
        radius_turn=v_cruise**(2)/(9.81*np.sqrt(n**(2)-1))
        t_turn_ea=radius_turn/v_cruise*np.pi
        t_circle_ea=t_turn_ea*2
        t_straight_ea = 304.8/v_cruise
        t_lap_one=t_turn_ea*2 + t_circle_ea + t_straight_ea
        t_lap_ea = t_turn_ea*2 + t_circle_ea + 2*t_straight_ea
        t_total=t_lap_one+t_lap_ea*(laps-1)

        # Calculate energies
        energy_turn_ea=power_turn_batt*t_turn_ea
        energy_circle_ea=power_turn_batt*t_circle_ea
        energy_straight_ea=power_cruise_batt*t_straight_ea
        energy_total_turning=energy_turn_ea*laps*2 + energy_circle_ea*laps
        energy_total=energy_straight_ea + 2*energy_straight_ea*(laps-1) + energy_total_turning
        energy_total=energy_total/3600
    
        return energy_total, t_total, cl_cruise

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
    mu=1.458*10**(-6)*(temp**(1.5)/(110.4+temp))
    re = rho * airspeed * length_wet / mu
    return re
