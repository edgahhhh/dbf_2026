import numpy as np
import matplotlib.pyplot as plt
import yaml
from scipy.interpolate import interp1d

class Aircraft():
    def __init__(self, aircraft_design_parameters, mission_parameters, course_parameters, initial_guess, banner_data, banner_stowing_config=1):
        """ Initialize aircraft model, initialize model each time you size aircraft """
        # Aircraft Parameters
        self.AR = aircraft_design_parameters['aspect ratio']
        self.b = aircraft_design_parameters['wing span']
        self.S = self.b**2/self.AR
        self.CLmax = aircraft_design_parameters['CLmax']
        self.CDmin_no_fus = aircraft_design_parameters['CDmin_no_fuselage']
        self.k = 1/(np.pi*self.AR*0.85)
        self.fineness_ratio = aircraft_design_parameters['fineness ratio']

        # Mission parameters, instead of laps, cruising speed is used for simplification
        self.cargo = mission_parameters['cargo']
        self.ducks = mission_parameters['ducks']
        self.V_m2_cruise = mission_parameters['mission two cruise speed']
        self.l_banner = mission_parameters['banner length']
        self.AR_banner = mission_parameters['banner aspect ratio']
        self.V_m3_cruise = mission_parameters['mission three cruise speed']
        self.S_banner = self.l_banner**(2)/self.AR_banner
        self.rac = 0.05*self.b*3.281+0.75
        self.banner_stowing_config = banner_stowing_config  # 1=longitudinal cylinder, 2=and so on....

        # Atmospheric conditions
        self.alt = aircraft_design_parameters['conditions']['altitude asl']
        self.cruise_alt = aircraft_design_parameters['conditions']['cruise altitude agl']
        self.rho = self.density_from_alt(height_asl=(self.alt+self.cruise_alt)) # use density at cruising altitude since most of mission is there

        # Course parameters, generally won't change
        self.ground_run = course_parameters['ground run']
        self.phi_m2 = aircraft_design_parameters['m2 bank angle']
        self.phi_m3 = aircraft_design_parameters['m3 bank angle']
        self.time_limit = course_parameters['time limit']
        self.length_straight = course_parameters['length of straights']

        self.n_m2 = 1/np.cos(self.phi_m2)  # load factor during m2 turns
        self.n_m3 = 1/np.cos(self.phi_m3)  # load factor during m3 turns

        # Initial mass build up
        self.mass_fuselage = 0.0
        self.mass_motor_initial_guess = initial_guess['motor mass']  # Initial motor mass guess, kg
        self.mass_landing_gear = 0.500    # From UCLA landing gear, kg, each aircraft has this
        self.mass_ducks =  self.ducks * 0.7 / 35.274    # mass of ducks, kg,
        self.mass_cargo = self.cargo * 6 / 35.274  # mass of cargo, kg
        self.mass_wing = self.wing_mass_from_area()
        self.mass_tail = self.tail_mass_from_wing_area()
        self.mass_m2_battery_initial_guess = initial_guess['m2 battery mass']  # Initial battery mass guess, kg
        self.mass_m3_battery_inital_guess = initial_guess['m3 battery mass']
        self.mass_banner = self.banner_mass_from_area()
        self.banner_stowing_electronics = 0.02  # I think 20 grams is reasonable, thats like three servos
        self.mass_propulsion_electronics = 0.5  # estimation from 2020 planes, including prop, spinner, reciever, RX battery, etc...
        
        # Banner drag data
        with open(banner_data, 'r') as file:
            self.banner_data= yaml.safe_load(file)
        self.banner_data_l = self.banner_data['length']
        self.banner_data_cd = self.banner_data['CD']
        self.banner_cd_interp_fun = interp1d(self.banner_data_l, self.banner_data_cd, fill_value='extrapolate')

        # Systems parameters, more or less placeholders for the time being, estimating power flow in take off and cruise
        self.eff_propeller_tof = 0.4
        self.eff_motor_tof = 0.9
        self.eff_esc_tof = 0.9
        self.eff_battery_tof = 0.95
        self.eff_system_tof = self.eff_propeller_tof * self.eff_motor_tof * self.eff_esc_tof * self.eff_battery_tof
        self.eff_propeller_cruise = 0.4
        self.eff_motor_cruise = 0.9
        self.eff_esc_cruise = 0.9
        self.eff_battery_cruise = 0.95
        self.eff_system_cruise = self.eff_propeller_cruise * self.eff_motor_cruise * self.eff_esc_cruise * self.eff_battery_cruise

        # Sizing values
        self.iteration_m2 = []
        self.iteration_all = []
        self.mass_empty_list = []
        self.mass_m2_gross_list = []
        self.mass_m3_gross_list = []
        self.P_m2_motor = []
        self.P_all_motor = []    
        self.E_battery_m2_list = []
        self.E_battery_m3_list = []
        self.mass_m2_battery = []
        self.mass_m3_battery = []
        self.mass_motor = []

        # Initialize variables that are mainly going to be used once the aircraft is sized
            # What matters here is the actual final parameters of aircraft, can have a function that populates this after sizing
        # performance
        self.wing_loading_m2 = 0.0
        self.power_loading_m2 = 0.0
        self.wing_loading_m3 = 0.0
        self.power_loading_m3 = 0.0 

        # Weight buildup
        self.mass_m1_gross = 0.0
        self.mass_m2_gross = 0.0
        self.mass_m3_gross = 0.0

        # Wing sizing results, delete these!!! elegantly
        self.wing_area = float(self.S)
        self.wing_aspect_ratio = float(self.AR)
        self.wing_wingspan = float(self.b)
        self.wing_chord = float(self.b/self.AR)
        self.wing_mass = float(self.mass_wing)

        # Tail sizing results
        self.htail_area = 0.0
        self.htail_vol_m1 = 0.0
        self.htail_vol_m2 = 0.0
        self.htail_vol_m3 = 0.0
        self.vtail_area = 0.0
        self.vtail_vol_m1 = 0.0
        self.vtail_vol_m2 = 0.0
        self.vtail_vol_m3 = 0.0
        self.tail_mass = 0.0

        # Fuselage sizing results
        self.d = 0.0
        self.l = 0.0
        self.S_fuselage = 0.0

        # Propulsion sizing results
        self.P_motor = 0.0  # constant, can't change the motor
        self.motor_mass = 0.0   # same as before
        self.E_battery_m2 = 0.0
        self.E_battery_m3 = 0.0
        self.battery_mass_m2 = 0.0
        self.battery_mass_m3 = 0.0

        # Mission Two Performance
        self.V_m2_st = 0.0
        self.V_m2_tof = 0.0
        self.V_m2_v = 0.0
        self.CL_m2_tof = 0.0
        self.CD_m2_fus_lof = 0.0
        self.CD_m2_fus_cruise = 0.0
        self.CL_m2_cruise = 0.0
        self.CD_m2_min_tof = 0.0
        self.CD_m2_min_cruise = 0.0
        self.CD_m2_tof = 0.0
        self.CD_m2_cruise = 0.0
        self.P_m2_tof_avail = 0.0 # T*V at take off roll
        self.P_m2_total_climb_batt = 0.0   # Total power at battery during climb 
        self.P_m2_avail_climb = 0.0  # Power at battery going to thrust during climb
        self.P_m2_exc_climb = 0.0     # Excess power at battery during climb, directly going to climb rate
        self.time_m2_tof = 0.0
        self.time_m2_climb = 0.0
        self.time_m2_cruise = 0.0
        self.time_m2_total = 0.0
        self.E_m2_tof = 0.0
        self.E_m2_climb = 0.0
        self.E_m2_cruise = 0.0
        self.L_m2_turn = 0.0    # use load factor instead
        self.laps_m2 = 0.0
        self.mission_two_score = 0.0

        # Mission Three Performance
        self.V_m3_st = 0.0
        self.V_m3_tof = 0.0
        self.V_m3_v = 0.0
        self.CL_m3_tof = 0.0
        self.CD_m3_min_tof = 0.0
        self.CD_m3_fus_lof = 0.0
        self.CD_m3_banner_tof = 0.0
        self.CD_m3_tof = 0.0
        self.CL_m3_cruise = 0.0
        self.CD_m3_min_cruise = 0.0
        self.CD_m3_fus_cruise = 0.0
        self.CD_m3_banner_cruise = 0.0
        self.CD_m3_cruise = 0.0
        self.P_m3_tof_avail = 0.0 # T*V at take off roll
        self.P_m3_total_climb_batt = 0.0   # Total power at battery during climb 
        self.P_m3_avail_climb = 0.0  # Power at battery going to thrust during climb
        self.P_m3_exc_climb = 0.0     # Excess power at battery during climb, directly going to climb rate
        self.time_m3_tof = 0.0
        self.time_m3_climb = 0.0
        self.time_m3_cruise = 0.0
        self.time_m3_total = 0.0
        self.E_m3_tof = 0.0
        self.E_m3_climb = 0.0
        self.E_m3_cruise = 0.0
        self.L_m3_turn = 0.0    # use load factor instead
        self.laps_m3 = 0.0
        self.mission_three_score = 0.0

        self.ground_mission_score = 0.0
        self.total_mission_score = 0.0


    def density_from_alt(self, height_asl):
        """ calculate density from altitude """
        theta = 1 + (-0.000022558) * height_asl
        rho = 1.225 * (theta**4.2561)
        return rho

    def reynolds_number(self, freestream_speed, length_wet):
        """ calculate the reynolds number in air from freestream """
        mu = 18*10**(-6) # Pa*s, need to get a better approximation for this as well
        re = self.rho * freestream_speed * length_wet / mu
        return re

    def fuselage_mass_from_empty_no_fus(self, mass_empty_no_fuselage):
        """ calculate mass of fuselage in kg
        @param mass_empty_no_fuselage: empty mass - fuselage mass """
        me = mass_empty_no_fuselage/(1-0.22)
        return(me*0.22)

    def calculate_fuselage_drag(self, reynolds_number, S_wet):
        """ calculate profile drag of fuselage using Horner estimation for streamline shapes """
        cf_fuselage = 1.328 / np.sqrt(reynolds_number)
        CD_fuselage_wet = cf_fuselage*(1+(1/self.fineness_ratio)**(3/2)) + 0.11*(1/self.fineness_ratio)**2
        CD_fuselage= CD_fuselage_wet * S_wet / self.S
        return CD_fuselage

    def battery_mass_from_energy(self, battery_energy):
        """ estimate mass of battery from total energy capacitance """
        energy_density = 0.00613805 # kg/Whr
        return (abs(battery_energy*energy_density))

    def battery_volume_from_mass(self, battery_mass):
        """ Calculate battery volume from mass """
        battery_mass_density = 1988.915 # kg/m^3
        return(battery_mass / battery_mass_density)

    def motor_mass_from_power(self, max_continuous_shaft_power):
        """ estimate motor mass from power required to drive shaft/propeller at take off
        assuming linear increase in power and motor mass, could use a better model (FRESHMAN!!!!!) 
        """
        motor_power_loading = 13.9194188   # W/g
        return (abs(max_continuous_shaft_power / motor_power_loading / 1000))

    def motor_volume_from_mass(self, motor_mass):
        """ Calculate battery mass from power"""
        motor_mass_density = 3310.70094  # kg/m^3
        return (motor_mass/motor_mass_density)
    
    def wing_mass_from_area(self):
        """ Calculate wing mass in kg from area in m^2 """
        return(self.S * 1.3247)
    
    def tail_mass_from_wing_area(self):
        """ calculate tail mass in kg from wing area in m^2 """
        return(self.S * 0.1996)

    def banner_mass_from_area(self):
        """ calculate banner mass from area in kg form banner area in m^2 """
        return(self.S_banner * 0.0359)

    def calculate_time_to_tof(self, Vlof):
        """ Calculate time to take off """
        acceleration = Vlof**2 / 2 / self.ground_run # m/s^2
        return (Vlof/acceleration)  # seconds

    def calculate_climb_rate(self, mass, CDmin, total_power_available, stall_speed, gravity=9.81):
        """ calculate rate of climb from best climb airspeed derived empirically from Gudmundsson's """  
        # wing_loading_lb_ft_2 = mass*gravity/self.S * 0.0208854  # wing loading in lbf/ft^2
        # V_opt = 0.5144*(43.591 + 2.2452*wing_loading_lb_ft_2)   # optimal airspeed during climb in m/s
        # IF OMITTING EMPIRICAL SOLUTION, climb airspeed ~ 1.2*Vstall
        V_climb = 1.2*stall_speed
        cd_climb = CDmin + self.k*(2*mass*gravity/self.rho/V_climb**(2)/self.S)**2    # assume climb angle is less than 15 degrees
        thrust = cd_climb*1/2*self.rho*V_climb**2*self.S
        available_power = thrust * V_climb
        excess_power = total_power_available - available_power # Excess power available, not at battery
        return (excess_power/mass/gravity)  # climb rate, m/s

    def calculate_time_to_complete_climb(self, climb_rate):
        """ calculate time to finish climb sequence """
        return (self.cruise_alt / climb_rate)

    def calculate_take_off_power(self, mass, CDtof, CLtof, mu=0.05, gravity = 9.81):
        """ estimate power at battery to achieve take off, assuming a constant acceleration """
        wing_loading = mass*gravity/self.S 
        thrust = (1.21/gravity/self.rho/self.CLmax/self.ground_run*wing_loading + 0.605/self.CLmax*(CDtof-mu*CLtof) + mu)*mass*gravity
        available_power = thrust * np.sqrt(2*wing_loading/self.rho/self.CLmax)/np.sqrt(2)
        return available_power  # Watts, if everything is in SI
    
    # def calculate_power_needed_for_level_turn(self, mass, airspeed, CDmin, load_factor, gravity=9.81):
    #     """ calculate available power needed to complete a level turn at some angle and airspeed """
    #     wing_loading = mass*gravity/self.S
    #     load_factor = 1/load_factor
    #     p_dynamic = 1/2*self.rho*airspeed**2
    #     thrust = p_dynamic*self.S*(CDmin + self.k*(load_factor*wing_loading/p_dynamic)**(2))
    #     return(thrust*airspeed)

    def calculate_power_needed_for_cruise(self, mass, airspeed, CD_cruise, CDmin, gravity=9.81):
        """ calculate available power needed to cruise at steady state (ignore trim drag) """
        wing_loading = mass*gravity/self.S
        p_dynamic = 1/2*self.rho*airspeed**2
        thrust = mass*gravity*(p_dynamic*CDmin/wing_loading + self.k/p_dynamic*wing_loading)

        return(thrust*airspeed)

    def calculate_cruise_performance(self, level_power_batt, turn_power_batt, mass, cruise_speed, time_tof, time_climb, phi, gravity=9.81):
        """ 
        calculate laps, time, and energy needed for cruise 
            @return laps [num]
            @return time_cruise_total [t]
            @return E_cruise [Whr]
        """
        n = 1/np.cos(phi)
        # radius_turn = mass * cruise_speed**2 / (1/np.cos(phi)*mass*gravity)/np.sin(phi)
        radius_turn = cruise_speed**2 / (gravity *np.sqrt(n**(2)-1))  # m
        time_turn_ea = radius_turn/cruise_speed*(np.pi)
        time_circle_ea = time_turn_ea*2
        time_straights_ea = self.length_straight/cruise_speed   # same for each unique cruising speed
        time_lap_1_temp = time_tof + time_climb + 2*time_turn_ea + time_straights_ea + time_circle_ea
        time_after_lap_1 = self.time_limit - time_climb - time_tof - time_lap_1_temp  # How much time left after first lap in seconds
        time_ea_lap = 2*time_straights_ea + 2*time_turn_ea + time_circle_ea
        laps = np.floor(1+time_after_lap_1/time_ea_lap)
        n_m2s = laps*4
        n_straights = laps*2 - 1
        time_cruise_total = time_after_lap_1 + time_ea_lap*(laps-1)
        E_cruise = (turn_power_batt*n_m2s*time_turn_ea + level_power_batt*n_straights*time_straights_ea)/3600  # fixed change here
        return laps, time_cruise_total, E_cruise

    def calculate_ground_mission_score(self):
        """" calculate ground mission score """
        t_per_duck = 2
        t_per_cargo = 2.5
        t_banner = 5
        t_best = t_per_duck*3 + t_per_cargo + t_banner
        t_gm = t_per_duck*self.ducks + t_per_cargo*self.cargo + t_banner
        self.ground_mission_score = max(0, min(t_best/t_gm, 1))
        return(self.ground_mission_score)

    def calculate_mission_two_score(self):
        """ calculate total mission two score """
        best_net = 3400
        income = self.ducks*(6 + 2*self.laps_m2) + self.cargo*(10 + 8*self.laps_m2)
        cost = self.E_battery_m2/100*self.laps_m2*(10+self.ducks*0.5+self.cargo*2)
        self.m2_net_income = income - cost
        self.mission_two_score = 1 + max(0, min((self.m2_net_income)/best_net, 1))  # was 2300
        return(self.mission_two_score)

    def calculate_mission_three_score(self):
        """ calculate mission three score """
        best_cost = 1720
        length_in = self.l_banner*39.37
        b_ft = self.b/3.281
        m3_cost = length_in*self.laps_m3/(0.05*b_ft+0.75)
        self.m3_cost = m3_cost
        self.mission_three_score =2 + max(0, min(m3_cost/best_cost, 1))
        return(self.mission_three_score)

    def calculate_total_mission_score(self):
        """ wrapper to calculate total mission score """
        self.total_mission_score = 1 + self.calculate_mission_two_score() + self.calculate_mission_three_score() + self.calculate_ground_mission_score()

    def calculate_banner_stowed_drag(self):
        """ calculate banner stowed drag based off stowing method """
        if self.banner_stowing_config == 1:
            """ cylinder with axial air flow (i.e longitudinally placed)"""
            self.CD_m3_banner_tof = 0.2 # Horner 3-12 fig 21 
    
    def size_m2_fuselage(self, battery_volume, motor_volume):
        """ wrapper to size fuselage from mission two parameters """
        volume_m2_electronics = (battery_volume + motor_volume)
        volume_m2_payload = self.ducks*(2.5/39.37)**3 + self.cargo * np.pi/4*(3/39.37)**(2) # m^3
        volume_m2_gross = volume_m2_electronics + volume_m2_payload
        self.d = (volume_m2_gross/self.fineness_ratio)**(1/3) # m
        self.l = self.fineness_ratio * self.d # m
        if self.l > 3:  # if length is above 10 ft
            self.l = 3
            self.d = np.sqrt(volume_m2_gross / self.l)
        elif self.d < 0.08: # if base is less than 3 in
            self.d = 0.08
            self.l = volume_m2_gross / self.d**2
        self.S_fuselage = 4*self.d*self.l + 2*self.d**2 #m^2
        # take-off drag
        Re_fus_tof = self.reynolds_number(self.V_m2_tof, self.l)
        self.CD_m2_fus_lof = self.calculate_fuselage_drag(Re_fus_tof, self.S_fuselage)
        # cruise drag
        Re_fus_cruise = self.reynolds_number(self.V_m2_cruise, self.l)  #m^2/s^2/Pa
        self.CD_m2_fus_cruise = self.calculate_fuselage_drag(Re_fus_cruise, self.S_fuselage)

    def size_m3_fuselage(self):
        """ calculate drag coefficients of fuselage relevant to mission three parameters """
        # take-off drag
        Re_fus_tof = self.reynolds_number(self.V_m3_tof, self.l)
        self.CD_m3_fus_lof = self.calculate_fuselage_drag(Re_fus_tof, self.S_fuselage)
        # cruise drag
        Re_fus_cruise = self.reynolds_number(self.V_m3_cruise, self.l)  #m^2/s^2/Pa
        self.CD_m3_fus_cruise = self.calculate_fuselage_drag(Re_fus_cruise, self.S_fuselage)

    def size_aircraft_all_missions(self, show_summary=False, debug_text=False):
        """ Size aircraft to meet parameters for mission two and mission three. 
        philosophy: 
            motor: sized to meet greatest power demanded (between takeoff and turning at 45 degree for both missions)
            battery: sized for both mission two and mission three
            """
        sizing_error = 1
        sizing_error_m2_batt = []
        sizing_error_m3_batt = []
        sizing_error_motor = []
        for i in range(10):
            self.iteration_all.append(i)
            if i == 0:
                """ initial guess """
                self.mass_m2_battery.append(self.mass_m2_battery_initial_guess)
                self.mass_m3_battery.append(self.mass_m3_battery_inital_guess)
                self.mass_motor.append(self.mass_motor_initial_guess)
            elif i > 0:
                """ refining guess """
                self.mass_m2_battery.append(self.battery_mass_from_energy(self.E_battery_m2_list[i-1]))
                self.mass_m3_battery.append(self.battery_mass_from_energy(self.E_battery_m3_list[i-1]))
                self.mass_motor.append(self.motor_mass_from_power(self.P_all_motor[i-1]))
            # Weight estimation
            mass_empty_no_fuselage = self.mass_motor[i] + self.mass_wing + self.mass_tail + self.mass_landing_gear + self.mass_propulsion_electronics
            self.mass_fuselage = self.fuselage_mass_from_empty_no_fus(mass_empty_no_fuselage)
            self.mass_empty_list.append(mass_empty_no_fuselage+self.mass_fuselage)
            self.mass_m2_gross_list.append(self.mass_empty_list[i]+self.mass_m2_battery[i]+self.mass_ducks+self.mass_cargo)
            self.mass_m3_gross_list.append(self.mass_empty_list[i]+self.mass_m3_battery[i]+self.mass_banner+self.banner_stowing_electronics)
            # Estimate stalling speed and take off speeds from mass
            self.V_m2_st = np.sqrt(2*self.mass_m2_gross_list[i]*9.81/self.rho/self.S/self.CLmax)
            self.V_m2_tof = 1.1*self.V_m2_st
            self.V_m3_st = np.sqrt(2*self.mass_m3_gross_list[i]*9.81/self.rho/self.S/self.CLmax)
            self.V_m3_tof = 1.1*self.V_m3_st

            # -- Mission Two -- #
            # Goal, size fuselage from volume of internal payload. Size propulsion system
            if self.mass_m2_battery[i] > self.mass_m3_battery[i]:
                """ Size fuselage for m2 battery """
                volume_battery = self.battery_volume_from_mass(self.mass_m2_battery[i])*1.1  # x1.1 to account for structure
            else:
                """ Size fuselage for m3 battery """
                volume_battery = self.battery_volume_from_mass(self.mass_m3_battery[i])*1.1  # x1.1 to account for structure
            volume_motor = self.motor_volume_from_mass(self.mass_motor[i])
            self.size_m2_fuselage(volume_battery, volume_motor)    # Populating fuselage dimensions and drag coefficients

            # take off performance
            self.CL_m2_tof = self.mass_m2_gross_list[i]*9.81*2/self.rho/self.V_m2_tof**2/self.S
            self.CD_m2_min_tof = self.CDmin_no_fus + self.CD_m2_fus_lof
            self.CD_m2_tof = self.CD_m2_min_tof + self.k*self.CL_m2_tof**2
            self.drag_m2_tof = 1/2*self.CD_m2_tof*self.rho*self.V_m2_tof**(2)*self.S

            # cruise performance
            self.CL_m2_cruise = self.mass_m2_gross_list[i]*9.81*2/self.rho/self.V_m2_cruise**2/self.S
            self.CD_m2_min_cruise = self.CDmin_no_fus + self.CD_m2_fus_cruise
            self.CD_m2_cruise = self.CD_m2_min_cruise + self.k*self.CL_m2_cruise**2
            self.drag_m2_cruise = 1/2*self.CD_m2_cruise*self.rho*self.V_m2_cruise**(2)*self.S
            p_dyn_m2_turn = 1/2*self.rho*self.V_m2_cruise**(2)
            self.drag_m2_turn = p_dyn_m2_turn*self.S*(self.CD_m2_min_cruise + self.k*(self.n_m2*self.CL_m2_cruise)**(2))

            # Constraint analysis between take off and constant velocity level turn
            self.P_m2_tof_avail = self.calculate_take_off_power(self.mass_m2_gross_list[i], self.CD_m2_tof, self.CL_m2_tof, mu=0.05)
            P_m2_tof_batt = self.P_m2_tof_avail/self.eff_system_tof
            P_m2_turn_avail = self.drag_m2_turn * self.V_m2_cruise
            P_m2_turn_batt = P_m2_turn_avail/self.eff_system_cruise
            if P_m2_tof_batt > P_m2_turn_batt:
                """ size climb to take off power """
                self.P_m2_total_climb_batt = P_m2_tof_batt
                P_m2_motor = self.P_m2_tof_avail / self.eff_propeller_tof
                P_sizing_batt_m2 = P_m2_tof_batt
            else:
                """ size climb to turning power """
                self.P_m2_total_climb_batt = P_m2_turn_batt
                P_m2_motor = P_m2_turn_avail / self.eff_propeller_cruise
                P_sizing_batt_m2 = P_m2_turn_batt

            # Energy needed for take off
            self.time_m2_tof = self.calculate_time_to_tof(self.V_m2_tof)
            self.E_m2_tof = self.P_m2_tof_avail/self.eff_system_tof*self.time_m2_tof/3600

            # Energy needed for climb, assuming climb angle is relatively small so V_gamma = V_horizontal
            P_m2_climb_avail = self.P_m2_total_climb_batt * self.eff_system_tof
            self.V_m2_v = self.calculate_climb_rate(self.mass_m2_gross_list[i], self.CD_m2_fus_lof, P_m2_climb_avail, self.V_m2_st)
            self.time_m2_climb = self.calculate_time_to_complete_climb(self.V_m2_v)
            self.E_m2_climb = self.P_m2_total_climb_batt*self.time_m2_climb/3600

            # Energy needed for cruise segment
            P_m2_level_batt = self.drag_m2_cruise*self.V_m2_cruise / self.eff_system_cruise
            self.laps_m2, self.time_m2_cruise, self.E_m2_cruise = self.calculate_cruise_performance(P_m2_level_batt, P_m2_turn_batt, self.mass_m2_gross_list[i], self.V_m2_cruise, self.time_m2_tof, self.time_m2_climb, self.phi_m2)

            # -- Mission Three -- #
            # Goal: Size propulsion system for mission three by accounting for banner stowed and free drag
            self.size_m3_fuselage()

            # take off performance
            self.CL_m3_tof = self.mass_m3_gross_list[i]*9.81*2/self.rho/self.V_m3_tof**2/self.S
            self.CD_m3_min_tof = self.CDmin_no_fus + self.CD_m3_fus_lof
            self.calculate_banner_stowed_drag() # populate CD_m3_banner_tof from chosen configuration
            self.CD_m3_tof = self.CD_m3_min_tof + self.k*self.CL_m3_tof**2 + self.CD_m3_banner_tof
            self.drag_m3_tof = self.CD_m3_tof*1/2*self.rho*self.V_m3_tof**(2)*self.S

            # cruise performance
            self.CL_m3_cruise = self.mass_m3_gross_list[i]*9.81*2/self.rho/self.V_m3_cruise**(2)/self.S
            self.CD_m3_min_cruise = self.CDmin_no_fus + self.CD_m3_fus_cruise
            # S_banner = self.l_banner**2 / self.AR_banner

            temp_CD = self.banner_cd_interp_fun(self.l_banner)

            self.CD_m3_banner_cruise = max( min(self.banner_data_cd), min(temp_CD, max(self.banner_data_cd))) * self.l_banner**(2)/self.AR_banner/self.S


            # self.CD_m3_banner_cruise = 0.02 * ((self.l_banner**(2)/self.AR_banner)/self.S) # (0.579*np.exp(-18.505*S_banner)+0.0655)*(S_banner/self.S) # 0.033*((self.l_banner**2)/(self.AR_banner))**(2)/self.S 
            self.CD_m3_cruise = self.CD_m3_min_cruise + self.k*self.CL_m3_cruise**2 + self.CD_m3_banner_cruise
            p_dyn_m3_cruise = 1/2*self.rho*self.V_m3_cruise**2
            self.drag_m3_cruise = p_dyn_m3_cruise*self.S*(self.CD_m3_min_cruise + self.CD_m3_banner_cruise + self.k*(self.CL_m3_cruise)**(2))
            # self.CD_m3_cruise*1/2*self.rho*self.V_m3_cruise**(2)*self.S_fuselage
            self.drag_m3_turn = p_dyn_m3_cruise*self.S*(self.CD_m3_min_cruise + self.CD_m3_banner_cruise + self.k*(self.n_m3*self.CL_m3_cruise)**(2))

            # Constraint analysis between take off and constant velocity level turn
            self.P_m3_tof_avail = self.calculate_take_off_power(self.mass_m3_gross_list[i], self.CD_m3_tof, self.CL_m3_tof, mu=0.05)
            P_m3_tof_batt = self.P_m3_tof_avail/self.eff_system_tof
            P_m3_turn_avail = self.drag_m3_turn * self.V_m3_cruise
            P_m3_turn_batt = P_m3_turn_avail/self.eff_system_cruise
            if P_m3_tof_batt > P_m3_turn_batt:
                """ size to take off power """
                self.P_m3_total_climb_batt = P_m3_tof_batt
                P_m3_motor = self.P_m3_tof_avail / self.eff_propeller_tof
                P_sizing_batt_m3 = P_m3_tof_batt
            else:
                """ size to turning power """
                self.P_m3_total_climb_batt = P_m3_turn_batt
                P_m3_motor = P_m3_turn_avail / self.eff_propeller_cruise
                P_sizing_batt_m3 = P_m3_turn_batt

            # Energy needed for take off
            self.time_m3_tof = self.calculate_time_to_tof(self.V_m3_tof)
            self.E_m3_tof = self.P_m3_tof_avail / self.eff_system_tof * self.time_m2_tof / 3600

            # Energy needed for climb
            P_m3_climb_avail = self.P_m3_total_climb_batt * self.eff_system_tof
            self.V_m3_v = self.calculate_climb_rate(self.mass_m3_gross_list[i], self.CD_m3_fus_lof, P_m3_climb_avail, self.V_m3_st)
            self.time_m3_climb = self.calculate_time_to_complete_climb(self.V_m3_v)
            self.E_m3_climb = self.P_m3_total_climb_batt*self.time_m3_climb/3600

            # Energy needed for cruise segment
            P_m3_level_batt = self.drag_m3_cruise*self.V_m3_cruise/self.eff_system_cruise
            self.laps_m3, self.time_m3_cruise, self.E_m3_cruise = self.calculate_cruise_performance(P_m3_level_batt, P_m3_turn_batt, self.mass_m3_gross_list[i], self.V_m3_cruise, self.time_m3_tof, self.time_m3_climb, self.phi_m3)

            # ----- Propulsion Sizing ------ #
            # Motor, sized for most power needed from M2 or M3
            if P_m3_motor > P_m2_motor:
                """ Size motor to mission three """
                self.P_all_motor.append(P_m3_motor)
                P_sizing_batt = P_sizing_batt_m3
            else:
                self.P_all_motor.append(P_m2_motor)
                P_sizing_batt = P_sizing_batt_m2
            mass_motor = self.motor_mass_from_power(self.P_all_motor[i])
            error_motor = abs(mass_motor-self.mass_motor[i])

            # Battery, sized for M2 and M3
            self.E_battery_m2_list.append(1/0.75*(self.E_m2_tof+self.E_m2_climb+self.E_m2_cruise))
            mass_m2_battery = self.battery_mass_from_energy(self.E_battery_m2_list[i])
            self.E_battery_m3_list.append(1/0.75*(self.E_m3_tof+self.E_m3_climb+self.E_m3_cruise))
            mass_m3_battery = self.battery_mass_from_energy(self.E_battery_m3_list[i])
            error_battery_m2 = abs(mass_m2_battery - self.mass_m2_battery[i])
            error_battery_m3 = abs(mass_m3_battery - self.mass_m3_battery[i])

            sizing_error_motor.append(error_motor)
            sizing_error_m2_batt.append(error_battery_m2)
            sizing_error_m3_batt.append(error_battery_m3)
            sizing_error = max(error_battery_m2, error_battery_m3, error_motor)
            if sizing_error < 0.01:
                self.mass_m1_gross = self.mass_empty_list[i]
                self.mass_m2_gross = self.mass_m2_gross_list[i]  # I may need a better way to name these things...
                self.mass_m3_gross = self.mass_m3_gross_list[i]
                self.wing_loading_m2 = self.mass_m2_gross * 9.81 / self.S
                self.wing_loading_m3 = self.mass_m3_gross * 9.81 / self.S

                self.P_m2_batt = P_sizing_batt_m3
                self.P_m2_motor = P_m2_motor
                self.P_m3_batt = P_sizing_batt_m3
                self.P_m3_motor = P_m3_motor
                self.P_sizing_batt = P_sizing_batt
                self.P_motor = self.P_all_motor[i]                
                self.power_loading_m2 = P_sizing_batt_m2 / (self.mass_m2_gross*9.81)
                self.power_loading_m3 = P_sizing_batt_m3 / (self.mass_m3_gross*9.81)

                self.E_battery_m2 = self.E_battery_m2_list[i]
                self.E_battery_m3 = self.E_battery_m3_list[i]

                break

        if debug_text is True:
            # Debugging
            print('\n$$$$$$$$$$ Debugging $$$$$$$$$$')
            eff_batt_to_motor_tof = self.eff_system_tof/self.eff_propeller_tof
            eff_batt_to_motor_cruise = self.eff_system_cruise/self.eff_propeller_cruise
            print('\n=== mission two stuffs ===')
            print(f'Drag buildup... \n> cd_fus_tof: {self.CD_m2_fus_lof}, cd_fus_cruise: {self.CD_m2_fus_cruise}')
            print(f'> cd_tof: {self.CD_m2_tof}, cd_cruise: {self.CD_m2_cruise}')
            print(f'> D_tof: {self.drag_m2_tof}, D_cruise: {self.drag_m2_cruise}, D_turn: {self.drag_m2_turn}')
            print(f'Power buildup... \n> P_motor_tof: {P_m2_tof_batt*eff_batt_to_motor_tof}, P_motor_cruise: {P_m2_level_batt*eff_batt_to_motor_cruise}, P_motor_turn: {P_m2_turn_avail/self.eff_propeller_cruise}')
            print(f'climb rate: {self.V_m2_v}')
            print(f'Course params... \n> t_tof: {self.time_m2_tof}, t_climb: {self.time_m2_climb}, t_cruise: {self.time_m2_cruise}')
            print('\n=== mission three stuff ===')
            print(f'Drag buildup... \n> cd_fus_tof: {self.CD_m3_fus_lof}, cd_fus_cruise: {self.CD_m3_fus_cruise}')
            print(f'> cd_tof: {self.CD_m3_tof}, cd_cruise: {self.CD_m3_cruise}')
            print(f'> cd_banner_tof: {self.CD_m3_banner_tof}, cd_banner_cruise: {self.CD_m3_banner_cruise}')
            print(f'> D_tof: {self.drag_m3_tof}, D_cruise: {self.drag_m3_cruise}, D_turn: {self.drag_m3_turn}')
            print(f'Power buildup... \n> P_motor_tof: {P_m3_tof_batt*eff_batt_to_motor_tof}, P_motor_cruise: {P_m3_level_batt*eff_batt_to_motor_cruise}, P_motor_turn: {P_m3_turn_avail/self.eff_propeller_cruise}')
            print(f'> climb rate: {self.V_m3_v}')
            print(f'Course params... \n> t_tof: {self.time_m3_tof}, t_climb: {self.time_m3_climb}, t_cruise: {self.time_m3_cruise}')
            
            print('\n=== Extras ===')
            print(f'> Vol_batt: {volume_battery} m^2, Vol_motor: {volume_motor} m^2')
            print(f'\n test: P_m3_level_batt; {P_m3_level_batt}, P_m3_turn_batt: {P_m3_turn_batt}, P_m3_turn_avail: {P_m3_turn_avail}')
            print('\n\n\n')

        if show_summary is True:
            """ call summary for m3 function """
            self.calculate_total_mission_score()
            self.all_sizing_summary()


    def m2_sizing_summary(self):
        """ temp way to get an okay summary if needed """
        print("\n====== Parameters =======")
        print(f"wingspan: {self.wing_wingspan} m, chord: {self.wing_chord} m")
        print(f"area: {self.wing_area} m^2, aspect ratio: {self.wing_aspect_ratio}")
        print(f"ducks: {self.ducks}, pucks: {self.cargo}, m2 laps {self.laps_m2}")
        print(f"stall speed: {self.V_m2_st} m/s, lift off speed: {self.V_m2_tof} m/s, cruise speed: {self.V_m2_cruise} m/s")
        print("======= Performance =======")
        print(f"total energy needed for m2: {self.E_battery_m2} W*hr")
        print(f"wing loading: {self.wing_loading_m2}, power loading: {self.power_loading_m2}")
        print(f"total mission score: {self.total_mission_score}, m2: {self.mission_two_score}, gm: {self.ground_mission_score}")
        plt.figure()
        plt.title('battery energy per self.iteration_m2 in W*hr')
        plt.plot(self.iteration_m2, self.E_battery_m2_list)
        plt.figure()
        plt.title('motor power per self.iteration_m2 in W')
        plt.plot(self.iteration_m2, self.P_m2_motor)
        plt.figure()
        plt.title('charter config gross mass per self.iteration_m2 in kg')
        plt.plot(self.iteration_m2, self.mass_m2_gross_list)

    def all_sizing_summary(self):
        plt.rcParams['axes.grid'] = True
        """ get summary for all sizing """
        print("====== Parameters =======")
        print(f"wingspan: {self.wing_wingspan} m, chord: {self.wing_chord} m")
        print(f"area: {self.wing_area} m^2, aspect ratio: {self.wing_aspect_ratio}")
        print(f"ducks: {self.ducks}, pucks: {self.cargo}, m2 lap: {self.laps_m2}")
        print(f"banner length: {self.l_banner}, m3 laps: {self.laps_m3} ")
        print("====== Mission two performance =======")
        print(f" Mission two score: {self.mission_two_score}, mass: {self.mass_m2_gross} kg, net_income: {self.m2_net_income}")
        print(f"total energy needed for m2: {self.E_battery_m2} W*hr")
        print(f"wing loading: {self.wing_loading_m2} N/m^2, power loading: {self.power_loading_m2} W/N")
        print(f"stall speed: {self.V_m2_st} m/s, lift off speed: {self.V_m2_tof} m/s, cruise speed: {self.V_m2_cruise} m/s")
        print("======= Mission three performance =======")
        print(f" Mission three score: {self.mission_three_score}, mass: {self.mass_m3_gross} kg, m3_cost: {self.m3_cost}")
        print(f"total energy needed for m3: {self.E_battery_m3} W*hr")
        print(f"wing loading: {self.wing_loading_m3} N/m^2, power loading: {self.power_loading_m3} W/N")
        print(f"stall speed: {self.V_m3_st} m/s, lift off speed: {self.V_m3_tof} m/s, cruise speed: {self.V_m3_cruise} m/s")
        print("======= Overall performance =======")
        print(f"motor power: {self.P_motor} W, battery power: {self.P_sizing_batt} W")
        print(f"total mission score: {self.total_mission_score}, gm: {self.ground_mission_score}\n load_factor_m2: {self.n_m2}, load_factor_m3: {self.n_m3}")
        plt.figure()
        plt.title("M2 battery energy per iteration")
        plt.plot(self.iteration_all, self.E_battery_m2_list)
        plt.figure()
        plt.title("M3 battery energy per iteration")
        plt.plot(self.iteration_all, self.E_battery_m3_list)
        plt.figure()
        plt.title('motor power per iteration in W')
        plt.plot(self.iteration_all, self.P_all_motor)
        plt.figure()
        plt.title('charter config gross mass per iteration in kg')
        plt.plot(self.iteration_all, self.mass_m2_gross_list)
        plt.figure()
        plt.title('banner config gross mass per iteration in kg')
        plt.plot(self.iteration_all, self.mass_m3_gross_list)

    def mass_buildup_summary(self):
        """ print summary of mass build up of aircraft"""