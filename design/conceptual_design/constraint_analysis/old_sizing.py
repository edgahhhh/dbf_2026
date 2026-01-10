import numpy as np
import matplotlib.pyplot as plt

class aircraft():
    def __init__(self, aircraft_design_parameters, mission_parameters, course_parameters, initial_guess):
        """ Initialize aircraft model """
        # Aircraft Parameters
        self.AR = aircraft_design_parameters['aspect ratio']
        self.b = aircraft_design_parameters['wing span']
        self.S = self.b**2/self.AR
        self.CLmax = aircraft_design_parameters['CLmax']
        self.CDmin_no_fus = aircraft_design_parameters['CDmin_no_fuselage']
        self.k = 1/(np.pi*self.AR*0.85)
        self.fineness_ratio = aircraft_design_parameters['fineness ratio']

        # Mission parameters, instead of laps, cruising speed is used for simplication
        self.cargo = mission_parameters['cargo']
        self.ducks = mission_parameters['ducks']
        self.V_m2_cruise = mission_parameters['mission two cruise speed']

        # Atmospheric conditions
        self.alt = aircraft_design_parameters['conditions']['altitude asl']
        self.cruise_alt = aircraft_design_parameters['conditions']['cruise altitude agl']
        self.rho = self.density_from_alt(self.cruise_alt)

        # Course parameters, generally won't change
        self.ground_run = course_parameters['ground run']
        self.phi = aircraft_design_parameters['turning bank angle']
        self.time_limit = course_parameters['time limit']
        self.length_straight = course_parameters['length of straights']

        # Initial guesses
        # generic
        self.mass_motor_initial_guess = initial_guess['motor mass']  # Initial motor mass guess, kg
        self.mass_landing_gear = 0.500    # From UCLA landing gear, kg, each aircraft has this
        self.mass_ducks =  self.ducks * 0.7 / 35.274    # mass of ducks, kg,
        self.mass_cargo = self.cargo * 6 / 35.274  # mass of cargo, kg
        self.mass_wing = self.S / 3  # kg, replace w/ wing weight to wingspan ratio or something similar
        self.mass_tail = self.mass_wing / 4   # kg
        # mission two
        self.mass_m2_battery_initial_guess = initial_guess['m2 battery mass']  # Initial battery mass guess, kg

        # Systems parameters, more or less placeholders for the time being, estimating power flow in take off and cruise
        self.eff_propeller_tof = 0.4
        self.eff_motor_tof = 0.9
        self.eff_esc_tof = 0.9
        self.eff_battery_tof = 0.95
        self.eff_system_tof = self.eff_propeller_tof * self.eff_motor_tof * self.eff_esc_tof * self.eff_battery_tof
        self.eff_propeller_cruise = 0.5
        self.eff_motor_cruise = 0.9
        self.eff_esc_cruise = 0.9
        self.eff_battery_cruise = 0.95
        self.eff_system_cruise = self.eff_propeller_cruise * self.eff_motor_cruise * self.eff_esc_cruise * self.eff_battery_cruise

        # Sizing values
        self.iteration = []
        self.mass_m2_gross_list = []
        self.P_m2_motor = []    # motor sized to take off power
        self.E_battery_m2_list = []

        # Initialize variables that are mainly going to be used once the aircraft is sized
            # What matters here is the actual final parameters of aircraft, can have a function that populates this after sizing
        # performance
        self.wing_loading_m2 = 0.0
        self.power_loading_m2 = 0.0
        self.wing_loading_m3 = 0.0
        self.power_loading_m3 = 0.0 

        # Weight buildup
        self.gross_mass_m1 = 0.0
        self.mass_m2_gross = 0.0
        self.gross_mass_m3 = 0.0

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
        self.motor_power = 0.0  # constant, can't change the motor
        self.motor_mass = 0.0   # same as before
        self.battery_energy_m1 = 0.0
        self.battery_mass_m1 = 0.0
        self.E_battery_m2 = 0.0
        self.battery_mass_m2 = 0.0
        self.battery_energy_m3 = 0.0
        self.battery_mass_m3 = 0.0

        # Mission Two Performance
        self.V_m2_st = 0.0
        self.V_m2_tof = 0.0
        self.V_m2_v = 0.0
        self.CD_m2_fus_lof = 0.0
        self.CD_m2_fus_cruise = 0.0
        self.CL_m2_tof_m2 = 0.0
        self.CL_m2_cruise = 0.0
        self.CD_m2_min_tof = 0.0
        self.CD_m2_min_cruise = 0.0
        self.CD_m2_tof = 0.0
        self.CD_m2_cruise = 0.0
        self.P_m2_tof_avail = 0.0 # T*V at take off roll
        self.P_m2_total_climb = 0.0   # Total power at battery during climb 
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
        self.mission_three_score = 0.0
        self.ground_mission_score = 0.0
        self.total_mission_score = 0.0
        self.banner_length = 0.0

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

    def calculate_fuselage_drag(self, reynolds_number, S_wet):
        """ calculate profile drag of fuselage using Horner estimation for streamline shapes """
        cf_fuselage = 1.328 / np.sqrt(reynolds_number)
        CD_fuselage_wet = cf_fuselage*(1+(1/self.fineness_ratio)**(3/2)) + 0.11*(1/self.fineness_ratio)**2
        CD_fuselage= CD_fuselage_wet * S_wet / self.S
        return CD_fuselage

    def battery_mass_from_energy(self, battery_energy):
        """ estimate mass of battery from total energy capacitance 
        This can be done two ways
        1. energy density constant
        2. energy/mass density where it increases as energy increase
        """
        energy_density = 0.01952 # kg/Whr
        return (abs(battery_energy*energy_density))

    def battery_volume_from_mass(self, battery_mass):
        """ Calculate battery volume from mass """
        battery_mass_density = 1988.915 # kg/m^3
        return(battery_mass / battery_mass_density)

    def motor_mass_from_power(self, max_continuous_shaft_power):
        """ estimate motor mass from power required to drive shaft/propeller at take off
        assuming linear increase in power and motor mass, could use a better model (FRESHMAN!!!!!) 
        """
        motor_power_loading = 5.6   # W/g
        return (abs(max_continuous_shaft_power / motor_power_loading / 1000))

    def motor_volume_from_mass(self, motor_mass):
        """ Calculate battery mass from power"""
        motor_mass_density = 68.45  # kg/m^3
        return (motor_mass/motor_mass_density)

    def calculate_take_off_power(self, mass, CDtof, CLtof, mu=0.05, gravity = 9.81):
        """ estimate power at battery to achieve take off, assuming a constant acceleration """
        wing_loading = mass*gravity/self.S 
        thrust = (1.21/gravity/self.rho/self.CLmax/self.ground_run*wing_loading + 0.605/self.CLmax*(CDtof-mu*CLtof) + mu)*mass*gravity
        available_power = thrust * np.sqrt(2*wing_loading/self.rho/self.CLmax)
        return available_power  # Watts, if everything is in SI

    def calculate_time_to_tof(self, Vlof):
        """ Calculate time to take off """
        acceleration = Vlof**2 / 2 / self.ground_run # m/s^2
        return (Vlof/acceleration)  # seconds

    def calculate_climb_rate(self, total_power, available_power, mass, gravity=9.81):
        excess_power = total_power - available_power
        return (excess_power/mass/gravity)  # climb rate, m/s?

    def calculate_time_to_complete_climb(self, climb_rate):
        """ calculate time to finish climb sequence """
        return (self.cruise_alt / climb_rate)

    def calculate_power_needed_for_level_turn(self, mass, airspeed, CDmin, gravity=9.81):
        """ calculate available power needed to complete a level turn at some angle and airspeed """
        wing_loading = mass*gravity/self.S
        load_factor = 1/np.cos(self.phi)
        p_dynamic = 1/2*self.rho*airspeed**2
        thrust = mass*gravity*p_dynamic*(CDmin/wing_loading + self.k*((load_factor/p_dynamic)**2) * wing_loading)
        return(thrust*airspeed)

    def calculate_power_needed_for_cruise(self, mass, airspeed, CDmin, gravity=9.81):
        """ calculate available power needed to cruise at steady state (ignore trim drag) """
        wing_loading = mass*gravity/self.S
        p_dynamic = 1/2*self.rho*airspeed**2
        thrust = mass*gravity*(p_dynamic*CDmin/wing_loading + self.k/p_dynamic*wing_loading)
        return(thrust*airspeed)

    def calculate_cruise_performance(self, level_power_batt, turn_power_batt, mass, cruise_speed, time_tof, time_climb, gravity=9.81):
        """ 
        calculate laps, time, and energy needed for cruise 
        @return laps
        @return time_cruise_total
        @return E_cruise
        """
        radius_m2_turn = mass * cruise_speed**2 / (1/np.cos(self.phi)*mass*gravity)/np.sin(self.phi)
        time_turn_ea = radius_m2_turn/cruise_speed*(180*np.pi/180)
        time_circle_ea = time_turn_ea*2
        time_straights_ea = self.length_straight/cruise_speed   # same for each unique cruising speed
        time_lap_1_temp = time_tof + time_climb + 2*time_turn_ea + time_straights_ea + time_circle_ea
        time_after_lap_1 = self.time_limit - time_climb - time_tof - time_lap_1_temp  # How much time left after first lap in seconds
        time_ea_lap = 2*time_straights_ea + 2*time_turn_ea + time_circle_ea
        laps = np.floor(1+time_after_lap_1/time_ea_lap)
        n_turns = laps*4
        n_straights = laps*2 - 1
        time_cruise_total = time_after_lap_1 + time_ea_lap*(laps-1)
        E_cruise = (turn_power_batt*n_turns*time_turn_ea + level_power_batt*n_straights*time_straights_ea)/3600
        return laps, time_cruise_total, E_cruise

    def calculate_ground_mission_score(self):
        """" calculate ground mission score """
        t_gm = self.ducks*0.5 + np.ceil(self.cargo/10)   # not accounting for cargo
        self.ground_mission_score = max(0, min(20/t_gm, 1))
        return(self.ground_mission_score)

    def calculate_mission_two_score(self):
        """ calculate total mission two score """
        income = self.ducks*(6 + 2*self.laps_m2) + self.cargo*(10 + 8*self.laps_m2)
        cost = self.E_battery_m2/100*self.laps_m2*(10+self.ducks*0.5+self.cargo*2)
        # print(income-cost)
        self.mission_two_score = 1 + max(0, min((income-cost)/2500, 1))
        return(self.mission_two_score)

    def calculate_mission_three_score(self):
        """ calculate mission three score """
        

    def calculate_total_mission_score(self):
        """ wrapper to calculate total mission score """
        self.total_mission_score = 1 + self.calculate_mission_two_score() + self.calculate_ground_mission_score() + 2.5

    def size_aircraft_mission_two(self):
        """ Size aircraft for mission two """
        # Initialize some lists for calculations (maybe for plotting as well?), maybe won't need these all?
        sizing_error = 1
        mass_m2_battery = []
        mass_motor = []
        sizing_error_batt = []
        sizing_error_motor = []
        for i in range(50):
            self.iteration.append(i)
            if i == 0:
                """ initial guess """
                mass_m2_battery.append(self.mass_m2_battery_initial_guess)  # Initial battery mass guess, kg
                mass_motor.append(self.mass_motor_initial_guess)  # Initial motor mass guess, kg
            elif i > 0:
                """ refined mass build up """
                mass_m2_battery.append(self.battery_mass_from_energy(self.E_battery_m2_list[i-1]))
                mass_motor.append(self.motor_mass_from_power(self.P_m2_motor[i-1]))
            # weight estimation
            self.mass_m2_gross_list.append(mass_m2_battery[i] + mass_motor[i] + self.mass_landing_gear + self.mass_ducks + self.mass_cargo + self.mass_wing + self.mass_tail)  # kg
            # Estimate stalling speed and take off speed from mass
            self.V_m2_st = np.sqrt(2*self.mass_m2_gross_list[i]*9.81/self.rho/self.S/self.CLmax)
            self.V_m2_tof = 1.1*self.V_m2_st

            # ----- Fuselage Sizing ----- #
            volume_m2_battery = self.battery_volume_from_mass(mass_m2_battery[i])
            volume_motor = self.motor_volume_from_mass(mass_motor[i])
            volume_m2_electronics = (volume_m2_battery + volume_motor)*1.1
            volume_m2_payload = (self.ducks*(2.5/39.37)**3 + self.cargo * np.pi/4/39.37) # m^3
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

            # ---- Estimate aerodynamics ----- #
            # take off configuration
            self.CL_m2_tof = self.mass_m2_gross_list[i]*9.81*2/self.rho/self.V_m2_tof**2/self.S
            self.CD_m2_min_tof = self.CDmin_no_fus + self.CD_m2_fus_lof
            self.CD_m2_tof = self.CD_m2_min_tof + self.k*self.CL_m2_tof**2
            # cruise configuration
            self.CL_m2_cruise = self.mass_m2_gross_list[i]*9.81*2/self.rho/self.V_m2_cruise**2/self.S
            self.CD_m2_min_cruise = self.CDmin_no_fus + self.CD_m2_fus_cruise
            self.CD_m2_cruise = self.CD_m2_min_cruise + self.k*self.CL_m2_cruise**2

            # ------ Estimate performance ----- #ADD A FUNCTION TO CALCULATE POWER AT ALL POINTS OF THE SYSTEM, so I dont populate this too much
            # Take off
            self.P_m2_tof_avail = self.calculate_take_off_power(self.mass_m2_gross_list[i], self.CD_m2_tof, self.CL_m2_tof, mu=0.05)
            self.time_m2_tof = self.calculate_time_to_tof(self.V_m2_tof)
            self.E_m2_tof = self.P_m2_tof_avail/self.eff_system_tof*self.time_m2_tof/3600
            # Climb
            self.P_m2_total_climb = self.P_m2_tof_avail/self.eff_system_tof # power at battery at full throttle
            self.P_m2_avail_climb = self.CD_m2_tof/2*self.rho*self.V_m2_tof**3*self.S/self.eff_system_tof   # review this, maybe make it a function so its more clearly documented
            self.P_m2_exc_climb = self.P_m2_total_climb-self.P_m2_avail_climb  # Not critical, but nice to store
            self.V_m2_v = self.calculate_climb_rate(self.P_m2_total_climb, self.P_m2_avail_climb, self.mass_m2_gross_list[i]) # m/s
            self.time_m2_climb = self.calculate_time_to_complete_climb(self.V_m2_v)
            self.E_m2_climb = self.P_m2_total_climb*self.time_m2_climb/3600
            # Cruise 
            P_m2_turn_batt = self.calculate_power_needed_for_level_turn(self.mass_m2_gross_list[i], self.V_m2_cruise, self.CD_m2_min_cruise)/self.eff_system_cruise
            P_m2_level_batt = self.calculate_power_needed_for_cruise(self.mass_m2_gross_list[i], self.V_m2_cruise, self.CD_m2_min_cruise)/self.eff_system_cruise
            self.laps_m2, self.time_m2_cruise, self.E_m2_cruise = self.calculate_cruise_performance(P_m2_level_batt, P_m2_turn_batt, self.mass_m2_gross_list[i], self.V_m2_cruise, self.time_m2_tof, self.time_m2_climb)

            # ----- Propulsion Sizing ------ #
            # Motor, sized for most power needed from M2 or M3
            self.P_m2_motor.append(self.P_m2_tof_avail/self.eff_propeller_tof)
            mass_motor_temp = self.motor_mass_from_power(self.P_m2_motor[i])
            error_motor = abs(mass_motor_temp-mass_motor[i])
            # Battery, sized for M2 and M3
            self.E_battery_m2_list.append(1.3*(self.E_m2_tof+self.E_m2_climb+self.E_m2_cruise))
            mass_battery_temp = self.battery_mass_from_energy(self.E_battery_m2_list[i])

            error_battery = abs(mass_battery_temp-mass_m2_battery[i])
            sizing_error_batt.append(error_battery)
            sizing_error_motor.append(error_motor)
            sizing_error = max(error_battery, error_motor)
            if sizing_error < 0.001:
                self.mass_m2_gross = self.mass_m2_gross_list[i]  # I may need a better way to name these things...
                self.wing_loading_m2 = self.mass_m2_gross_list[i] * 9.81 / self.S
                self.power_loading_m2 = self.P_m2_tof_avail / self.eff_system_tof / (self.mass_m2_gross*9.81)
                self.E_battery_m2 = self.E_battery_m2_list[i]
                self.n_turn = 1/np.cos(self.phi)
                break

    def size_aircraft_all_mission(self, roc=10):
        """ Size aircraft to meet parameters for mission two and mission three. 
        philosophy: 
            motor: sized to meet greatest power demanded 
            battery: sized for each mission
            
            """



    def m2_sizing_summary(self):
        """ temp way to get an okay summary if needed """
        print("====== Parameters =======")
        print(f"wingspan: {self.wing_wingspan} m, chord: {self.wing_chord} m")
        print(f"area: {self.wing_area} m^2, aspect ratio: {self.wing_aspect_ratio}")
        print(f"ducks: {self.ducks}, pucks: {self.cargo}, m2 laps {self.laps_m2}")
        print(f"stall speed: {self.V_m2_st} m/s, lift off speed: {self.V_m2_tof} m/s, cruise speed: {self.V_m2_cruise} m/s")
        print("======= Performance =======")
        print(f"total energy needed for m2: {self.E_battery_m2} W*hr")
        print(f"wing loading: {self.wing_loading_m2}, power loading: {self.power_loading_m2}")
        print(f"total mission score: {self.total_mission_score}, m2: {self.mission_two_score}, gm: {self.ground_mission_score}")
        plt.figure()
        plt.title('battery energy per self.iteration in W*hr')
        plt.plot(self.iteration, self.E_battery_m2_list)
        plt.figure()
        plt.title('motor power per self.iteration in W')
        plt.plot(self.iteration, self.P_m2_motor)
        plt.figure()
        plt.title('charter config gross mass per self.iteration in kg')
        plt.plot(self.iteration, self.mass_m2_gross_list)