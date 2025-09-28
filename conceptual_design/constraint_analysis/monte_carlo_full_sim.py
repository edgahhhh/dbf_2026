from sizing import aircraft
import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class MonteCarlo:
    def __init__(self, samples):
        self.samples = samples
        self.converged_parameters = []
        self.parameters = []
        self.wingspan = []
        self.wing_area = []
        self.aspect_ratio = []
        self.ducks = []
        self.cargo = []
        self.cruise_speed_m2 = []
        self.laps_m2 =[]
        self.l_banner = []
        self.cruise_speed_m3 = []
        self.laps_m3 = []
        self.aircraft_list = []
        self.mission_two_score = []
        self.mission_three_score = []
        self.ground_mission_score = []
        self.mission_two_and_ground_mission_score = []
        self.total_mission_score = []
        self.power_motor = []
        self.battery_energy_m2 = []
        self.wing_loading_m2 = []
        self.power_loading_m2 = []
        self.mass_m2_gross = []
        self.battery_energy_m3 = []
        self.wing_loading_m3 = []
        self.power_loading_m3 = []
        self.mass_m3_gross = []
        self.m2_net_income = []
        self.m3_cost = []

    @staticmethod
    def generate_monte_carlo_inputs():
        """ return: [wingspan, aspect_ratio, cargo, ducks, m2_cruise_speed, banner_length, m3_cruise_speed]"""
        wing_span = random.uniform(0.9144, 1.524)
        # wing_span = 1.524   # m
        aspect_ratio = random.uniform(3, 12)
        cargo = np.floor(random.uniform(1, 30))
        ducks = np.floor(random.uniform(3*cargo, 120))
        m2_cruise_speed = random.uniform(10, 40)
        l_banner = random.uniform(0.3, 12)
        m3_cruise_speed = random.uniform(10, 40)
        return (np.array([wing_span, aspect_ratio, cargo, ducks, m2_cruise_speed, l_banner, m3_cruise_speed]))
    
    def run_simulation(self):
        """ run monte carlo simulation """
        conditions = { 'altitude asl':397,'cruise altitude agl':10}
        aircraft_design_parameters = {'aspect ratio' : 4,'wing span':1.5,'CLmax':0.9,'CDmin_no_fuselage':0.03,'fineness ratio':6,'turning bank angle':np.pi/4}
        aircraft_design_parameters['conditions'] = conditions
        mission_parameters = {'cargo':15,'ducks':45,'mission two cruise speed': 25,'banner length':12,'mission three cruise speed':20,'banner aspect ratio':5}
        course_parameters = {'ground run': 30,'time limit': 300,'length of straights': 304.8 }
        initial_guess = {'motor mass': 0.1,'m2 battery mass': 0.8,'m3 battery mass':0.5}
        for i in range(self.samples):
            """ Run monte carlo """
            parameters_temp = self.generate_monte_carlo_inputs()
            self.parameters.append(parameters_temp)

            aircraft_design_parameters['wing span'] = self.parameters[i][0]
            aircraft_design_parameters['aspect ratio'] = self.parameters[i][1]
            mission_parameters['cargo'] = self.parameters[i][2]
            mission_parameters['ducks'] = self.parameters[i][3]
            mission_parameters['mission two cruise speed'] = self.parameters[i][4]
            mission_parameters['banner length'] = self.parameters[i][5]
            mission_parameters['mission three cruise speed'] = self.parameters[i][6]

            temp_aircraft = aircraft(aircraft_design_parameters, mission_parameters, course_parameters, initial_guess)
            with np.errstate(over='ignore', invalid='ignore'):  # don't log errors
                temp_aircraft.size_aircraft_all_missions()
            if not np.isfinite(temp_aircraft.V_m2_st) or temp_aircraft.E_battery_m2>99.99 or temp_aircraft.E_battery_m3>99.99:  # throw out models that don't converge or exceed 100 Whr
                continue
            elif temp_aircraft.mass_m2_gross<0.001 or temp_aircraft.mass_m3_gross<0.001:
                continue
            elif temp_aircraft.P_motor > 4000:
                continue
            else: 
                temp_aircraft.calculate_total_mission_score()
                self.converged_parameters.append(parameters_temp)   # for matrix sorting if necessary
                self.wingspan.append(self.parameters[i][0])
                self.wing_area.append(temp_aircraft.wing_area)
                self.aspect_ratio.append(self.parameters[i][1])
                self.cargo.append(self.parameters[i][2])
                self.ducks.append(self.parameters[i][3])
                self.cruise_speed_m2.append(self.parameters[i][4])
                self.l_banner.append(self.parameters[i][5])
                self.cruise_speed_m3.append(self.parameters[i][6])

                self.laps_m2.append(temp_aircraft.laps_m2)
                self.laps_m3.append(temp_aircraft.laps_m3)
                self.aircraft_list.append(temp_aircraft)
                self.mission_two_score.append(temp_aircraft.mission_two_score)
                self.mission_three_score.append(temp_aircraft.mission_three_score)
                self.ground_mission_score.append(temp_aircraft.ground_mission_score)
                self.mission_two_and_ground_mission_score.append(temp_aircraft.mission_two_score + temp_aircraft.ground_mission_score)
                self.total_mission_score.append(temp_aircraft.total_mission_score)

                self.power_motor.append(temp_aircraft.P_motor)
                self.battery_energy_m2.append(temp_aircraft.E_battery_m2)
                self.battery_energy_m3.append(temp_aircraft.E_battery_m3)

                self.wing_loading_m2.append(temp_aircraft.wing_loading_m2)
                self.wing_loading_m3.append(temp_aircraft.wing_loading_m3)

                self.power_loading_m2.append(temp_aircraft.power_loading_m2)
                self.power_loading_m3.append(temp_aircraft.power_loading_m3)

                self.mass_m2_gross.append(temp_aircraft.mass_m2_gross)
                self.mass_m3_gross.append(temp_aircraft.mass_m3_gross)

                self.m2_net_income.append(temp_aircraft.m2_net_income)
                self.m3_cost.append(temp_aircraft.m3_cost)

    def plot_m2_cargo_ducks_aspd(self):
        """ plot m2 scoring by cargo, ducks, and cruise speed """
        plt.figure()
        cargo_np = np.asarray(self.cargo)
        ducks_np = np.asarray(self.ducks)
        m2_np = np.asarray(self.mission_two_score)
        vc_m2 = np.asarray(self.cruise_speed_m2)
        ax = plt.gcf().add_subplot(111, projection='3d')
        ax.set_title('M2 score by cargo, ducks, and Vcruise')
        sc = ax.scatter(ducks_np, cargo_np, vc_m2, c=m2_np, cmap='magma')
        cb = plt.colorbar(sc, ax=ax, shrink=0.6, location='left')
        cb.set_label('m2')
        ax.set_xlabel('ducks')
        ax.set_ylabel('cargo')
        ax.set_zlabel('cruise speed (m/s)')

    def plot_m2_gm_ducks_per_cargo_laps_ws(self):
        """ plot m2 and gm score by ducks/cargo, laps, and wing loading """
        plt.figure()
        cargo_np = np.asarray(self.cargo)
        ducks_np = np.asarray(self.ducks)
        ducks_cargo_np = ducks_np/cargo_np
        laps_m2_np = np.asarray(self.laps_m2)
        m2_gm_np = np.asarray(self.mission_two_and_ground_mission_score)
        ws_np = np.asarray(self.wing_loading_m2)
        ax = plt.gcf().add_subplot(111, projection='3d')
        ax.set_title('M2 and GM score by ducks/cargo, laps, and W/S')
        sc = ax.scatter(ducks_cargo_np, laps_m2_np, ws_np, c=m2_gm_np, cmap='magma')
        cb = plt.colorbar(sc, ax=ax, shrink=0.6, location='left')
        cb.set_label('m2 + gm')
        ax.set_xlabel('ducks/cargo')
        ax.set_ylabel('laps m2')
        ax.set_zlabel('W/S (N/m^2)')

    def plot_m2_gm_ducks_per_cargo_ebatt_ws(self):
        """ plot m2 and gm score by ducks/cargo, laps, and wing loading """
        plt.figure()
        cargo_np = np.asarray(self.cargo)
        ducks_np = np.asarray(self.ducks)
        ducks_cargo_np = ducks_np/cargo_np
        e_np = np.asarray(self.battery_energy_m2)
        m2_gm_np = np.asarray(self.mission_two_and_ground_mission_score)
        ws_np = np.asarray(self.wing_loading_m2)
        ax = plt.gcf().add_subplot(111, projection='3d')
        ax.set_title('M2 and GM score by ducks/cargo, ebatt m2, and W/S')
        sc = ax.scatter(ducks_cargo_np, e_np, ws_np, c=m2_gm_np, cmap='magma')
        cb = plt.colorbar(sc, ax=ax, shrink=0.6, location='left')
        cb.set_label('m2 + gm')
        ax.set_xlabel('ducks/cargo')
        ax.set_ylabel('Ebatt m2')
        ax.set_zlabel('W/S (N/m^2)')

    def plot_m2_gm_cargo_ducks_aspd(self):
        """ plot m2 and gm score by cargo, ducks, and cruise speed """
        plt.figure()
        cargo_np = np.asarray(self.cargo)
        ducks_np = np.asarray(self.ducks)
        m2_gm_np = np.asarray(self.mission_two_and_ground_mission_score)
        vc_np = np.asarray(self.cruise_speed_m2)
        ax = plt.gcf().add_subplot(111, projection='3d')
        ax.set_title('M2 and GM score by cargo, ducks, and Vcruise')
        sc = ax.scatter(cargo_np, ducks_np, vc_np, c=m2_gm_np, cmap='magma')
        ax.set_xlabel('cargo')
        ax.set_ylabel('ducks')
        ax.set_zlabel('cruise speed (m/s)')
        cb = plt.colorbar(sc, ax=ax, shrink=0.6, location='left')
        cb.set_label('m2 + gm')

    def plot_m2_gm_pw_ws_ebatt(self):
        """ plot m2 and gm score by power loading, wing loading, and battery energy """
        plt.figure()
        ws_np = np.asarray(self.wing_loading_m2)
        pw_np = np.asarray(self.power_loading_m2)
        m2_gm_np = np.asarray(self.mission_two_and_ground_mission_score)
        e_np = np.asarray(self.battery_energy_m2)
        ax = plt.gcf().add_subplot(111, projection='3d')
        ax.set_title('M2 and GM score by m2 P/W, m2 W/S, and m2 Ebatt')
        sc = ax.scatter(ws_np, pw_np, e_np, c=m2_gm_np, cmap='magma')
        cb = plt.colorbar(sc, ax=ax, shrink=0.6, location='left')
        cb.set_label('m2 + gm')
        ax.set_xlabel('WS m2 (N/m^2)')
        ax.set_ylabel('PW m2(W/N)')
        ax.set_zlabel('Ebatt m2 (W*hr)')

    def plot_m2_gm_b_ar_mass(self):
        plt.figure()
        b_np = np.asarray(self.wingspan)
        ar_np = np.asarray(self.aspect_ratio)
        m2_gm_np = np.asarray(self.mission_two_and_ground_mission_score)
        m_np = np.asarray(self.gross_mass_m2)
        ax = plt.gcf().add_subplot(111, projection='3d')
        sc = ax.scatter(b_np, ar_np, m_np, c=m2_gm_np, cmap='magma')
        cb = plt.colorbar(sc, ax=ax, shrink=0.6, location='left')
        ax.set_title('M2 and GM score by wingspan, aspect ratio, and mass')
        cb.set_label('m2 + gm')
        ax.set_xlabel('wingspan (m)')
        ax.set_ylabel('AR')
        ax.set_zlabel('m2 gross mass (kg)')

    def plot_m2_gm_ws_ar_ebatt(self):
        plt.figure()
        ws_np = np.asarray(self.wing_loading_m2)
        ar_np = np.asarray(self.aspect_ratio)
        m2_gm_np = np.asarray(self.mission_two_and_ground_mission_score)
        e_np = np.asarray(self.battery_energy_m2)
        ax = plt.gcf().add_subplot(111, projection='3d')
        sc = ax.scatter(ws_np, ar_np, e_np, c=m2_gm_np, cmap='magma')
        cb = plt.colorbar(sc, ax=ax, shrink=0.6, location='left')
        cb.set_label('m2 + gm')
        ax.set_title('M2 and GM score by W/S, AR, m2 Ebatt')
        ax.set_xlabel('WS (N/m^2)')
        ax.set_ylabel('AR')
        ax.set_zlabel('Ebatt m2 (W*hr)')

    def plot_m3_lbanner_ar_m3speed(self):
        plt.figure()
        m3 = np.asarray(self.mission_three_score)
        l = np.asarray(self.l_banner)
        ar = np.asarray(self.aspect_ratio)
        v = np.asarray(self.cruise_speed_m3)
        ax = plt.gcf().add_subplot(111, projection='3d')
        sc = ax.scatter(l, ar, v, c=m3, cmap='magma')
        cb = plt.colorbar(sc, ax=ax, shrink=0.6, location='left')
        cb.set_label('m3 score')
        ax.set_title('M3 Score by l_banner, AR, and V_cruise_m3')
        ax.set_xlabel('l_banner (m)')
        ax.set_ylabel('AR')
        ax.set_zlabel('Vcruise_m3 (m/s)')

    def plot_total_score_lbanner_ar_m3laps(self):
        plt.figure()
        m = np.asarray(self.total_mission_score)
        length = np.asarray(self.l_banner)
        ar = np.asarray(self.aspect_ratio)
        laps = np.asarray(self.laps_m3)
        ax = plt.gcf().add_subplot(111, projection='3d')
        sc = ax.scatter(length, ar, laps, c=m, cmap='magma')
        cb = plt.colorbar(sc, ax=ax, shrink=0.6, location='left')
        ax.set_title('Total Score by l_banner, AR, and V_cruise_m3')
        ax.set_xlabel('l_banner (m)')
        ax.set_ylabel('AR')
        ax.set_zlabel('Vcruise_m3 (m/s)')

    def plot_total_score_lbanner_ar_ebatt(self):
        plt.figure()
        m = np.asarray(self.total_mission_score)
        length = np.asarray(self.l_banner)
        ar = np.asarray(self.aspect_ratio)
        e = np.asarray(self.battery_energy_m3)
        ax = plt.gcf().add_subplot(111, projection='3d')
        sc = ax.scatter(length, ar, e, c=m, cmap='magma')
        cb = plt.colorbar(sc, ax=ax, shrink=0.6, location='left')
        ax.set_title('Mission Score by l_banner, AR, and m3 Ebatt')
        ax.set_xlabel('l_banner (m)')
        ax.set_ylabel('AR')
        ax.set_zlabel('Ebatt m3 (W*hr)')

    def plot_total_score_em2_em3_pm(self):
        plt.figure()
        tm = np.asarray(self.total_mission_score)
        e2 = np.asarray(self.battery_energy_m2)
        e3 = np.asarray(self.battery_energy_m3)
        pm = np.asarray(self.power_motor)
        ax = plt.gcf().add_subplot(111, projection='3d')
        sc = ax.scatter(e2, e3, pm, c=tm, cmap='magma')
        cb = plt.colorbar(sc, ax=ax, shrink=0.6, location='left')
        cb.set_label('mission score')
        ax.set_title('Mission Score by E_m2, E_m3, and P_motor')
        ax.set_xlabel('M2, W*hr')
        ax.set_ylabel('M3, W*hr')
        ax.set_zlabel('motor, W')


plt.rcParams['axes.grid'] = True

samples = 20000
sim = MonteCarlo(samples)
sim.run_simulation()
print(f'\n\n\n\nconverged: {len(sim.converged_parameters)}')

print(f"max m2, net income {max(sim.m2_net_income)}")
print(f"max m3, cost {max(sim.m3_cost)}")
print('\n\n')

max_score = max(sim.total_mission_score)
best_index = sim.total_mission_score.index(max_score)
best_plane = sim.aircraft_list[best_index]
best_plane.all_sizing_summary()

# sim.plot_m2_cargo_ducks_aspd()
sim.plot_m2_gm_cargo_ducks_aspd()
sim.plot_m2_gm_pw_ws_ebatt()
sim.plot_total_score_em2_em3_pm()
sim.plot_m3_lbanner_ar_m3speed()
sim.plot_total_score_lbanner_ar_ebatt()
sim.plot_total_score_em2_em3_pm()
plt.show()

