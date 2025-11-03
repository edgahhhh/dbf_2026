""" Actuating Disk Analysis """
import numpy as np
import matplotlib.pyplot as plt

# pylint: disable=invalid-name

class ActuatingDisk():
    """ actuating disk model """
    def __init__(self, diam, alt):
        """
        @param diam: diameter [m]
        @param alt: altitude ASSL [m]
        """
        self.diameter = diam
        self.rho = self.density_from_altitude(alt)

    @ staticmethod
    def density_from_altitude(altitude ,g=9.81):
        """
        Calculate density from altitude
        """
        R = 287.05 # specific gas constant for dry air [J/kg*K]
        T_0 = 288.15 # sea level temperature [K]
        P_0 = 101325 # sea level pressure [Pa]
        L = -0.0065 # temperature lapse rate [K/m]

        T = T_0 + L * altitude # temperature at given altitude
        P = P_0 * (T/T_0)**(-g/(R*L)) # pressute at given altitude
        rho = P / (R*T) # density at given altitude

        return rho

    def power_req(self, thrust, velocity):
        """ 
        Find power required to drive propeller to meet given condition
        @param thrust [N]
        @param velocity [N]
        """
        kp = 1
        P_req = 1/kp*(thrust*velocity+thrust**(1.5)/np.sqrt(2*self.rho*np.pi/4*self.diameter**(2)))
        return P_req

diameters = 1/39.37*np.linspace(4, 30, 34)
models_list = []
p_req_list = []
for d in diameters:
    model = ActuatingDisk(diam=d, alt=400)
    models_list.append(model)
    p_temp = model.power_req(23.445, 12.511)
    p_req_list.append(p_temp)
p_req = np.asarray(p_req_list)

# Plotting the results
plt.figure(figsize=(10, 6))
plt.plot(diameters * 39.37, p_req, label='Power Required', color='blue')
plt.title('Power Required vs Propeller Diameter')
plt.xlabel('Propeller Diameter (inches)')
plt.ylabel('Power Required (W)')
plt.grid()
plt.show()
