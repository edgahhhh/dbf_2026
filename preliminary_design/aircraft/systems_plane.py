""" systems plane 
    preliminary design aircraft object """
# import numpy as np
# import matplotlib.pyplot as plt
from base_aircraft import BaseAircraft

# pylint: disable=redefined-outer-name

# def drag_buildup():
#     """ function to develop a drag buildup for an aircraft model
#     this may look like something that estimates cf for the wetted area of the plane
#     and then estimate drag for components from gudmundsson 
#     This may return a np.array with drag coefficients in reference to wing area """
#     return None

# Temp stuff for testing

class Aircraft(BaseAircraft):
    """ Aircraft class for preliminary design use """
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
            altitude=altitude
            )
        
        def _mission_two_performance(self):
            """ calculate mission two performance """
            # The overall goal I THINK is to find the optimal cruise speed for mission two
            # Also want to find stuff like lap times 

        
        def calculate_performance(self,):
            """ this should calculate performance of the aircraft in all missions"""
    

    