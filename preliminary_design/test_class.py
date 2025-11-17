import numpy as np


 
class Car():
    """ object for finding adequate props """
    def __init__(self, wheel):
        """ Initialize object for each propeller """
        self.wheel = wheel
        print(f'wheel size: {self.wheel}')
    
    def function(self, new_size):
        self.wheel = new_size

    @staticmethod
    def function(this_new_wheel):



this_car = Car(16)
this_car.function()
this_other_car=Car(18)