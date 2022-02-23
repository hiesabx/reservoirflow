#%% 1. Import Statements: 
from openresim.base import Base
import numpy as np


#%% 2. Well Class: 
class Well(Base):
    '''
    Well class to create a well.
    '''
    name = 'Well'
    def __init__(self, i=None, q=None, s=None, r=None):
        self.set_properties(i, q, s, r)

    def set_location(self, i):
        self.i = i
    set_loc = set_location


    def set_rate(self, q):
        self.rate = self.q = q
    set_q = set_rate


    def set_skin(self, s):
        self.skin = self.s = s
    set_s = set_skin


    def set_radius(self, r):
        self.radius = self.r = r
    set_r = set_radius


    def set_properties(self, i=None, q=None, s=None, r=None):
        if i:
            self.set_location(i)
        if q:
            self.set_rate(q)
        if s:
            self.set_skin(s)
        if r:
            self.set_radius(r)


if __name__ == '__main__':
    well = Well(i=4, q=-600, s=1.5, r=3.5)
    print(well)