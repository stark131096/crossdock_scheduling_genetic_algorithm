# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 14:47:44 2020

@author: StarkPC
"""

import re
import pandas as pd
import numpy as np
import math
import random
import matplotlib.pyplot as plt
# %matplotlib inline 
#!pip install mpld3
#import mpld3
#mpld3.enable_notebook() 
import datetime

class Sim_Params:
    def __init__(self):
        # demand and supply parameters
        self.num_supply      = 2         # suppliers
        self.num_cd          = 2         # cross docks
        self.num_retail      = 5         # retailers
        self.num_products    = 2         # products
        self.num_in_vehicle  = 2         # inbound vehicles in each crossdock
        self.num_out_vehicle = 2         # outbound vehicles in each crossdock
        
        # test parameters
        # have used deterministic values instead of probabilistic
        # will change the values in future
        self.lambda_p   = 0.5
        self.csp_11     = 3.5
        self.cik_1      = 3.5
        self.w_i_111    = 550
        self.w_sp_11    = 750
        self.h_ip       = 0.025
        self.o_s_111    = 55
        self.o_rr_1     = 55
        self.o_is_11    = 75
        self.o_ir_1     = 75
        self.l_iksp_1   = 0.075
        self.u_ikp_1    = 0.075
        self.l_ikp      = 0.075
        self.u_ikrp     = 0.075
        self.d_rp       = 100
        
        self.f_ik_1     = 1000
        self.f_ik       = 1000
        self.w_ik_1     = 150
        self.w_ik       = 150

class Product:
    def __init__(self):
        self._id         = -1
        self._supplier   = -1
        self._retailer   = -1
        self._cd         = -1
        self._invehicle  = -1
        self._outvehicle = -1
        self._inrank     = -1
        self._outrank    = -1
        self._picked_up  = -1
        self._delivered  = -1
        self._waittime   = -1
        
        