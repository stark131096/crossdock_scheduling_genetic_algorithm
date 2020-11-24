# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 21:07:40 2020

@author: StarkPC
"""

import pandas as pd
import numpy as np
import math
import random
import matplotlib.pyplot as plt
import mpld3
# from gurobipy import Model, GRB, quicksum
import collections
from collections import Counter
# import parameters
import time

start_time = time.time()
random.seed(210)

class Simulation:

    def Module_Initialization(self):
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
        
    ###############################################################################
    def Module_Create_Dataframes(self):
        # creating list of suppliers
        self.list_suppliers = list(range(1, self.num_supply + 1))
        self.dict_supplier_product_id = {new_list: [] for new_list in range(1, self.num_supply + 1)}
                
        # creating dictionary for crossdock-product assignment
        self.list_crossdocks = list(range(1, self.num_cd + 1))
        
        
        # creating dictionary for invehicle-product id assignment
        self.list_invehicles = list(range(1, self.num_in_vehicle*self.num_cd + 1))
        
        
        # creating product id (1,2,3,4,....) vs product type (1,2) dictionary
        # {1: 1, 2: 2, 3: 1, 4: 2, 5: 1, 6: 2, 7: 1, 8: 2, 9: 1, 10: 2}
        self.list_product_id = list(range(1, self.num_products * self.num_retail+1))
        self.list_type = list(range(1, self.num_products + 1)) * self.num_retail
        self.dict_prod_type = dict(zip(self.list_product_id, self.list_type))
        print (self.dict_prod_type)
        
        # creating a retailer-product id dictionary (key - retailer, value - product id)
        # {1: [1, 2], 2: [3, 4], 3: [5, 6], 4: [7, 8], 5: [9, 10]}
        self.list_retailers = list(range(1, self.num_retail + 1))
        self.grouped_product_ids = [self.list_product_id[i:i + self.num_products] for i in 
                               range(0, len(self.list_product_id), self.num_products)]       
        self.dict_retail_prod_id = dict(zip(self.list_retailers, self.grouped_product_ids))
        print(self.dict_retail_prod_id)
                
    ###############################################################################    
    def Module_Gene_Decode(self, gene_sequence):
        pass
        self.seq_1 = gene_sequence[0]   # retailer-supplier product assignment
        self.seq_2 = gene_sequence[1]   # assignment of products to cross-docks & inbound routing
        self.seq_3 = gene_sequence[2]   # consolidation decisions & outbound routing
        self.seq_4 = gene_sequence[3]   # departure times of inbound vehicles from cross-docks
    
    ###############################################################################
    def Module_Supplier_Assignment(self):
        # creating a supplier-product id dictionary (key - supplier, value - product id) from seq_1
        # {1: [1, 4, 5], 2: [2, 3, 6]}
        l1 = self.seq_1
        for i in range(len(l1)):
            for s_id in self.list_suppliers:
                if (l1[i] == s_id):
                    self.dict_supplier_product_id[s_id].append(i+1)
        print (self.dict_supplier_product_id)
        
        
    
    ###############################################################################
    def Module_Calculate_Cost(self, gene_sequence):
        # Initialize all the parameters
        self.Module_Initialization()
        self.Module_Create_Dataframes()
        self.Module_Gene_Decode(gene_sequence)
        self.Module_Supplier_Assignment()
        
gene_seq = [[1,2,2,1,1,2,1,1,2,1],2,3,4]

sim = Simulation()
sim.Module_Calculate_Cost(gene_seq)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        