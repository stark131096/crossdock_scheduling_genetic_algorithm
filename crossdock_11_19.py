# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 15:02:46 2020

@author: StarkPC
"""

import pandas as pd
import numpy as np
import math
import random
import matplotlib.pyplot as plt
import mpld3
import collections
from collections import Counter
import time

start_time = time.time()
random.seed(210)

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
        
class Crossdock_In_Vehicle:
    def __init__(self):
        self._id = -1
        self._wait_time = -1
        self._cd_id = -1
        self._route = [0]
        self._prod_order = [0]
    
class Crossdock_Out_Vehicle:
    def __init__(self):
        self._id = -1
        self._wait_time = -1
        self._cd_id = -1
        self._route = [0]
        self._prod_order = [0]
        
class Simulation:
    #product = Product()
    
    def Module_Initialization(self):
        
        self.supplier_dist = [[0,46],[46,0]]
        self.cd_supplier_dist = [[44,25],[71,11]]
        self.cd_retailer_dist = [[53,22,39,68,35],[86,79,80,55,98]]
        self.retailer_dist = [[0,21,16,21,47],[21,0,36,35,25],[16,36,0,12,81],
                              [21,35,12,0,76],[47,25,81,76,0]]
        
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
        self.cik        = 3.5 
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
        self.list_product_ids = list(range(1, self.num_retail * self.num_products + 1 ))
        self.list_products = []
        self.list_in_vehicles = []
        self.list_out_vehicles = []
        
    ###############################################################################    
    def Module_Create_Products(self):
        # total products = retailers * products
        for i in range(len(self.list_product_ids)):
            new_prod = Product()
            new_prod._id = i+1
            new_prod._retailer = i//self.num_products + 1
            self.list_products.append(new_prod)
            
    ###############################################################################
    def Module_Create_Vehicles(self):
        self.total_in_vehicles = self.num_cd*self.num_in_vehicle
        for i in range(self.total_in_vehicles):
            new_vehicle = Crossdock_In_Vehicle()
            new_vehicle._id = i+1
            new_vehicle._cd_id = i//2 + 1
            self.list_in_vehicles.append(new_vehicle)
            
        self.total_out_vehicles = self.num_cd*self.num_out_vehicle
        for i in range(self.total_out_vehicles):
            new_vehicle = Crossdock_Out_Vehicle()
            new_vehicle._id = i+1
            new_vehicle._cd_id = i//2 + 1
            self.list_out_vehicles.append(new_vehicle)
        
    ###############################################################################   
    def Module_Gene_Decode_1(self, gene_sequence):
        self.seq_1 = gene_sequence[0]   # retailer-supplier product assignment
        for i in range(len(self.seq_1)):
            supply_id = self.seq_1[i]
            self.list_products[i]._supplier = supply_id
            
    ###############################################################################
    def Module_Gene_Decode_2(self, gene_sequence):
        self.seq_2 = gene_sequence[1]   # assignment of products to cross-docks & inbound routing        
        # there will be (num_cd * num_vehicles - 1) zero indices in total
        
        # assigning cross-dock and vehicle number
        for i in range(len(self.seq_2)):
            prod_id = self.seq_2[i]
            if (prod_id == 0):
                continue
            else:
                list_short = self.seq_2[0:i+1]
                num_zeros = list_short.count(0)
                cd_id = (num_zeros//self.num_in_vehicle) + 1
                supplier_id = self.seq_1[prod_id-1]
                invehicle_id = (num_zeros%self.num_in_vehicle) + 1
                self.list_products[prod_id-1]._cd = cd_id
                self.list_products[prod_id-1]._invehicle = invehicle_id
                
                
                index = 2*(cd_id-1) + 1*(invehicle_id-1)
                self.list_in_vehicles[index]._route.append(supplier_id)
                self.list_in_vehicles[index]._prod_order.append(prod_id)
        
        # closing the inbound truck loop (0,2,2,0)
        print ("inbound loops")
        for i in range(len(self.list_in_vehicles)):
            #print ("truck:", i+1)
            self.list_in_vehicles[i]._route.append(0)
            self.list_in_vehicles[i]._prod_order.append(0)
            print (self.list_in_vehicles[i]._route, self.list_in_vehicles[i]._prod_order)
        
        # assigning the rank for pickup
        temp_rank = 0  
        for i in range(len(self.seq_2)):
            prod_id = self.seq_2[i]
            
            if (prod_id == 0):
                temp_rank = 0
            else:
                temp_rank += 1
                self.list_products[prod_id-1]._inrank = temp_rank
                    
    ###############################################################################
    def Module_Gene_Decode_3(self, gene_sequence):
        self.seq_3 = gene_sequence[2]   # consolidation decisions & outbound routing
        cd_items_list = [[] for i in range(self.num_cd)]
        
        for i in range(len(self.seq_3)):
            prod_id = self.seq_3[i]
            if (prod_id <= 0): 
                continue
            else:
                cd_id = self.list_products[prod_id-1]._cd
                cd_items_list[cd_id-1].append(prod_id)
        
        for i in range(len(cd_items_list)):
            sub_string = cd_items_list[i]
            for j in range(len(sub_string)):
                prod_id = sub_string[j]
                index_in_seq = self.seq_3.index(prod_id)
                cd_id = self.list_products[prod_id-1]._cd
                index_divide = self.seq_3.index(-cd_id)
                
                if (index_in_seq < index_divide):
                    self.list_products[prod_id-1]._outvehicle = 1
                elif (index_in_seq >= index_divide):
                    self.list_products[prod_id-1]._outvehicle = 2
        
        out_rank_list = [[1,1] for i in range(self.num_cd)]           
        flat_list = [item for sublist in cd_items_list for item in sublist]
        
        for i in range(len(flat_list)):
            prod_id = flat_list[i]
            cd_id = self.list_products[prod_id-1]._cd
            outvehicle_id = self.list_products[prod_id-1]._outvehicle
            out_rank = out_rank_list[cd_id-1][outvehicle_id-1]
            self.list_products[prod_id-1]._outrank = out_rank
            out_rank_list[cd_id-1][outvehicle_id-1] += 1
            
            retailer_id = (prod_id-1)//self.num_products + 1
            index = 2*(cd_id-1) + 1*(outvehicle_id-1)
            self.list_out_vehicles[index]._route.append(retailer_id)
            self.list_out_vehicles[index]._prod_order.append(prod_id) 
        
        # closing the outbound truck loop (0,2,2,0)
        print ("outbound loops")
        for i in range(len(self.list_out_vehicles)):
            #print ("truck:", i+1)
            self.list_out_vehicles[i]._route.append(0)
            self.list_out_vehicles[i]._prod_order.append(0)
            print (self.list_out_vehicles[i]._route, self.list_out_vehicles[i]._prod_order)
            
    ###############################################################################
    def Module_Gene_Decode_4(self, gene_sequence):
        self.seq_4 = gene_sequence[3]   # departure times of inbound vehicles from cross-docks
        for i in range(len(self.seq_4)):
            wait_time = self.seq_4[i]
            self.list_in_vehicles[i]._wait_time = wait_time
    
    ###############################################################################
    def Module_Print_Essentials(self):
        print ("id, retail, supply, cd, invehicle, inrank, outvehicle, outrank")
        for i in range(len(self.list_products)):
            print (self.list_products[i]._id, self.list_products[i]._retailer, 
                   self.list_products[i]._supplier, self.list_products[i]._cd,
                   self.list_products[i]._invehicle, self.list_products[i]._inrank,
                   self.list_products[i]._outvehicle, self.list_products[i]._outrank)
            
        
    ###############################################################################
    def Module_In_Fixed_Cost(self):
        list_cd_vehicle = []
        for i in range(len(self.list_products)):
            cd_id = self.list_products[i]._cd
            vehicle_id = self.list_products[i]._invehicle
            list_cd_vehicle.append((cd_id, vehicle_id))
            
        num_vehicles = len(set(list_cd_vehicle))
        in_fixed_cost = num_vehicles * self.f_ik_1
        #print (num_vehicles, in_fixed_cost) 
        return in_fixed_cost
    
    ###############################################################################
    def Module_In_Trans_Cost(self):
        total_dist = 0
        for k in range(len(self.list_in_vehicles)):
            route_dist = 0
            cd_id = self.list_in_vehicles[k]._cd_id
            route = self.list_in_vehicles[k]._route
            simple_route = list(dict.fromkeys(route))
            
            #print ("cd_id = ", cd_id)
            simple_route.append(0)
            print (simple_route)
            for i in range(len(simple_route)-1):
                current_pt = simple_route[i]
                next_pt = simple_route[i+1]
                if (current_pt == 0):
                    dist = self.cd_supplier_dist[cd_id-1][next_pt-1]
                elif (next_pt == 0):
                    dist = self.cd_supplier_dist[cd_id-1][current_pt-1]
                else:
                    dist = self.supplier_dist[current_pt-1][next_pt-1]
                #print ("dist = ", dist)
                route_dist += dist
                total_dist += dist
            #print (route_dist)
        print ("inbound distance = ", total_dist)
        total_cost = total_dist*self.cik_1
        print ("inbound cost = ", total_cost)
        return total_cost
    
    ###############################################################################
    def Module_Out_Fixed_Cost(self):
        list_cd_vehicle = []
        for i in range(len(self.list_products)):
            cd_id = self.list_products[i]._cd
            vehicle_id = self.list_products[i]._outvehicle
            list_cd_vehicle.append((cd_id, vehicle_id))
            
        num_vehicles = len(set(list_cd_vehicle))
        out_fixed_cost = num_vehicles * self.f_ik
        print (num_vehicles, out_fixed_cost) 
        return out_fixed_cost
    
    ###############################################################################
    def Module_Out_Trans_Cost(self):
        total_dist = 0
        for k in range(len(self.list_out_vehicles)):
            route_dist = 0
            cd_id = self.list_out_vehicles[k]._cd_id
            route = self.list_out_vehicles[k]._route
            simple_route = list(dict.fromkeys(route))
            
            #print ("cd_id = ", cd_id)
            simple_route.append(0)
            print (simple_route)
            for i in range(len(simple_route)-1):
                current_pt = simple_route[i]
                next_pt = simple_route[i+1]
                if (current_pt == 0):
                    dist = self.cd_retailer_dist[cd_id-1][next_pt-1]
                elif (next_pt == 0):
                    dist = self.cd_retailer_dist[cd_id-1][current_pt-1]
                else:
                    dist = self.retailer_dist[current_pt-1][next_pt-1]
                #print ("dist = ", dist)
                route_dist += dist
                total_dist += dist
            print (route_dist)
        print ("outbound distance = ", total_dist)
        total_cost = total_dist*self.cik
        print ("outbound cost = ", total_cost)
        return total_cost
    
    ###############################################################################
    def Module_Inbound_Cost(self):
        fixed_cost = self.Module_In_Fixed_Cost()
        trans_cost = self.Module_In_Trans_Cost()
        inbound_cost = fixed_cost + trans_cost
        return inbound_cost

    ###############################################################################
    def Module_Outbound_Cost(self):
        pass
        fixed_cost = self.Module_Out_Fixed_Cost()
        trans_cost = self.Module_Out_Trans_Cost()
        outbound_cost = fixed_cost + trans_cost
        return outbound_cost
    
    def Module_Storage_Cost(self):
        
        return 0
    ###############################################################################
    def Module_Calculate_Cost(self, gene_sequence):
        # Initialize all the parameters
        self.Module_Initialization()
        self.Module_Create_Dataframes()
        self.Module_Create_Products()
        self.Module_Create_Vehicles()
        self.Module_Gene_Decode_1(gene_sequence)
        self.Module_Gene_Decode_2(gene_sequence)
        self.Module_Gene_Decode_3(gene_sequence)
        self.Module_Gene_Decode_4(gene_sequence)
        #self.Module_Print_Essentials()
        self.inbound_cost = self.Module_Inbound_Cost()
        self.outbound_cost = self.Module_Outbound_Cost()
        self.storage_cost = self.Module_Storage_Cost()
        
        self.total_cost = self.inbound_cost + self.outbound_cost + self.storage_cost
        print (self.total_cost)
###############################################################################
# ========================================================
# ============ START THE EVOLUTIONARY PROCESS ============
# ========================================================

gene_seq = [[1,2,2,1,1,2,1,1,2,1],[2,3,0,4,1,0,5,9,8,0,6,7,10],[1,6,3,-1,8,10,4,2,-2,5,7,9],[50,220,0,90]]
sim = Simulation()
sim.Module_Calculate_Cost(gene_seq)

