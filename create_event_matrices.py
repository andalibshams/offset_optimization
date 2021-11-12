# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 18:04:00 2021

@author: ashams
"""

#%%
import pandas as  pd 
import numpy as np
from scipy import stats


#%%

def mat_set_global_variables(cyc_len, num_inter, sim_start, sim_end, input_df):
    global cycle_length, number_of_intersection, simulation_start, simulation_time
    cycle_length = cyc_len
    number_of_intersection = num_inter
    simulation_start = sim_start
    simulation_time = sim_end
    
    global number_of_cycle, cycle_remainder
    number_of_cycle = simulation_time//cycle_length
    cycle_remainder = simulation_time%cycle_length
    
    global input_parameters
    input_parameters = input_df

def mat_get_global_variable():
    print("Cycle Length: ", cycle_length)
    print("Total Number of Intersection: ", number_of_intersection)
    print("Simulation Start: ", simulation_start)
    print("Simulation End Time: ", simulation_time)
    
#%%

def create_probability_of_green(signal, shifter = 0):
    green_time = np.zeros((number_of_intersection,2, (number_of_cycle+1)*cycle_length))

    for link_ in range(len(input_parameters)):
        intersection_num = input_parameters["Intersection Number"][link_]
        intersection_ = input_parameters["intersection_id"][link_]

        sg_ = input_parameters.sg_id[link_]

        signal_ = signal.loc[(signal.SC==intersection_num) & (signal.SG==sg_)]

        if(sg_==2): bound = 0
        else: bound  = 1

        for step in range(0, len(signal_), 2):
            bog = signal_.iloc[step, 5]

            if(step<len(signal_) - 1): eog = signal_.iloc[step+1, 4]
            else: eog = simulation_time

            assert(bog+shifter<eog) ### checks if the shifter value is not too high. 
            
            green_time[intersection_ -1,bound, bog+shifter :eog] = 1
        
    green_distribution = green_time.reshape(number_of_intersection,2, (number_of_cycle+1), cycle_length)
    
    cycle_sum = np.ones((cycle_length, 1)) * number_of_cycle
    cycle_sum[:cycle_remainder, 0] += 1

    return np.sum(green_distribution, axis = 2)/cycle_sum.T, green_time


#%%

def create_eog_split_matrix(signal):
    eog_matrix = np.zeros((number_of_intersection, 2, number_of_cycle+10))
    split_matrix = np.zeros((number_of_intersection, 2, number_of_cycle+10))

    eog_cycle_num = np.zeros((number_of_intersection, 2), dtype = np.int16)
    split_cycle_num = np.zeros((number_of_intersection, 2), dtype = np.int16)
    
    for link_ in range(len(input_parameters)):
        intersection_num = input_parameters["Intersection Number"][link_]
        intersection_ = input_parameters["intersection_id"][link_]

        sg_ = input_parameters.sg_id[link_]

        signal_ = signal.loc[(signal.SC==intersection_num) & (signal.SG==sg_)]

        if(sg_==2): bound = 0
        else: bound  = 1

        for step in range(0, len(signal_), 2):
            bog = signal_.iloc[step, 5]

            if(step<len(signal_) - 1):
                eog = signal_.iloc[step+1, 4]

                if(eog>0.1): 
                    eog_matrix[intersection_ - 1][bound][eog_cycle_num[intersection_ - 1][bound]] = eog%cycle_length
                    eog_cycle_num[intersection_ - 1][bound] += 1

            else: eog = simulation_time

            if(bog>0.1 and eog!=simulation_time):
                split_matrix[intersection_ - 1][bound][split_cycle_num[intersection_ - 1][bound]] = eog - bog
                split_cycle_num[intersection_ - 1][bound] += 1    
    
    return eog_matrix, split_matrix, eog_cycle_num, split_cycle_num


#%%

def create_green_matrix_for_bw(split_matrix, split_cycle_num, eog_matrix, eog_cycle_num):
    expcted_split = np.round(np.sum(split_matrix, axis = 2)/ split_cycle_num, 0)
    expcted_split = np.minimum(expcted_split, cycle_length)

    final_green_time = np.zeros((number_of_intersection,2, cycle_length))

    for i in range(number_of_intersection):
        for j in range(2):
            # split_current = split_matrix[i, j, :split_cycle_num[i][j]]
            # expcted_split2 = int(np.round(stats.mode(split_current)[0], 0))
    
            eog = int(stats.mode(eog_matrix[i,j, :eog_cycle_num[i][j]])[0])
            bog = (eog - int(expcted_split[i][j]))%cycle_length
    
            if(eog>bog): 
                final_green_time[i, j, bog:eog] = 1
            else: 
                final_green_time[i, j, bog:cycle_length] = 1
                final_green_time[i, j, 0:eog] = 1

    #check if everything is right
    assert(np.sum(np.sum(final_green_time, axis = 2) - expcted_split)==0)
                           
    return final_green_time

#%%

def align_green_with_travel_time(green_matrix):
    link_ = 0
    for i in range(number_of_intersection):
        for j in range(2):
            tt = int(input_parameters["travel time"][link_])

            placeHolder_green  = np.hstack((green_matrix[i, j, cycle_length - tt: cycle_length], \
                                            green_matrix[i, j, 0: cycle_length - tt]))

            green_matrix[i, j, :] = placeHolder_green[:]

            link_ += 1

    return np.array(green_matrix)

#%%
def create_capactiy_matrix(green_aog):
    capacity = np.zeros((number_of_intersection,2, cycle_length))

    for link_ in range(len(input_parameters)):
        intersection_ = input_parameters["intersection_id"][link_]
        sg_ = input_parameters.sg_id[link_]
        vpb = input_parameters.vpb[link_]

        if(sg_==2): bound = 0
        else: bound  = 1
        capacity[intersection_ -1,bound, :] = green_aog[intersection_ -1,bound, :] * vpb * number_of_cycle

    return  np.array(capacity)

#%%

def create_vehicle_distribution_matrix(travel_time):
    vehicle_arrival = np.zeros((number_of_intersection,2, (number_of_cycle+1)*cycle_length))

    for link_ in range(len(input_parameters)):

        if(input_parameters.loc[link_, "sg_id"]==2): bound = 0
        else: bound  = 1

        intersection_ = input_parameters.loc[link_, "intersection_id"] -1

        if(input_parameters.loc[link_, "LinkFromInt"]==-1): 
            vehicle_arrival[intersection_,bound, :] = 0
            continue
                
        travel_time_intersection = input_parameters.loc[link_, "travel_time"]
        travel_time_id = input_parameters.loc[link_, "tt_id"]

        veh_tt = travel_time.loc[travel_time.tt_id == travel_time_id].copy()    
        veh_tt.loc[:, "timestamp"] = np.round((veh_tt.time - veh_tt.tt + travel_time_intersection), 0)
        veh_tt.timestamp = veh_tt["timestamp"].astype(np.int64)

        unique, counts = np.unique(veh_tt.timestamp, return_counts=True)
        vehicle_arrival[intersection_,bound, unique] = counts

    vehicle_arrival_reshape = vehicle_arrival.reshape(number_of_intersection,2, (number_of_cycle+1), cycle_length) 
    
    return np.sum(vehicle_arrival_reshape, axis = 2)

#%%

def create_CV_distribution_matrix(travel_time):
    vehicle_arrival = np.zeros((number_of_intersection,2, (number_of_cycle+1)*cycle_length))

    for link_ in range(len(input_parameters)):

        if(input_parameters.loc[link_, "sg_id"]==2): bound = 0
        else: bound  = 1

        intersection_ = input_parameters.loc[link_, "intersection_id"] -1

        if(input_parameters.loc[link_, "LinkFromInt"]==-1): 
            vehicle_arrival[intersection_,bound, :] = 0
            continue
                
        travel_time_intersection = input_parameters.loc[link_, "travel_time"]
        travel_time_id = input_parameters.loc[link_, "tt_id"]
        upstream_tt_id = input_parameters.loc[link_, "Prev_tt_id"]
        
        vehicles_upstream_intersection =  np.array(travel_time.loc[travel_time.tt_id == upstream_tt_id, "vehNo"])
        veh_tt = travel_time.loc[(travel_time.tt_id == travel_time_id) & (travel_time['vehNo'].isin(vehicles_upstream_intersection))].copy()    
        
        veh_tt.loc[:, "timestamp"] = np.round((veh_tt.time - veh_tt.tt + travel_time_intersection), 0)
        veh_tt.timestamp = veh_tt["timestamp"].astype(np.int64)

        unique, counts = np.unique(veh_tt.timestamp, return_counts=True)
        vehicle_arrival[intersection_,bound, unique] = counts

    vehicle_arrival_reshape = vehicle_arrival.reshape(number_of_intersection,2, (number_of_cycle+1), cycle_length) 
    
    return np.sum(vehicle_arrival_reshape, axis = 2)
