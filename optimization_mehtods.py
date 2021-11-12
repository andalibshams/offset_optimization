# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 14:54:42 2021

@author: ashams
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#%%

def set_global_variables(cyc_len, num_inter, sim_start, sim_end, input_df):
    global cycle_length, number_of_intersection, simulation_start, simulation_time
    cycle_length = cyc_len
    number_of_intersection = num_inter
    simulation_start = sim_start
    simulation_time = sim_end
    
    global input_parameters
    input_parameters = input_df

def get_global_variable():
    print("Cycle Length: ", cycle_length)
    print("Total Number of Intersection: ", number_of_intersection)
    print("Simulation Start: ", simulation_start)
    print("Simulation End Time: ", simulation_time)
    
#%%
def set_fixed_global_matrices(fin_prog, fin_veh_dist, fin_green_time, fin_band_prog, final_capacity):
    global final_probability_of_green, final_vehicle_distribution, final_green_time, green_PrOG, fin_capacity
    
    final_probability_of_green = fin_prog
    final_vehicle_distribution = fin_veh_dist 
    final_green_time = fin_green_time
    green_PrOG = fin_band_prog
    fin_capacity = final_capacity

def set_variable_global_matrices():
    global green_aog, vehicle_distribution, capacity, green, PrOG
    
    green_aog = np.array(final_probability_of_green)
    vehicle_distribution = np.array(final_vehicle_distribution)
    capacity = np.array(fin_capacity)
    green = np.array(final_green_time) 
    PrOG = np.array(green_PrOG)
    
def reset_matrices():
    green_aog[:, :, :] = np.array(final_probability_of_green)
    vehicle_distribution[:, :, :] = np.array(final_vehicle_distribution)
    capacity[:, :, :] = np.array(fin_capacity)
    green[:, :, :] = np.array(final_green_time) 
    PrOG[:, :, :] = np.array(green_PrOG)
    
def get_matrix(matrix_type): 
    if(matrix_type == 0): return green_aog
    if(matrix_type == 1): return vehicle_distribution
    if(matrix_type == 2): return capacity
    if(matrix_type == 3): return green
    if(matrix_type == 4): return PrOG

#%%

def make_adjustments(bin_GAdj, bin_VAdj, offsetAdj):
    link_ = 0

    for intersection in range(number_of_intersection):
        for j in range(2):
            if(input_parameters.loc[link_, "sg_id"]==2): bound = 0
            else: bound  = 1

            next_Intersection = input_parameters.loc[link_, "nextIntersection"]

            next_intersection_valid = 0
            if(next_Intersection>=0 and next_Intersection< number_of_intersection):
                    next_intersection_valid = 1

            shift = int(offsetAdj[intersection])


            if(next_intersection_valid):
                placeHolder_VA  = np.hstack(( bin_VAdj[next_Intersection, bound, cycle_length - shift: cycle_length], \
                                              bin_VAdj[next_Intersection, bound, 0: cycle_length - shift] ))

                bin_VAdj[next_Intersection, bound, :] = placeHolder_VA[:]


            placeHolder_GA  = np.hstack(( bin_GAdj[intersection, bound, cycle_length - shift: cycle_length], \
                                      bin_GAdj[intersection, bound, 0: cycle_length - shift]))

            bin_GAdj[intersection, bound, :] = placeHolder_GA[:]

            link_ += 1

#%%

def plot_distributions_matrices(probability_of_green, vehicle_distribution, fig_size):
    signal_seq = np.arange(0,cycle_length, 1)
    veh_seq = np.arange(0,cycle_length, 1)

    cols = (number_of_intersection)
    rows = 2

    fig, axes = plt.subplots(nrows = rows, ncols = cols,sharey=True, sharex = True, figsize=fig_size)

    fig.tight_layout()

    for link_ in range(len(input_parameters)):
        i = input_parameters["intersection_id"][link_]-1

        if(input_parameters.loc[link_, "sg_id"]==2): j = 1
        else: j  = 0

        ax = axes[j][i]
        ax.set_title("Intersection Number: " + str(i+1) + " SG: " +str(input_parameters.loc[link_, "sg_id"]))
        ax.grid(which = 'both', linestyle = 'dashed')
        ax.grid(b=True, which='minor', alpha=0.2)
        ax.minorticks_on()
        ax.set_axisbelow(True)

        try: ax.bar(signal_seq,probability_of_green[i][1-j], width = 1, color = 'green', alpha = 0.8)
        except: breakpoint()
        
        ax2=ax.twinx()
        ax2.bar(veh_seq,vehicle_distribution[i][1-j], width = 1, color = 'black', alpha = 0.6)


#%%
def green_adjustments(bin_GAdj, offsetAdj):
    for i in range(number_of_intersection):
        offset_value = int(offsetAdj[i])
        for j in range(2):
            bin_GAdj[i, j, :] = np.hstack(( bin_GAdj[i, j, cycle_length - offset_value: cycle_length], \
                                            bin_GAdj[i, j, 0: cycle_length - offset_value]))


#%% Arrival on green
def arrival_on_green(bin_GAdj, bin_VAdj): 

    bin_AOG = np.multiply(bin_GAdj, bin_VAdj) 
    sys_AOG = np.sum(bin_AOG)

    return sys_AOG

#%% badnwidth

def bandwidth(bin_GAdj):
    return np.sum(np.prod(bin_GAdj, axis = 0))

#%% estimate delay

def delay(capacity_matrix, bin_VAdj):
### If (overall) volume exceeds capacity this method will not work

    queue= np.zeros(bin_VAdj.shape)

    veh_left = bin_VAdj - capacity_matrix

    ### first iteration    
    for i in range(cycle_length):
            queue[:, :, i] = np.maximum(queue[:, :,i-1] + veh_left[:, :, i], 0)


### second iteration
    for i in range(cycle_length):
            queue[:, :, i] = np.maximum(queue[:, :,i-1] + veh_left[:, :, i], 0)

    delay = np.sum(queue)/3600

    return delay

#%% def 

def number_of_stops(capacity_matrix, bin_VAdj):
    ### If (overall) volume exceeds capacity this method will not work
    
        queue= np.zeros(bin_VAdj.shape)

        veh_left = bin_VAdj - capacity_matrix

        ### first iteration
        for i in range(cycle_length):
                queue[:, :, i] = np.maximum(queue[:, :, i-1] + veh_left[:, :, i], 0)
    
    ### second iteration
        total_stops = 0
    
        for i in range(cycle_length):
                queue[:, :, i]   = np.maximum(queue[:, :, i-1] + veh_left[:, :, i], 0)

                stops = np.maximum(bin_VAdj[:, :, i] - np.maximum(capacity_matrix[:, :, i] - queue[:, :, i-1], 0), 0)

                total_stops += np.sum(stops)

        return total_stops


#%%

def performance(offsetAdj, method_to_use):
    ## 1: arrival profile
    ## 2: delay model
    ## 3: Bandwidth maximization (mode_eog - mean_bog)
    ## 4: PrOG-maximize
    ## 5: Number of stops
    ## 6: delay + 20*number of stops

    if(method_to_use==1):
    # Reset matrix
        green_aog[:, :, :] = np.array(final_probability_of_green)
        vehicle_distribution[:, :, :] = np.array(final_vehicle_distribution)

        # Offset Adjustment
        make_adjustments(green_aog, vehicle_distribution, offsetAdj)

        return -arrival_on_green(green_aog, vehicle_distribution)

    if(method_to_use==2):
    # Reset matrix
        capacity[:, :, :] = np.array(fin_capacity)
        vehicle_distribution[:, :, :] = np.array(final_vehicle_distribution)

        # Offset Adjustment
        make_adjustments(capacity, vehicle_distribution, offsetAdj)

        return delay(capacity, vehicle_distribution)

    if(method_to_use==3):
        green[:, :, :] = np.array(final_green_time) 
        green_adjustments(green, offsetAdj)

        return  -bandwidth(green)

    if(method_to_use==4):

        PrOG[:, :, :] = np.array(green_PrOG)
        green_adjustments(PrOG, offsetAdj)

        return  -bandwidth(PrOG)

    if(method_to_use==5):
        # Reset matrix
        capacity[:, :, :] = np.array(fin_capacity)
        vehicle_distribution[:, :, :] = np.array(final_vehicle_distribution)

        # Offset Adjustment
        make_adjustments(capacity, vehicle_distribution, offsetAdj)

        return number_of_stops(capacity, vehicle_distribution)

    if(method_to_use==6):
        # Reset matrix
        capacity[:, :, :] = np.array(fin_capacity)
        vehicle_distribution[:, :, :] = np.array(final_vehicle_distribution)

        # Offset Adjustment
        make_adjustments(capacity, vehicle_distribution, offsetAdj)

        return delay(capacity, vehicle_distribution) + 20/3600 * number_of_stops(capacity, vehicle_distribution)

#%%

### Hill-Climbing Algorithm

# to get faster speed stop the plotting

def mutate_solution(solution, neighbors): 
    index = np.random.randint(0,number_of_intersection)
    
    change = neighbors[np.random.randint(0,len(neighbors))]
    solution[index] =  np.mod(solution[index]+ change, cycle_length)


def hill_Climbing_algorithm(offsetAdj, method_to_use, neighbors,number_of_iteration, number_of_random_start, plotting):

    fig, axes = plt.subplots(1, 1)
    axes.set_xlabel("Iteration Number")
    axes.set_ylabel("Performance Index")

    axes.grid(which = 'both', linestyle = 'dashed')
    axes.grid(b=True, which='minor', alpha=0.2)
    axes.minorticks_on()
    #     axes.set_xlim([0, number_of_iteration])

    np.random.seed(1234)
    best = np.zeros((number_of_random_start, number_of_intersection), dtype = np.int32)
    best_score = np.zeros(number_of_random_start)

    ys = []
    xs = []

    i = 0
    j = 0
    rand_start = 0

    #     offsetAdj = np.zeros(number_of_intersection, dtype = np.int32)

    best[rand_start, :] = offsetAdj
    best_score[rand_start] = performance(offsetAdj, method_to_use)

    while(1):

        mutate_solution(offsetAdj, neighbors)

        score = performance(offsetAdj, method_to_use)

        if(plotting): 
            axes.scatter(i, score, color = 'black', s = 2)
            if(rand_start%2): axes.scatter(i, best_score[rand_start], color = 'red', s = 2)
            else: axes.scatter(i, best_score[rand_start], color = 'blue', s = 2)
            fig.canvas.draw()
        #         print(best)

        if score < best_score[rand_start]:
            best[rand_start, :] = offsetAdj
            best_score[rand_start] = score

            j = 0

            xs.append(i)
            ys.append(best_score[rand_start])


        offsetAdj[:] = best[rand_start, :] 
        #         print(i, " ", j, " ", offsetAdj[:])
        i += 1
        j += 1

        if(j>=number_of_iteration):
            rand_start += 1
            j = 0
            if(rand_start>=number_of_random_start): break

            offsetAdj[:] = np.random.randint(number_of_intersection)

            best[rand_start, :] = offsetAdj
            best_score[rand_start] = performance(offsetAdj, method_to_use)


    if(plotting==0): axes.plot(xs, ys)
    print(best_score)
    
    max_idx = np.argmin(best_score)
    
    offsetAdj[:] = best[max_idx, :]

#%%

def plot_distribution(parent, fig_size):
    ### if original = 1; print parent matrix otherwise modified one
    if(parent): plot_distributions_matrices(final_probability_of_green, final_vehicle_distribution, fig_size)
    else: plot_distributions_matrices(green_aog, vehicle_distribution, fig_size)
    
#%%

def offset_all_matrices(offsetAdj):
    make_adjustments(green_aog, vehicle_distribution, offsetAdj)
    
    placeHolder_veh_dist = np.array(vehicle_distribution)
    make_adjustments(capacity,placeHolder_veh_dist , offsetAdj)
    
    green_adjustments(green, offsetAdj)
    green_adjustments(PrOG, offsetAdj)

    
#%%
def print_output(offsetAdj): 
    reset_matrices()
    
    offset_all_matrices(offsetAdj)
    
    sys_aog = arrival_on_green(green_aog, vehicle_distribution)
    sys_delay = delay(capacity, vehicle_distribution)
    sys_Number_Stops = number_of_stops(capacity, vehicle_distribution)

    sys_bandwidth = bandwidth(green)
    sys_pog = bandwidth(PrOG)
    
    print("System AOG: %.2f" %sys_aog)
    print("System POG: %.2f" %sys_pog)
    print("System Delay: %.2f" %sys_delay)
    print("System Bandwidth: %.0f" %sys_bandwidth)
    print("System Number of Stops: %.2f" %sys_Number_Stops)
    