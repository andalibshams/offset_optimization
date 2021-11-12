### Link-pivot algorithm

def link_pivot_algorithm(offsetAdj, method_to_use):
        delta = np.zeros(number_of_intersection, dtype = np.int16)

    ### check all the indexes

        for i in range(number_of_intersection): 
                for j in range(i):
                        offsetAdj[j] = 0
                        for k in range(i-1, j-1, -1):
                                offsetAdj[j] += delta[k]

                        offsetAdj[j] = offsetAdj[j]%cycle_length


        bestMOE = performance(offsetAdj, method_to_use)
        bestLPI = 0

        for j in range(cycle_length):
                offsetAdj[i] = j

                for k in range(i):
                        if (j>0):
                                offsetAdj[k] = (offsetAdj[k]+1)%cycle_length
                sys_performance = performance(offsetAdj, method_to_use)

                if(sys_performance<bestMOE): 
                        bestLPI = j
                        delta[i] = j
                        bestMOE = sys_performance


        for i in range(number_of_intersection):
                offsetAdj[i] = 0
                for k in range(number_of_intersection-1, i-1, -1): 
                        offsetAdj[i] += delta[k]

        offsetAdj[i] = offsetAdj[i]%cycle_length