#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 16 23:52:36 2018

@author: dimitra
"""
from __future__ import division
import numpy as np
from scipy.spatial import distance
import scipy as sc









def solve_per_neuron(idn,N, spiketimes, taff,ISIs, num_interv, num_events):
    
    ksp = 20 # estimated maximum number of cross-events per inter-spike interval
    # in case of larger number of cross- events create Bag_of_times again
    rep = 0
    while rep < 10:
        try:
            #create empty bag of times
            Bag_of_times =  [np.zeros((N,ksp)) for mi in range(num_interv)] 
            
            """"
            Bag_of_times = [W_i_1, .......,   W_i_si]     ,where s1 is the number of ISI of neuron i
                                                          ,and W_1_1 is the array with the relative times of spikes of the other neurons
                                                          that occured during the 1st ISI of neuron i
                                                          W is an array of size N x ksp , where ksp is the maximum number of times a neuron
                                                          fired during the ISI of another neuron
            """
            #actual max number of cross interval events 
            num_of_IIE = 0 
            k_per_interval = []            
        
            start = 0 # starting point of interval in spiketimes array
            stop = 1 # end point
            pivot =  [0 for _ in range(N) ] # points where I have stoped reading for every neuron
            while stop <= num_interv:
                
                for j in range(N):
                    if ( j != idn ):
                        taf = taff[j]
                        count_j_spikes = 0 # antistoixei sto ksp tou pinaka
                        while start == 0 and (spiketimes[idn][start] > spiketimes[j][pivot[j]]+taf) :      #=
                            pivot[j] = pivot[j] +1
                        #if ( len(spiketimes[j]) > pivot[j] ) and spiketimes[idn][start] == spiketimes[j][pivot[j]] +taf :
                        #    pivot[j] += 1
                        while (( len(spiketimes[j]) > pivot[j] ) and ( spiketimes[idn][stop] > spiketimes[j][pivot[j]] +taf ) and ( spiketimes[idn][start] <= spiketimes[j][pivot[j]] +taf ) ):
                            Bag_of_times[start][j,count_j_spikes] = spiketimes[j][pivot[j]] - spiketimes[idn][start] 
                            count_j_spikes = count_j_spikes + 1
                            pivot[j] = pivot[j] +1
                        
                        if len(k_per_interval)-1 < start:
                            k_per_interval.append( count_j_spikes )
                        else:
                            if k_per_interval[start] < count_j_spikes:
                                k_per_interval[start] = count_j_spikes
                        
                        if num_of_IIE < count_j_spikes  :
                            num_of_IIE = count_j_spikes 
                start = start + 1
                stop = stop + 1   
                
        except IndexError:
           rep += 1
           ksp += 20
           print('{:#<3} {:^50} {:#>3}'.format('','Oooops!!! More than %d cross-intervals per interspike- interval! Retrying...'%ksp,''))
           continue  
        break

    Bag_of_times3 =  [ timesarray[:,:num_of_IIE] for timesarray in  Bag_of_times ]
  
    del Bag_of_times
    #### Flatten - create vectors   
    # the two list comprehensions with enummerate are used to acces the np.arrayws W and flatten them,
    # and append used to add the ISI length 
    Flat_times =[np.append(b.flatten('F'),ISIs[inx2]) for inx2,b in enumerate(Bag_of_times3)]   
    # has the same structure as Bag of times 1st index refers to the id of ISI of the neuron 
    # 2nd index to the relative time of spiking of jth neuron
    
    
    #### Delete columns that corespond to the current unit    
    #create indexes for deletion of columns that correspond to the current unit i
    inde = [ idn+j*N  for j in range(num_of_IIE) ]    
    inde = inde[::-1]      
    Flat_times = np.delete(Flat_times, inde, 1) # delete columns thst correspond to the current unit i

    # calculated distances among events
    D = distance.cdist(Flat_times, Flat_times, 'euclidean') 
    # find event with minimum distance to all other events
    medoid =  np.argmin(np.mean(D,axis= 1 ))   
    
    #Select N*num_of_IIE points for every  equation    
    indexes_for_solving = [] 
    sorted_points = np.argsort(D[medoid])
    #indexes_for_solving[i] = sorted_points[1:(N*num_of_IIE[i])+1+additional]    #leave out the 1st because this is the medoid
    indexes_for_solving = sorted_points[1:] 
    
    k_per_interval = np.array(k_per_interval)    
    

    indx = indexes_for_solving[:num_events]
    k = np.max(k_per_interval[indx])
    
    new_Flat_times = Flat_times[:,:k*(N-1)]
    Flat_times = np.concatenate((new_Flat_times , (Flat_times[:,-1].reshape(-1,1))),axis  =1)
    
    del new_Flat_times
    
    ### Linearization center
    Linearization_center = np.array( Flat_times[medoid] )  
    Linearization_center.shape = (1,Linearization_center.shape[0] )
    
    # Calculate differences between reference event and other events
    WDT = Flat_times[indx,:] - Linearization_center[:,:]
    
    ### Cross interval matrix
    W  = WDT[:,:-1]
    
    ### Inter spike interval vector
    DT = WDT[:,-1]

    try:
        b = np.dot(np.linalg.pinv(W),DT)
    except np.linalg.linalg.LinAlgError as err:
        b = np.dot(sc.linalg.pinv(W),DT)





    results = b.reshape((-1,N-1))
    #### Add the estimated gradiends across k spiking profiles 
    result = np.sum(results, axis = 0) 
    
    end_result = np.concatenate( [ result[:idn] , [0] , result[idn:] ]   )
    
    return end_result






















if __name__ == "__main__":
    
    print('{:#<58}'.format(''))
    print('{:#<3} {:^50} {:#>3}'.format('','Reconstructing connectivity from spiketimes',''))
    # import spike data (list of lists)
    # first dimension corresponds to neuron id, second dimension to spike id
    spiketimes = [[]]
    
    
    #number of observed neurons
    N = len(spiketimes)
    #convert spiketriains to interspike intervals
    ISIs = map(np.diff, spiketimes)
    #number of recorded intervals per neuron
    num_interv_per_neuron = map(len,ISIs)
    
    
    # estimated connectivity
    # first dimension pre-synaptic neuron
    # second dimension post-synaptic
    Jhat = np.zeros((N,N))
    
    # interaction delays same structure as connectivity matrix
    # imteraction delay may be identified by optimasation procedure
    De = np.zeros((0,0))
    # number of events that will be employed in the solution
    num_events = 500 
    for i in range(N):
        
        Jhat[:,i] = solve_per_neuron(i,N, spiketimes, De[:,i], ISIs[i], num_interv_per_neuron[i], num_events)
    