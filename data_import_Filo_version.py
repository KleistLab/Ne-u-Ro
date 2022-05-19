#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 12:15:00 2022

@author: carolinabarata
"""

import numpy as np
import pandas as pd
import os
import re
import matplotlib.pyplot as plt
from os import path


plt.style.use('seaborn-darkgrid')


# Filopodia data import, based on the code and data found in the repository:
# https://github.com/vkleist/Filo

# The data structure created has the format:
    # Data_Structure[mutant][time][type]
    # mutant: 'DLar' or 'LiprinA' or 'Syd1' or 'Trio' or 'WT'
    # time: 'P40' or 'P60'
    # type: 'sF' or 'ellF'  (short- and long-lived filopodia)
    # Example: Data_Structure[DLar][P40][sF] contains all data of short-lived filopodia at
    # P40 for the DLar knockout mutant. The data here (and for all examples)
    # includes information about the start times, end times, lifetimes of filopodia,
    # and information regarding which growth centre does each filopdia belong to:
    # fields 'Start Time', 'End Time', 'Lifetime' and GC', respectively.    
   


Data_Structure = { 'WT': 
        { 'P40' :
         { 'sF' : { }, 'ellF': {}
         },
          'P60' :
         { 'sF' : { }, 'ellF': {}
         } 
        },
         'DLar':  
        { 'P40' :
         { 'sF' : { }, 'ellF': {}
         },
          'P60' :
         { 'sF' : { }, 'ellF': {}
         } 
        },
         'LiprinA':  
        { 'P40' :
         { 'sF' : { }, 'ellF': {}
         },
          'P60' :
         { 'sF' : { }, 'ellF': {}
         } 
        },
         'Syd1':  
        { 'P40' :
         { 'sF' : { }, 'ellF': {}
         },
          'P60' :
         { 'sF' : { }, 'ellF': {}
         } 
        },
         'Trio':  
        { 'P40' :
         { 'sF' : { }, 'ellF': {}
         },
          'P60' :
         { 'sF' : { }, 'ellF': {}
         } 
        },    
    }
    
    
mutants = ['DLar', 'LiprinA', 'Syd1', 'Trio', 'WT']  
    
#  import file assuming it is in csv format

# import file assuming it is in xlsx format

basepath = path.dirname(__file__)

path =  basepath + '/Filo_Data/FiloData/'

# this threshold enables to classify filopodia as unstable (lifetime < classification_threshold) or
# stable (lifetime >= classification_threshold)
classification_threshold = 8

for mutant in mutants:

    path_mutant = path + mutant
    
    files_in_path = os.listdir(path_mutant)
    
    if (mutant == 'WT'):
        growth_cone_id = 'gc'
        time_id = 'P'
        files_to_consider = [x for x in files_in_path if re.search('csv', x)]
        
    else:
        growth_cone_id = 'GC'
        time_id = 'fast'
        files_to_consider = [x for x in files_in_path if re.search('xlsx', x)]
        
    
    data_P40_sF = pd.DataFrame(columns=['Life Time', 'Start Time', 'End Time', 'GC'])
    data_P40_ellF = pd.DataFrame(columns=['Life Time', 'Start Time', 'End Time', 'GC'])
    data_P60_sF = pd.DataFrame(columns=['Life Time', 'Start Time', 'End Time', 'GC'])
    data_P60_ellF = pd.DataFrame(columns=['Life Time', 'Start Time', 'End Time', 'GC'])
    
    nr_filos_per_minute_sF_40_all_GC = []
    nr_filos_per_minute_ellF_40_all_GC = []
    nr_filos_per_minute_sF_60_all_GC = []
    nr_filos_per_minute_ellF_60_all_GC = []
    
    fig, axs = plt.subplots(2, 3, figsize=(12,6), dpi=400, sharex = 'all', sharey = 'col')
    time_plot = np.arange(0,60)
    
    fig_hist, axs_hist = plt.subplots(2, 2, figsize=(12,6), dpi=400, sharex = 'all', sharey = 'all')

    for i in range(len(files_to_consider)):
        
        name_file = files_to_consider[i]
        
        pos_GC =  name_file.find(growth_cone_id)
        pos_time =  name_file.find(time_id)
        
        GC_number = name_file[pos_GC+2:pos_GC+3]
        
        if (time_id == 'P'):
            time = name_file[pos_time:pos_time+3]
    
        else:
            time_name = name_file[pos_time+4:pos_time+5]
            if (time_name == '1'):
                time = 'P40'
            else: 
                time = 'P60'
                
        
        # for xlsx files 
        if "xlsx" in name_file:
            
            data = pd.read_excel(path_mutant + '/' + name_file, usecols = ['Life Time', 'Start Time', 'End Time'])
        
        # for csv files
        else:
            
            data = pd.read_csv(path_mutant + '/' + name_file, usecols = ['Life Time [min]', 'Start Time Step', 'End Time Step'])
            # rename columns
            data = data.rename(columns={'Life Time [min]': 'Life Time', 'Start Time Step': 'Start Time', 'End Time Step': 'End Time'})
        
        # make list as big as data frame with the same value for GC
        growth_cone = GC_number
        
        data['GC'] = growth_cone
        
        data_sF = data.loc[data['Life Time'] < classification_threshold]
        data_ellF = data.loc[data['Life Time'] >= classification_threshold]
        
        nr_filos_sF = len(data_sF)
        nr_filos_ellF = len(data_ellF)
        nr_filos_per_minute_sF = np.zeros ((60,1))
        nr_filos_per_minute_ellF = np.zeros ((60,1))
        # sF
        for counter in range(0,nr_filos_sF-1):
            start = data_sF.loc[data_sF.index[counter], 'Start Time']
            end = data_sF.loc[data_sF.index[counter], 'End Time']
            nr_filos_per_minute_sF[start:end + 1] = nr_filos_per_minute_sF[start:end + 1] + 1
        # ellF
        for counter in range(0,nr_filos_ellF-1):
            start = data_ellF.loc[data_ellF.index[counter], 'Start Time']
            end = data_ellF.loc[data_ellF.index[counter], 'End Time']
            nr_filos_per_minute_ellF[start:end + 1] = nr_filos_per_minute_ellF[start:end + 1] + 1
            
        
        if (time == 'P40'):
            data_P40_sF = pd.concat([data_P40_sF,data_sF])
            data_P40_ellF = pd.concat([data_P40_ellF, data_ellF])
            nr_filos_per_minute_sF_40_all_GC.append(np.transpose(nr_filos_per_minute_sF))
            nr_filos_per_minute_ellF_40_all_GC.append(np.transpose(nr_filos_per_minute_ellF))
            axs[0, 0].plot(time_plot, nr_filos_per_minute_sF)
            axs[0, 1].plot(time_plot, nr_filos_per_minute_ellF)

        else:
            data_P60_sF = pd.concat([data_P60_sF,data_sF])
            data_P60_ellF = pd.concat([data_P60_ellF, data_ellF])
            nr_filos_per_minute_sF_60_all_GC.append(np.transpose(nr_filos_per_minute_sF))
            nr_filos_per_minute_ellF_60_all_GC.append(np.transpose(nr_filos_per_minute_ellF))
            axs[1, 0].plot(time_plot, nr_filos_per_minute_sF)
            axs[1, 1].plot(time_plot, nr_filos_per_minute_ellF)

    
    # Add some statistics to the plots of the filopodia numbers per minute    
    
    means_sF_P40 = np.mean(nr_filos_per_minute_sF_40_all_GC, axis = 0)
    means_ellF_P40 = np.mean(nr_filos_per_minute_ellF_40_all_GC, axis = 0)
    means_sF_P60 = np.mean(nr_filos_per_minute_sF_60_all_GC, axis = 0)
    means_ellF_P60 = np.mean(nr_filos_per_minute_ellF_60_all_GC, axis = 0)
    
    axs[0, 0].plot(time_plot, np.transpose(means_sF_P40), color = 'k', linestyle='--')
    axs[0, 1].plot(time_plot, np.transpose(means_ellF_P40), color = 'k', linestyle='--')
    axs[1, 0].plot(time_plot, np.transpose(means_sF_P60), color = 'k', linestyle='--')
    axs[1, 1].plot(time_plot, np.transpose(means_ellF_P60), color = 'k', linestyle='--')
    
    ratio_P40 = means_ellF_P40/means_sF_P40
    ratio_P60 = means_ellF_P60/means_sF_P60
    
    axs[0, 2].plot(time_plot, np.transpose(ratio_P40), color = 'k')
    axs[0, 2].axhline(y=np.mean(ratio_P40), color='r', linestyle='--')
    axs[1, 2].plot(time_plot, np.transpose(ratio_P60), color = 'k')
    axs[1, 2].axhline(y=np.mean(ratio_P60), color='r', linestyle='--')
    
    axs[0, 0].set(ylabel='P40')
    axs[1, 0].set(xlabel = 'time (min)', ylabel='P60')
    axs[1, 1].set(xlabel = 'time (min)')
    axs[1, 2].set(xlabel = 'time (min)')
    
    axs[0, 0].set_title('Nr unstable filopodia (' + str(mutant) +')')
    axs[0, 1].set_title('Nr stable filopodia (' + str(mutant) +')')
    axs[0, 2].set_title('Ratio: nr stable to unstable filopodia (' + str(mutant) +')')
    
    fig.show()     
    
    
    # Plot histograms of filopodia numbers
    
    bin_edges = np.arange(-0.5, 20.5, 1)
    
    nr_filos_per_minute_sF_40_all_GC_for_histo = np.concatenate(nr_filos_per_minute_sF_40_all_GC).ravel()
    nr_filos_per_minute_ellF_40_all_GC_for_histo = np.concatenate(nr_filos_per_minute_ellF_40_all_GC).ravel()
    nr_filos_per_minute_sF_60_all_GC_for_histo = np.concatenate(nr_filos_per_minute_sF_60_all_GC).ravel()
    nr_filos_per_minute_ellF_60_all_GC_for_histo = np.concatenate(nr_filos_per_minute_ellF_60_all_GC).ravel()
    
    hist_val_sF_40, bin_edges = np.histogram(nr_filos_per_minute_sF_40_all_GC_for_histo, bins = bin_edges)
    hist_val_ellF_40, bin_edges = np.histogram(nr_filos_per_minute_ellF_40_all_GC_for_histo, bins = bin_edges)
    hist_val_sF_60, bin_edges = np.histogram(nr_filos_per_minute_sF_60_all_GC_for_histo, bins = bin_edges)
    hist_val_ellF_60, bin_edges = np.histogram(nr_filos_per_minute_ellF_60_all_GC_for_histo, bins = bin_edges)

    bar_edges = np.arange(-0.5, 19.5, 1)

    axs_hist[0,0].bar(bar_edges, hist_val_sF_40/np.sum(hist_val_sF_40), align = 'edge', width = 1)
    axs_hist[0,1].bar(bar_edges, hist_val_ellF_40/np.sum(hist_val_ellF_40), align = 'edge', width = 1)
    axs_hist[1,0].bar(bar_edges, hist_val_sF_60/np.sum(hist_val_sF_60), align = 'edge', width = 1)
    axs_hist[1,1].bar(bar_edges, hist_val_ellF_60/np.sum(hist_val_ellF_60), align = 'edge', width = 1)
   
    mean_sF_40 = np.mean(nr_filos_per_minute_sF_40_all_GC_for_histo)
    std_sF_40 = np.std(nr_filos_per_minute_sF_40_all_GC_for_histo)
    mean_ellF_40 = np.mean(nr_filos_per_minute_ellF_40_all_GC_for_histo)
    std_ellF_40 = np.std(nr_filos_per_minute_ellF_40_all_GC_for_histo)
    mean_sF_60 = np.mean(nr_filos_per_minute_sF_60_all_GC_for_histo)
    std_sF_60 = np.std(nr_filos_per_minute_sF_60_all_GC_for_histo)
    mean_ellF_60 = np.mean(nr_filos_per_minute_ellF_60_all_GC_for_histo)
    std_ellF_60 = np.std(nr_filos_per_minute_ellF_60_all_GC_for_histo)
    
    axs_hist[0, 0].axvline(x=mean_sF_40, color='r', linestyle='--')
    axs_hist[0, 0].axvline(x=mean_sF_40 - std_sF_40, color='r', linestyle=':')
    axs_hist[0, 0].axvline(x=mean_sF_40 + std_sF_40, color='r', linestyle=':')
    axs_hist[0, 1].axvline(x=mean_ellF_40, color='r', linestyle='--')
    axs_hist[0, 1].axvline(x=mean_ellF_40 - std_ellF_40, color='r', linestyle=':')
    axs_hist[0, 1].axvline(x=mean_ellF_40 + std_ellF_40, color='r', linestyle=':')
    axs_hist[1, 0].axvline(x=mean_sF_60, color='r', linestyle='--')
    axs_hist[1, 0].axvline(x=mean_sF_60 - std_sF_60, color='r', linestyle=':')
    axs_hist[1, 0].axvline(x=mean_sF_60 + std_sF_60, color='r', linestyle=':')
    axs_hist[1, 1].axvline(x=mean_ellF_60, color='r', linestyle='--')
    axs_hist[1, 1].axvline(x=mean_ellF_60 - std_ellF_60, color='r', linestyle=':')
    axs_hist[1, 1].axvline(x=mean_ellF_60 + std_ellF_60, color='r', linestyle=':')

    axs_hist[0, 0].set(ylabel='P40')    
    axs_hist[1, 0].set(xlabel = 'time (min)', ylabel='P60')
    axs_hist[1, 1].set(xlabel = 'time (min)')
    axs_hist[0, 0].set_title('Nr unstable filopodia (' + str(mutant) +')')
    axs_hist[0, 1].set_title('Nr stable filopodia (' + str(mutant) +')')

    axs_hist[0,0].tick_params(axis='both', labelleft = 'True', labelbottom = 'True')
    axs_hist[0,1].tick_params(axis='both', labelleft = 'True', labelbottom = 'True')
    axs_hist[1,0].tick_params(axis='both', labelleft = 'True', labelbottom = 'True')
    axs_hist[1,1].tick_params(axis='both', labelleft = 'True', labelbottom = 'True')

    
    axs_hist[0,0].set_xlim([-1, 20.5])
    axs_hist[0,1].set_xlim([-1, 20.5])
    axs_hist[1,0].set_xlim([-1, 20.5])
    axs_hist[1,1].set_xlim([-1, 20.5])
    
    axs_hist[0,0].set_ylim([0, 0.3])
    axs_hist[0,1].set_ylim([0, 0.3])
    axs_hist[1,0].set_ylim([0, 0.3])
    axs_hist[1,1].set_ylim([0, 0.3])
    
    fig_hist.show()  
    
    
    # Save the data fields in the data structure
    
    Data_Structure[mutant]['P40']['sF'] = { 'Life Time' : data_P40_sF['Life Time'], 
                                           'Start Time': data_P40_sF['Start Time'],
                                           'End Time': data_P40_sF['End Time'],
                                           'GC': data_P40_sF['GC']}
    
    Data_Structure[mutant]['P40']['ellF'] = { 'Life Time' : data_P40_ellF['Life Time'], 
                                           'Start Time': data_P40_ellF['Start Time'],
                                           'End Time': data_P40_ellF['End Time'],
                                           'GC': data_P40_ellF['GC']}
    
    Data_Structure[mutant]['P60']['sF'] = { 'Life Time' : data_P60_sF['Life Time'], 
                                           'Start Time': data_P60_sF['Start Time'],
                                           'End Time': data_P60_sF['End Time'],
                                           'GC': data_P60_sF['GC']}
    
    Data_Structure[mutant]['P60']['ellF'] = { 'Life Time' : data_P60_ellF['Life Time'], 
                                           'Start Time': data_P60_ellF['Start Time'],
                                           'End Time': data_P60_ellF['End Time'],
                                           'GC': data_P60_ellF['GC']}

        
        




