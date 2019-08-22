#!/usr/bin/env python
# coding: utf-8

import math, random, statistics
import xml.etree.ElementTree as xmlET
import numpy as np

np.random.seed(123)


# In[12]:


class SpotArray:
    # Used for bundling together data pertaining to a set of spots (donors or acceptors)
    # Normally contains four numpy arrays: 
    # .x, .y, ,.z for coordinates
    # .t_unav for containing the time until which the spot will be unavailable
    pass  

def build_spots_inside_ball(Nspots,rad,geometry='uniform',**options):
    # Builds 3D spots uniformly randomly distributed inside a ball. 
    # Algorithm: 
    # 1) generate spots uniformly randomly distributed inside a cube
    # 2) if there are any spots that don't fit our conditions, regenerate them
    # 3) repeat until all spots fit our conditions
    # NB: not very efficient, but still very fast for our purposes (~10^3-10^5 spots)
    #
    # Nspots - how many spots to build
    # rad - sphere radius
    # geometry:
    # * 'uniform' (default): spots inside a sphere of radius == rad
    # * 'shell-only': - spots inside a sphere of radius == rad with hollow center of radius = shell_thickness
    # * 'core-only':  - spots inside a sphere of radius == rad with hollow shell of radius = shell_thickness 
    
    spots = SpotArray()
    spots.x = np.zeros(Nspots)
    spots.y = np.zeros(Nspots)
    spots.z = np.zeros(Nspots)
    spots.t_unav = np.zeros(Nspots)
    bad_spots = np.full(Nspots, True)
    
    while np.any(bad_spots):
        N_bad_spots = np.sum(bad_spots)
        spots.x[bad_spots] = ((np.random.rand(N_bad_spots)-0.5) * 2) * rad
        spots.y[bad_spots] = ((np.random.rand(N_bad_spots)-0.5) * 2) * rad
        spots.z[bad_spots] = ((np.random.rand(N_bad_spots)-0.5) * 2) * rad
        r2 = np.power(spots.x,2) + np.power(spots.y,2) + np.power(spots.z,2)

        if geometry in ['uniform']:
            bad_spots = r2 > rad*rad
        elif geometry in ['shell-only']:
            shell_thickness = options.get['shell_thickness']
            bad_spots = np.logical_or( (r2 > rad*rad), (r2 < (rad-shell_thickness)*(rad-shell_thickness)) )
        elif geometry in ['core-only']:
            shell_thickness = options.get['shell_thickness']
            bad_spots = r2 > (rad-shell_thickness)*(rad-shell_thickness)
        else:
            raise Exception('Unknown shell geometry')

    return spots

def build_spots_spherical_surface(Nspots,rad):
    # Builds 3D spots uniformly randomly distributed on a spherical surface. 
    #
    # Nspots - how many spots to build
    # rad - sphere radius
    
    u = np.random.rand(Nspots)
    v = np.random.rand(Nspots)
    theta = 2*np.pi*u
    phi = np.arccos(2*v - 1)
    
    spots = SpotArray()
    spots.x = rad * np.cos(theta) * np.sin(phi)
    spots.y = rad * np.sin(theta) * np.sin(phi)
    spots.z = rad * np.cos(phi)
    spots.t_unav = np.zeros(Nspots)

    return spots


# In[26]:


def plot_geometry(donors,acceptors):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    get_ipython().run_line_magic('matplotlib', 'inline')
    
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(donors.x,donors.y,donors.z, c='b')
    ax.scatter(acceptors.x,acceptors.y,acceptors.z, c='r')


# In[30]:


def get_user_input(input_filename):
    # Reads init_xml.xml and returns a dictionary of input parameters
    # User can input what he wants or just press enter to use default value
    # Each input parameter is itself a dictionary with 'description', 'type' and 'value' keys
    print('Input parameters for simulation (use dot as decimal separator); press Enter if you want to keep default parameter value')

    input_xml_tree = xmlET.parse(input_filename)
    input_xml_tree_root = input_xml_tree.getroot()
    input_params = {}

    for child in input_xml_tree_root:
        if child.attrib['t'] == 'float':
            input_params[child.tag] = {'value': float(child.text), 'description': child.attrib['desc'], 'type': child.attrib['t']}
        elif child.attrib['t'] == 'int':
            input_params[child.tag] = {'value': int(child.text), 'description': child.attrib['desc'], 'type': child.attrib['t']}

        print('Input '+input_params[child.tag]['description']+' (default: '+str(input_params[child.tag]['value'])+'):')
        user_input = input('>>> ')
        if user_input != '':
            if input_params[child.tag]['type'] == 'float':
                input_params[child.tag]['value'] = float(user_input)
            elif input_params[child.tag]['type'] == 'int':
                input_params[child.tag]['value'] = int(user_input)

    return input_params, input_xml_tree


# In[119]:


def simulate_single_particle(input_params,sim_type,verbose_output=False):
    # assign input params to easily-readable variables
    rad = input_params['UCNP_radius']['value']
    dist = input_params['dye_ucnp_distance']['value']
    donorpercentage = input_params['emitter_doping_percentage']['value']
    nacceptors = input_params['number_of_dyes']['value']
    frad = input_params['FRET_radius']['value']
    eflux = input_params['exciton_flux']['value']
    donorlifetime = input_params['donor_lifetime']['value']
    acceptorlifetime = input_params['acceptor_lifetime']['value']
    totaltime = input_params['total_sim_time']['value']

    if (sim_type == 'h') or (sim_type == ''):
        dopedvolume = 4 / 3 * math.pi * math.pow(rad,3)                 # find UCNP volume
        Ndonors = int(donorpercentage / 100 * dopedvolume / 0.10744)    # find number of donor ions, constant is the volume of beta-NaYF4 unit cell in nm3
    elif (sim_type == 'c'):
        shell_thickness = input_params['undoped_shell_thickness']['value']
        dopedvolume = 4 / 3 * math.pi * math.pow(rad-shell_thickness,3)  # find how much UCNP volume is doped
        Ndonors = int(donorpercentage / 100 * dopedvolume / 0.10744)    # find number of donor ions, constant is the volume of beta-NaYF4 unit cell in nm3
    elif (sim_type == 's'):
        shell_thickness = input_params['doped_shell_thickness']['value']
        ucnpvolume = 4 / 3 * math.pi * math.pow(rad,3)                   # find UCNP volume
        undopedvolume = 4 / 3 * math.pi * math.pow(rad-shell_thickness,3) # find how much volume is undoped
        dopedvolume = ucnpvolume - undopedvolume                         # find how much UCNP volume is doped
        Ndonors = int(donorpercentage / 100 * dopedvolume / 0.10744)     # find number of donor ions, constant is the volume of beta-NaYF4 unit cell in nm3

    # -------------------------------- STEP 1: build the coords of donors and acceptors ----------------------------
    # build an array of randomly distributed donor spots (donor ions inside UCNP sphere/core/shell)
    print('Building donor coordinates: '+str(Ndonors)+' donors...')
    donors = []
    if (sim_type == 'h') or (sim_type == ''):
        donors = build_spots_inside_ball(Ndonors,rad,geometry='uniform')
    elif (sim_type == 'c'):
        donors = build_spots_inside_ball(Ndonors,rad,geometry='core-only',shell_thickness=shell_thickness)
    elif (sim_type == 's'):
        donors = build_spots_inside_ball(Ndonors,rad,geometry='shell-only',shell_thickness=shell_thickness)

    # build an array of randomly distributed acceptor spots (dyes around UCNP at a fixed distance from surface)
    print('Building acceptor coordinates: '+str(nacceptors)+' acceptors...')
    acceptors = build_spots_spherical_surface(nacceptors,rad+dist)

    # -------------------------------- STEP 2: build transfer probability matrix -----------------------------------
    print('Building transfer probability matrix...')
    # Probabilities of transfer for each donor-acceptor pair
    # are equal to R_fret^6 / (distance^6)
    frad6 = np.power(frad,6)  
    probmatrix = frad6 / (
        np.power(
            (   np.power(np.subtract.outer(donors.x, acceptors.x), 2)
              + np.power(np.subtract.outer(donors.y, acceptors.y), 2)
              + np.power(np.subtract.outer(donors.z, acceptors.z), 2)
            )
        , 2)
    )
    
    # -------------------------------- STEP 3: calculate exciton flux ----------------------------------------------
    print('Calculating exciton flux...')
    eflux = eflux           # temporary, rework for laser power, particle absorption and other parameters later

    # -------------------------------- STEP 4: generate exciton schedule -------------------------------------------
    print('Generating exciton schedule...')

    Nexcitons = math.floor(totaltime*eflux)
    exschedule = np.stack(
        (np.sort(np.random.rand(Nexcitons)) * totaltime,
         np.random.randint(0,Ndonors-1,size=Nexcitons))
        , axis=1
    )
    
    # -------------------------------- STEP 5: play over excitons --------------------------------------------------
    print('Starting playing excitons...')
    fretevents = 0                                  # initialize FRET event counter
    radevents = 0                                   # initialize radiative event counter
    inversedonorlifetime = 1 / donorlifetime        # to speed up the computations
    inverseacceptorlifetime = 1 / acceptorlifetime  # to speed up the computations
    for exciton in exschedule:
        current_time = exciton[0]                       # set time to current exciton time        
        current_donor = int(exciton[1])
        # refresh available donors:
        # count all donors that are still excited at this point in time        
        excdonorcounter = int(np.sum(donors.t_unav > current_time))
        # unexcite all donors which were excited but now are free
        donors.t_unav[donors.t_unav <= current_time] = 0
                
        # try to excite the donor which the current exciton hit 
        # if it's available - go through, otherwise this exciton is lost
        if (donors.t_unav[current_donor] == 0):        # if the donor we excite is available right now
            # compute when donor will release
            # initially no fret is probable
            totalprob = 0
            availableacceptorindices = []
            
            # deexcite all acceptors that should be available at this point in time
            acceptors.t_unav[acceptors.t_unav <= current_time] = 0
            totalprob = np.sum(probmatrix[current_donor][acceptors.t_unav==0])
            # store the indices of available acceptors
            # note: for a 1d array np.where returns a single-element tuple, so we must access its first element
            availableacceptorindices = np.where(acceptors.t_unav==0)[0]
            
            # calculate when the donor will release (see eq. 8 and 9 in Corry et al., 2005)
            actuallifetime = 1 / (inversedonorlifetime*(1+totalprob))
            releaseperiod = -1*actuallifetime*np.log(np.random.rand())
            # donor becomes unavailable for some time
            donors.t_unav[current_donor] = current_time + releaseperiod

            # now let's see where it transfers the energy: it can luminesce or it can do RET
            # (lifetimes that we put in already include the non-radiative pathways)

            # we can choose from radiation (-1) and any available acceptors
            choices = np.concatenate(
                ( np.array([-1]),
                 availableacceptorindices )
            )
            # we get weights for choices out of probmatrix, weights are normalized by Tt/Td (see eq. 10 in Corry et al.)
            Tt_Td = actuallifetime / donorlifetime                        
            weights = np.concatenate(
                ( np.array([Tt_Td]), 
                 probmatrix[current_donor][availableacceptorindices]*Tt_Td )
            )

            # choose between a radiative and FRET event with appropriate weights
            reschoice = int(np.random.choice(choices, p=weights))
            if (reschoice == -1):                                   # we have a radiative event!
                radevents += 1                                      # increase radiative counter
                toteff = fretevents / (fretevents + radevents)
                if verbose_output:
                    print(
                        ('{num:03d}'.format(num=excdonorcounter))+' donors excited; '+
                        'Total eff. so far: '+'{:.3f}'.format(toteff*100)+'%  '+
                        'Current: radiative event: from D #'+('{num:03d}'.format(num=current_donor))+
                        '             ', 
                        end='\r'
                    )
            else:                                                   # we have a FRET event!
                fretevents += 1                                     # increase FRET counter
                toteff = fretevents / (fretevents + radevents)
                # make acceptor unavailable
                acceptors.t_unav[reschoice] = current_time + (-1*inverseacceptorlifetime*np.log(np.random.rand()))     
                if verbose_output:
                    print(
                        ('{num:03d}'.format(num=excdonorcounter))+' donors excited; '+
                        'Total eff. so far: '+'{:.3f}'.format(toteff*100)+'%  '+
                        'Current: FRET event:      from D #'+('{num:03d}'.format(num=current_donor))+
                        ' to A #'+('{num:03d}'.format(num=reschoice)), 
                        end='\r'
                    )

    toteff = fretevents / (fretevents + radevents)
    print('Single simulation finished: FRET efficiency = '+'{:.3f}'.format(toteff*100)+'%')
    return toteff

def output_sim_results_to_file(avgtoteff,stdevtoteff,nsim):
    resfile = open('simulation_results.txt', 'w')
    resfile.write('Simulation ended successfully:')
    resfile.write('Total multidye FRET efficiency: '+'{:.2f}'.format(avgtoteff*100)+'%, stdev '+'{:.2f}'.format(stdevtoteff*100)+'% for '+('{num:01d}'.format(num=nsim))+' simulations\n')
    resfile.close()


# In[122]:


print('-----               UCNP Monte Carlo simulator, v1                   -----')
print('----- If you are using this program, please cite Dukhno et al., Nanoscale 2017 -----')
print('Choose type of particles to simulate: h (homogeneous, default), c (core-doped) or s (shell-doped)')
sim_type = input('>>> ')
input_params = {}
if (sim_type == 'h') or (sim_type == ''):
    # input_params, input_xml_tree = get_user_input('init_xml.xml')
    input_params, input_xml_tree = get_user_input('init_xml_tm.xml')
elif (sim_type == 'c'):
    input_params, input_xml_tree = get_user_input('init_coreshell_inner_xml.xml')
elif (sim_type == 's'):
    input_params, input_xml_tree = get_user_input('init_coreshell_outer_xml.xml')
else:
    quit()

# output a parameter file with user's parameters
input_xml_tree.write('last_simulation_parameters.xml')
# get number of single simulation repeats
nsim = input_params['simulation_repeats']['value']
# make a list to store fret efficiencies from all simulations
freteffs = []

for sim in range(0,nsim):
    print('----- Simulation #'+('{num:03d}'.format(num=sim))+' -----')
    toteff = simulate_single_particle(input_params,sim_type)
    freteffs.append(toteff)
    print('\nTotal multidye FRET efficiency: '+'{:.3f}'.format(toteff*100)+'% '+'for sim # '+('{num:03d}'.format(num=sim)))

avgtoteff = statistics.mean(freteffs)
stdevtoteff = statistics.stdev(freteffs)

# plot_geometry(donors,acceptors)
print('\n\nTotal multidye FRET efficiency: '+'{:.2f}'.format(avgtoteff*100)+'%, stdev '+'{:.2f}'.format(stdevtoteff*100)+'% for '+('{num:01d}'.format(num=nsim))+' simulations')
output_sim_results_to_file(avgtoteff,stdevtoteff,nsim)
input('\n\nAll simulations done, press Enter to exit.')


# In[126]:


get_ipython().run_line_magic('prun', 'simulate_single_particle(input_params,sim_type,verbose_output=False)')

