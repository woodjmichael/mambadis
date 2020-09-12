# mamba.py
# Fast economic dispatch simulation of PV-battery-generator microgrids

# python 3.7

__author__ = "Michael Wood"
__email__ = "michael.wood@mugrid.com"
__copyright__ = "Copyright 2020, muGrid Analytics"
__version__ = "6.8"

#
# Versions
#

#   6.8 - first peak shaving with battery charging from grid, plots tweak
#   6.7 - merge failed, re-comitting known-good v6.6
#   6.6 - import demand targets, debug printouts, dispatch plot tweaks, quickstart input data conventions
#   6.5 - new demand targets from linked results 090820
#   6.4 - demand debug, plot tweaks, no weekend arb in peak shaving, vectors output file now has datetime in col 1
#   6.3 - new output directory, change "output_" file to "resilience_", change "superloop_" filename to "resilience_superloop_", "add "sim_meta_" file, change which files are output when
#   6.2 - variable soc0 resilience simulation uses soc 35040 from previous utility sim (automatically checks for vectors file in ./Data/Output with appropriate name), update diesel fuel curves
#   6.1 - **change filename**, peak shaving includes monthly demand targets (hard coded), switch to arbitrage on weekends
# 6.0 - working peak shaving with annual target, new load line on plot
#   5.17 - split simulate_outage into _resilience and _utility_on (-sim arg), very simple utility_on dispatch, simulate_resilience almost direct from commit d36a4ca, no min/max soc, change program args, remove runtime from output file
#   5.16 - fixed soc_update and plots, some if/else dispatch logic, prevent crazy small bat.P vals, new debug and dispatch errors
#   5.15 - annual on-grid dispatch settings: RT batt efficiecny to 0.95, days to 365, simulation_interval to 8760
#   5.14 - soc_max and _min now available within Run options
#   5.13 - add grid_online switch and rework dispatch strategy for grid_online, batt charges only from PV for now, add smart_charging_on switch and keep off
#   5.12 - add battery hours arg
#   5.11 - adjust superloop output filename (significant digits on params)
#   5.10 - add dummy parameter to make output filenames unique during parallel runs
#   5.9 - fix version number
#   5.8 - bug: fuel curve A for 35-50 kW was 10x too high, add new fuel curve for propane 150 kW
#   5.7 - MG: load scaling factor, more superloop dims (independent battery power, generator size), custom sim interval, output file tweaks
#   5.6 - bug: "../" instead of "./"
#   5.5 - batt power varies 1:1 with energy unless -bp [kW], named output file for every loop in superloop, code version
#   5.4 - bug fix: for super small load-gen imbalance code would declare failure, microgrid class
#   5.3 - varying batt power doesn't do much, leave out for now
#   5.2 - for superloop batt power varies to always be 1h cap
#   5.1 - bug where superloop output fails
# 5.0 - superloop runs a matrix of battery energy and PV scale (oversize) factors, output and superloop files
#   4.6 - keep mambadis.py in main dir, output resilience conf and ttff
#   4.5 - fixed help, -sk skip ahead arg, new folder organization, new filenames for clarity
#   4.4 - help is broken (removed), add pv scaling factor
#   4.3 - help
#   4.2 - min on/off time was causing early microgrid failures
#   4.1 - change dispatch when gen tank is empty
# 4.0 - real CC dispatch including min on/off time, nonzero batt efficiency, need to clean up algo & "grid" aspect
#   3.11 - fix bug where simulation dispatch vector showed wrong P_pv
#   3.10 - simulation dispatch vector now output to csv
#   3.9 - ** PV 15 min interval data tested **, new filename
#   3.8 - can accept 15 min interval PV data (not fully tested)
#   3.7 - add _nf to numpy float and _ni to numpy int array variable names, _dt to datetime variables, remove some implicit casts
#   3.6 - simplify offsets, fuel curve lookup table, some small things
#   3.5 - measure runtime, command line run options overwrite defaults
#   3.4 - ** run from command line **, output to csv filename
#   3.3 - rename to mambadis.py ("fast snake" dispatch)
#   3.2 - now actually tested on python 3, pass sitename to load function, HR fire broken
#   3.1 - ** SOC_0 = 1.0 ***, load stats, tested on HR Adult Center (manual changes)
# 3.0 - manual changes to work with FISH data, PV timeseries wrap around 12/31
#   2.3 - switch to python 3.7 for dev, note that input data tested with HR Fire
#   2.2 - simulate_outage() works for many loops but "indexing" error, can't pre-allocate load. pv. bat. etc
#   2.1 - include gen fuel calc and tank size (not working now)
# 2.0 - "test_vX.py" changed to "lilcc_vX.py"
#   1.3 - no more "if LSimbalance>0" logic (not needed)
#   1.2 - "all" load/pv data as well as L-length data, tested on different t0's
#   1.1 - match up load data to PV vector, PV data in kW
# 1.0 - separate synthetic and real data
#   0.1 - load in real PV and load data
# 0.0 - first attempt with synthetic data (no rev number in file)

################################################################################
#
# Modules
#
################################################################################

import sys, os, platform
import numpy as np
import matplotlib.pyplot as plt
import csv
import datetime as dt
import time

################################################################################
#
# Classes
#
################################################################################

#
# Load and Solar Data (inputs)
#

class DataClass:
    def __init__(me, tstep, length):
        me.timestep = tstep              # units of seconds
        me.offset = 0                    # in case load and pv data don't begin at same time
        me.datetime = []
        me.P_kw_nf = np.zeros((length,), dtype=float)
        me.onlineTime_h_ni = np.zeros((length,),dtype=int)
        me.time_to_grid_import_h_nf = np.zeros((length,), dtype=float) # [h]
        me.soc_nf = np.zeros((length,), dtype=float) # [h]
        me.code_runtime_s = 0

    def clear(me):
        me.datetime = []
        me.P_kw_nf.fill(0)
        me.time_to_grid_import_h_nf.fill(0)

#
# Generator
#

class GenClass:
    def __init__(me, Pn_kw, a, b, tankCap, tstep, length):
        me.timestep = tstep                             # [s]
        me.dpph = 3600./tstep
        me.tstepHr = tstep/3600.
        me.Pn_kw = Pn_kw
        me.fuelcurve_Acoeff = a                         # [gal/h/kw]
        me.fuelcurve_Bcoeff = b                         # [gal/h]
        me.fuelcurve_Ainv = 1/a
        me.fuelTankCapacity = tankCap                   # [gal]
        me.P_kw_nf = np.zeros((length,), dtype=float)
        me.fuelConsumed = 0                              # [gal]
        me.min_on_time = 1                              # [timesteps]
        me.min_off_time = 1                             # [timesteps]
        me.on_time = me.min_on_time
        me.off_time = me.min_off_time
        me.prev_power = 0

    def clear(me):
        me.P_kw_nf.fill(0)
        me.fuelConsumed = 0


    def power_request(me, i, p_req):
        if (p_req > 0) & (p_req < 0.001):
            p_req = 0

        if (p_req > 0) & me.cool_down_complete():
            p_final = min(p_req,  me.Pn_kw, me.Pmax_tank())
        elif (p_req == 0) & me.warm_up_complete():
            p_final = 0.
        else:
            p_final = 0.
        me.on_off_counter(p_final)
        me.P_kw_nf[i] = np.array([p_final])
        me.fuel_calc(p_final)
        return p_final

    def Pmax_tank(me):
        fuel_remaining = max(me.fuelTankCapacity - me.fuelConsumed,0)
        if fuel_remaining:
            power = (fuel_remaining*me.dpph - me.fuelcurve_Bcoeff)*me.fuelcurve_Ainv
        else:
            power = 0
        return power

    def fuel_calc(me, power):
        if power:
            fuel_delta = (power*me.fuelcurve_Acoeff + me.fuelcurve_Bcoeff)*me.tstepHr
            me.fuelConsumed += fuel_delta

    def on_off_counter(me, now_power):
        if me.prev_power:
            if now_power:
                if (me.on_time < me.min_on_time): me.on_time = me.on_time + 1
            else:
                me.off_time = 1
        else:
            if now_power:
                me.on_time = 1
            else:
                if (me.off_time < me.min_off_time): me.off_time = me.off_time + 1
        me.prev_power = now_power

    def warm_up_complete(me):
        ret = me.on_time >=  me.min_on_time
        return ret

    def cool_down_complete(me):
        ret = me.off_time >=  me.min_off_time
        return ret

    def tank_empty(me):
        return  (me.fuelTankCapacity - me.fuelConsumed) < 0.001


#
# Grid
#

class GridClass:
    def __init__(me, Pn_kw, length):
        me.Pn_kw = Pn_kw
        me.P_kw_nf = np.zeros((length,), dtype=float)
        me.time_to_import = 0
        me.import_on = 0
        me.offlineCounter = 0

    def clear(me):
        me.P_kw_nf.fill(0)
        me.time_to_import = 0
        me.import_on = 0

    def power_request(me, index, p_req):
        me.timer_tick()
        p_final = 0.
        if p_req > 0:
            p_final = np.minimum(p_req,me.Pn_kw)
            me.import_on = 1
        elif p_req < 0:
            p_final = np.maximum(p_req, -me.Pn_kw)
        me.P_kw_nf[index] = p_final
        return p_final

    def timer_tick(me):
        if not me.import_on:
            me.time_to_import += 1

#
# Microgrid
#

class MicrogridClass:
    def __init__(me):
        me.failure = 0
        me.time_to_failure = 0
        me.demand_targets = []

    def clear(me):
        me.time_to_failure = 0
        me.failure = 0

    def timer_tick(me):
        if not me.failure:
            me.time_to_failure += 1

    def failed(me):
        me.failure = 1


#
# Battery
#

class BattClass:
    def __init__(me, Pn_kw, En_kwh, soc0, timestep, length):
        me.Pn_kw = Pn_kw                            # pos:dischg, neg:chg
        me.En_kwh = En_kwh
        me.soc0 = soc0
        me.soc_prev = soc0
        me.timestep = timestep                      # units of seconds
        me.P_kw_nf = np.zeros((length,), dtype=float)
        me.soc_nf = np.zeros((length,), dtype=float)
        me.soc_flag = 0                             # binary
        me.eff_chg = 0.95
        me.eff_dischg = 0.95

    def clear(me):
        me.soc_prev = me.soc0
        me.P_kw_nf.fill(0)
        me.soc_nf.fill(0)
        me.soc_flag = 0

    def power_request(me, index, p_req):
        if p_req > 0:
            p_final = min(p_req, me.Pn_kw, me.P_max_soc())
        elif p_req < 0:
            p_final = max(p_req, -me.Pn_kw, me.P_min_soc())
        else:
            p_final = 0.
        if abs(p_final) < 0.001: p_final = 0
        me.soc_nf[index] = me.soc_update(p_final)
        me.P_kw_nf[index] = p_final
        return p_final

    def soc_update(me,p_soc):
        if p_soc > 0:
            soc_new = me.soc_prev - (p_soc * me.timestep/3600.0 / me.eff_dischg)/me.En_kwh
        elif p_soc < 0:
            soc_new = me.soc_prev - (p_soc * me.timestep/3600.0 * me.eff_chg)/me.En_kwh
        else:
            soc_new = me.soc_prev

        me.soc_prev = soc_new
        return soc_new

    def P_max_soc(me):
        return me.soc_prev * me.En_kwh * (3600.0/me.timestep) / me.eff_dischg

    def P_min_soc(me):
        return -1*(1 - me.soc_prev) * me.En_kwh * (3600.0/me.timestep) * me.eff_chg

    def empty(me,i):
        return (me.soc_nf[i] < 0.001)

    def full(me,i):
        return (me.soc_nf[i] > 0.999)

    def over_half(me,i):
        return (me.soc_nf[i] > 0.5)

    def set_soc0(me,soc):
        me.soc0 = soc
        me.soc_prev = soc


#
# Faults
#

class FaultClass:
    def __init__(me):
        me.mainloop = 0
        me.energybalance = 0
        me.index = 0
        me.fuelCurveCoeffs = 0
        me.dispatch = 0

    def print_faults(me):

        #print("\033[1;31;40m Bright Green  \n")
        if err.mainloop:
            print('')
            print('\033[1;31;1m error main loop (qty {:d})'.format(err.mainloop))
        if err.energybalance:
            print('')
            print('\033[1;31;1m error energy balance (qty {:d})'.format(err.energybalance))
        if err.index:
            print('')
            print('\033[1;31;1m error indexing (qty {:d})'.format(err.index))
        if err.fuelCurveCoeffs:
            print('')
            print('\033[1;31;1m error fuel curve coeffs (qty {:d})'.format(err.fuelCurveCoeffs))
        if err.dispatch:
            print('')
            print('\033[1;31;1m error with dispatch (qty {:d})'.format(err.dispatch))

    def main_loop(me):
        me.mainloop += 1

    def energy_balance(me):
        me.energybalance += 1

    def indexing(me):
        me.index += 1

    def gen_fuel_coeffs(me):
        me.fuelCurveCoeffs += 1

    def dispatch(me):
        me.dispatch += 1


################################################################################
#
# Functions
#
################################################################################

#
# Find load-pv data offset
#

# find the offset between start of load data and PV data
def find_load_pv_offset():
    Load_start_of_data = load_all.datetime[0]
    Load_start_of_data_doy = Load_start_of_data.timetuple().tm_yday     # day of year
    Load_start_of_data_hoy = (Load_start_of_data_doy-1)*24 + Load_start_of_data.hour    # hour of year
    pv_all.offset = Load_start_of_data_hoy

#
# Fuel curve coefficients
#

def lookup_fuel_curve_coeffs(power, gen_fuel_propane):
    # diesel
    # coeffs = [fuel_curve_coeff_A, fuel_curve_coeff_B]
    if not gen_fuel_propane:
        if power < 25:      coeffs = [0.0680, 0.250] # 20 kW
        elif power < 35:    coeffs = [0.0720, 0.750] # 30 kW
        elif power < 50:    coeffs = [0.0810, 0.750] # 40 kW
        elif power < 67:    coeffs = [0.0660, 0.850] # 60 kW
        elif power < 87:    coeffs = [0.0656, 1.050] # 75 kW
        elif power < 120:   coeffs = [0.0644, 0.950] # 100 KW
        elif power < 280:   coeffs = [0.0648, 1.350] # 200 KW
        elif power < 520:   coeffs = [0.0655, 2.050] # 400 KW
        else:
            coeffs = [0.0655, 2.05]
            err.gen_fuel_coeffs()

    # propane
    elif gen_fuel_propane:
        if power < 180:
            coeffs = [0.1441,0.665]
        else:
            coeffs = [0.1441,0.665]
            err.gen_fuel_coeffs()
    else:
        coeffs = [0.0655, 2.050]
        err.gen_fuel_coeffs()

    return coeffs


#
# Synthetic data creation
#

def create_synthetic_data():
    pv_oneday = 100*np.array([0., 0., 0., 0., 0., 0., 1., 2., 3., 4., 5., 6., 7., 6., 5., 4., 3., 2., 1., 0., 0., 0., 0., 0.])
    load_oneday = 100*np.array([1., 1., 1., 1., 1., 1., 5., 3., 1., 1., 4., 2., 1., 1., 1., 1., 6., 8., 12., 8., 5., 3., 2., 1.])
    pv.P_kw_nf = np.concatenate([pv_oneday,pv_oneday,pv_oneday])
    load.P_kw_nf = np.concatenate([load_oneday,load_oneday,load_oneday])

#
# Import load data
#

def import_load_data(site, load_stats):
    filename = './Data/Load/' + site + '_load.csv'
    with open(filename,'r') as f:
        datacsv = list(csv.reader(f, delimiter=","))
        del datacsv[0]
        t = []
        p = []
        for line in datacsv:
            newtxt = line[0] #hr fire .split('P')[0][0:-1]
            newval = dt.datetime.strptime(newtxt, '%Y/%m/%d %H:%M')
            t.append(newval)
    load_all.datetime = t + t

    my_data = np.genfromtxt(filename, delimiter=',')
    md = load_scaling_factor * my_data[1:,1]
    load_all.P_kw_nf = np.concatenate((md,md),axis=0)

    if load_stats:
        print('max load [kw] = {:.1f}'.format(np.amax(load_all.P_kw_nf)))
        print('avg load [kw] = {:.1f}'.format(np.average(load_all.P_kw_nf)))
        print('min load [kw] = {:.1f}'.format(np.amin(load_all.P_kw_nf)))

#
# Import solar vector
#

def import_pv_data(site):

    if solar_data_inverval_15min:
        filename = './Data/Solar/' + site + '_solar_35040.csv'
    else:
        filename = './Data/Solar/' + site + '_solar.csv'

    with open(filename,'r') as f:
        datacsv = list(csv.reader(f, delimiter=","))
        del datacsv[0]
        t = []
        p = []

        for line in datacsv:
            if solar_data_inverval_15min:
                newtxt = line[0]
                newval = dt.datetime.strptime(newtxt, '%Y/%m/%d %H:%M')
            else:
                newtxt = line[1]
                newval = dt.datetime.strptime(newtxt, '%Y-%m-%d %H:%M:%S')
            t.append(newval)

    pv_all.datetime = t + t

    my_data = np.genfromtxt(filename, delimiter=',')
    pv = my_data[1:,-1]/1000.  # W to kW

    i=0
    for row in pv:
        if row > 0:
            pass
        else:
            pv[i] = 0
        i += 1

    pv_all.P_kw_nf = np.concatenate((pv,pv),axis=0)
    pv_all.P_kw_nf = pv_scaling_factor * pv_all.P_kw_nf

#
# Import 35040 of soc from annual simulation
#

def import_soc_35040(site):
    vals = np.genfromtxt(soc_filename, delimiter=',')
    dispatch_previous.soc_nf = np.concatenate((vals[1:],vals[1:]),axis=0)


#
# Import 12 demand target values
#

def import_demand_targets(site):
    filename = './Data/Demand targets/' + site + '_demand_targets.csv'
    vals = np.genfromtxt(filename, delimiter=',')
    demand_targets.monthly = vals[1:] / demand_throttle

#
# Simulate resilience (from commit d36a4ca)
#

def simulate_resilience(t_0,L):

    #
    # Slice data
    #

    # temporarily store small chunk "all load" vector in "load"

    # where in "all load" data this run will begin
    m_0 = t_0 + load_all.offset    # offset usually 0

    # where in "all load" data this run will end
    m_end = m_0 + L

    load.P_kw_nf = load_all.P_kw_nf[m_0:m_end]
    load.datetime = load_all.datetime[m_0:m_end]

    if debug_indexing:
        print('LOAD')
        print('m_0={:d}'.format(m_0) + ' m_end={:d}'.format(m_end))
        print(load.datetime[0])
        print(load.datetime[-1])
        print(load.P_kw_nf.size)
        print(len(load.datetime))
        print('')

    # temporarily store small chunk "all pv" vector in "pv"



    if solar_data_inverval_15min:
        n_0 = t_0 + pv_all.offset
        n_end = n_0 + L

    else:
        # where in "all PV" data this run will begin
        n_0 = t_0//4 + pv_all.offset      # offset usually 0

        # where in "all load" data this run will end
        n_end = n_0 + L//4                # offset usually 0



    pv.P_kw_nf = pv_all.P_kw_nf[n_0:n_end]
    pv.datetime = pv_all.datetime[n_0:n_end]

    if debug_indexing:
        print('PV')
        print('n_0={:d}'.format(n_0) + ' n_end={:d}'.format(n_end))
        print(pv.datetime[0])
        print(pv.datetime[-1])
        print(pv.P_kw_nf.size)
        print(len(pv.datetime))
        print('')

    # check indexing
    # beginning and ending date and hour should match between load and pv
    if (load.datetime[0].day    - pv.datetime[0].day):
        err.indexing()
    if (load.datetime[-1].day   - pv.datetime[-1].day):
        err.indexing()
    if (load.datetime[0].hour   - pv.datetime[0].hour):
        err.indexing()
    if (load.datetime[-1].hour  - pv.datetime[-1].hour):
        err.indexing()


    #
    # Algorithm
    #

    # begin battery with the soc from a previous year-long utility-on simulation
    if vary_soc: bat.set_soc0(dispatch_previous.soc_nf[t_0])

    if debug_res: print('soc={:.2f}'.format(bat.soc0))

    chg = 0


    for i in range(L):

        if solar_data_inverval_15min:
            # only increment i_pv every 4 i-increments
            i_pv = i
        else:
            i_pv = i//4

        LSimbalance = load.P_kw_nf[i]      -   pv.P_kw_nf[i_pv]   # load-solar imbalance

        if not chg and not gen.tank_empty():

            battpower = bat.power_request(i,LSimbalance)

            LSBimbalance =  LSimbalance     -   battpower       # load-solar-batt imbalance

            genpower = gen.power_request(i,LSBimbalance)

            # at low soc assume to turn gen on
            if bat.empty(i):
                LSGimbalance = LSimbalance - gen.Pn_kw
                battpower = bat.power_request(i,LSGimbalance)
                LSBimbalance = LSimbalance - battpower
                genpower = gen.power_request(i,LSBimbalance)
                chg = 1

        if chg and not gen.tank_empty():

            LSGimbalance = LSimbalance - gen.Pn_kw
            battpower = bat.power_request(i,LSGimbalance)
            LSBimbalance = LSimbalance - battpower
            genpower = gen.power_request(i,LSBimbalance)

            # turn off gen when batt full or half full and PV turns on
            if bat.full(i) or (bat.over_half(i) and (pv.P_kw_nf[i_pv] > 0)):
                battpower = bat.power_request(i,LSimbalance)
                LSBimbalance = LSimbalance - battpower
                genpower = gen.power_request(i,LSBimbalance)
                chg = 0

        if gen.tank_empty():
            chg = 0
            battpower = bat.power_request(i,LSimbalance)
            LSBimbalance = LSimbalance - battpower
            genpower = gen.power_request(i,0)

        LSBGimbalance = LSimbalance - battpower  -   genpower        # load-solar-batt-gen imbalance

        # check if load is fully served
        if(LSBGimbalance > 0.1):
            microgrid.failed()
        else:
            microgrid.timer_tick()

        gridpower = grid.power_request(i,LSBGimbalance)



        if gridpower <= 0:
            grid.offlineCounter += 1                        # time that microgrid services load


        # check energy balance
        if np.absolute((LSimbalance - bat.P_kw_nf[i] - gen.P_kw_nf[i] - grid.P_kw_nf.item(i))) > 0.001:
            err.energy_balance()

    time_to_grid_import = microgrid.time_to_failure/(3600./load.timestep)

    # vectors
    if output_vectors:
        filename = output_dir + '/vectors_{}.csv'.format(filename_param)
        with open(filename, 'w', newline='') as file:
            output = csv.writer(file)
            output.writerow(['time','load','pv','b_kw','b_soc','gen','grid','diff'])
            for i in range(L):
                if solar_data_inverval_15min:
                    i_pv = i
                else:
                    i_pv = i//4 # only increment i_pv every 4 i-increments
                l=load.P_kw_nf.item(i)
                p=pv.P_kw_nf.item(i_pv)
                b=bat.P_kw_nf.item(i)
                s=bat.soc_nf.item(i)
                g=gen.P_kw_nf.item(i)
                G=grid.P_kw_nf.item(i)
                d=l-p-b-g-G
                output.writerow([load.datetime[i],l,p,b,s,g,G,d])

    if debug: print('checksum: {:.1f}'.format(np.sum(bat.P_kw_nf)))

    return time_to_grid_import


#
# Simulate utility online (from branch soc_analysis1 + mw mods)
#

def simulate_utility_on(t_0,L):

    #
    # Slice data
    #


    # temporarily store small chunk "all load" vector in "load"

    # where in "all load" data this run will begin
    m_0 = t_0 + load_all.offset    # offset usually 0

    # where in "all load" data this run will end
    m_end = m_0 + L

    load.P_kw_nf = load_all.P_kw_nf[m_0:m_end]
    load.datetime = load_all.datetime[m_0:m_end]

    if debug_indexing:
        print('LOAD')
        print('m_0={:d}'.format(m_0) + ' m_end={:d}'.format(m_end))
        print(load.datetime[0])
        print(load.datetime[-1])
        print(load.P_kw_nf.size)
        print(len(load.datetime))
        print('')

    # temporarily store small chunk "all pv" vector in "pv"
    if solar_data_inverval_15min:
        n_0 = t_0 + pv_all.offset
        n_end = n_0 + L
    else:
        # where in "all PV" data this run will begin
        n_0 = t_0//4 + pv_all.offset      # offset usually 0

        # where in "all load" data this run will end
        n_end = n_0 + L//4                # offset usually 0



    pv.P_kw_nf = pv_all.P_kw_nf[n_0:n_end]
    pv.datetime = pv_all.datetime[n_0:n_end]

    if debug_indexing:
        print('PV')
        print('n_0={:d}'.format(n_0) + ' n_end={:d}'.format(n_end))
        print(pv.datetime[0])
        print(pv.datetime[-1])
        print(pv.P_kw_nf.size)
        print(len(pv.datetime))
        print('')

    # check indexing
    # beginning and ending date and hour should match between load and pv
    if (load.datetime[0].day    - pv.datetime[0].day):
        err.indexing()
    if (load.datetime[-1].day   - pv.datetime[-1].day):
        err.indexing()
    if (load.datetime[0].hour   - pv.datetime[0].hour):
        err.indexing()
    if (load.datetime[-1].hour  - pv.datetime[-1].hour):
        err.indexing()




    #
    # Algorithm
    #

    for i in range(L):

        if peak_shaving:
            month = load.datetime[i].month
            demand_target = demand_targets.monthly[month-1]

            day_of_week = load.datetime[i].weekday()

            if day_of_week > 4: # saturday=5 sunday=6
                weekend = 0 #1 # disable for now
            else:
                weekend = 0

        if solar_data_inverval_15min:
            # only increment i_pv every 4 i-increments
            i_pv = i
        else:
            i_pv = i//4

        LSimbalance = load.P_kw_nf[i]      -   pv.P_kw_nf[i_pv]   # load-solar imbalance

        # arbitrage
        if not peak_shaving or (peak_shaving and weekend):
            battpower = bat.power_request(i,LSimbalance)
            LSBimbalance = LSimbalance - battpower
            gpower = grid.power_request(i, LSBimbalance)

        # peak shaving with charging batt from grid
        elif peak_shaving and grid_charging:

            battpower = bat.power_request(i,LSimbalance - demand_target)
            LSBimbalance = LSimbalance - battpower
            gpower = grid.power_request(i,LSBimbalance)

        # fringe case: peak shaving without charging batt from grid
        elif peak_shaving and not grid_charging:

            # demand is too high: dispatch battery
            if LSimbalance > demand_target:
                battpower = bat.power_request(i,LSimbalance - demand_target)
                LSBimbalance = LSimbalance - battpower
                gpower = grid.power_request(i,LSBimbalance)

            # demand is below threshold, do nothing
            elif LSimbalance <= demand_target and LSimbalance > 0:
                battpower = bat.power_request(i,0)
                LSBimbalance = LSimbalance - battpower
                gpower = grid.power_request(i,LSBimbalance)

            # excess solar: charge battery from solar
            elif LSimbalance <= 0:
                battpower = bat.power_request(i,LSimbalance)
                LSBimbalance = LSimbalance - battpower
                gpower = grid.power_request(i,LSBimbalance)

        else:
            print('error: to peak shave or not?')
            quit()

        LSBGimbalance = LSimbalance - battpower  -   gpower        # load-solar-batt-grid/gen imbalance

        # check energy balance
        if np.absolute((LSimbalance - bat.P_kw_nf[i] - gen.P_kw_nf[i] - grid.P_kw_nf.item(i))) > 0.001:
            err.energy_balance()

    # vectors
    if output_vectors:
        filename = output_dir + '/vectors_{}.csv'.format(filename_param)
        with open(filename, 'w', newline='') as file:
            output = csv.writer(file)
            output.writerow(['time','load','pv','b_kw','b_soc','gen','grid','diff'])
            for i in range(L):
                if solar_data_inverval_15min:
                    i_pv = i
                else:
                    i_pv = i//4 # only increment i_pv every 4 i-increments
                l=load.P_kw_nf.item(i)
                p=pv.P_kw_nf.item(i_pv)
                b=bat.P_kw_nf.item(i)
                s=bat.soc_nf.item(i)
                g=gen.P_kw_nf.item(i)
                G=grid.P_kw_nf.item(i)
                d=l-p-b-g-G
                output.writerow([load.datetime[i],l,p,b,s,g,G,d])

    # print a single vector
    if batt_vector_print:
        print('battery power output vector [kW]')
        for val in bat.P_kw_nf:
            print(val)
        print('-----batt vector done-------')
        print('\n\n')

    # print some helpful debug stats
    if debug_energy:
        sum_batdis = np.sum(np.clip(bat.P_kw_nf,0,10000))/4
        sum_batchg = np.sum(-1*np.clip(bat.P_kw_nf,-10000,0))/4

        sum_load = np.sum(load.P_kw_nf)/4
        sum_pv = np.sum(pv.P_kw_nf)/4
        sum_grid = np.sum(grid.P_kw_nf)/4

        soc_final = bat.soc_nf.item(L-1)
        delta_soc = bat.En_kwh * (1 - soc_final)

        max_grid = np.max(grid.P_kw_nf)

        print('sum load: {:.0f}'.format(sum_load))
        print('sum pv: {:.0f}'.format(sum_pv))
        print('sum grid: {:.0f}'.format(sum_grid))

        print('sum bat dis: {:.0f}'.format(sum_batdis))
        print('sum bat chg: {:.0f}'.format(sum_batchg))

        print('soc final: {:.5f}'.format(soc_final))
        print('delta soc [kwh]: {:.0f}'.format(delta_soc))

        print('max grid [kw]: {:.0f}'.format(max_grid))

    # print some helpful debug stats
    if debug_demand:
        print('** debug demand **')

        # output some stats
        peak_demand = np.max(load.P_kw_nf - pv.P_kw_nf - bat.P_kw_nf)
        peak_demand_index = np.argmax(load.P_kw_nf - pv.P_kw_nf - bat.P_kw_nf)


        print('\nperiod start: {}'.format(load.datetime[0]))
        print('period end:   {}'.format(load.datetime[-1]))
        print('period peak demand [kW]: {:.1f}'.format(peak_demand))
        print('period peak demand time: {}'.format(load.datetime[peak_demand_index]))
        print('period peak demand index [tstep]: {:d}'.format(peak_demand_index))
        print('\ndemand targets [kW]: {}'.format(str(demand_targets.monthly)))

        demand_ratchet = 0
        month = 1
        demands = []
        for i in range(L):
            if (load.datetime[i].month > month) or (i==(L-1)):
                demands.append(demand_ratchet)
                demand_ratchet = 0
                month += 1
            if grid.P_kw_nf[i] > demand_ratchet:
                demand_ratchet = grid.P_kw_nf[i]

        results.demands = np.array(demands)
        print('\ndemands (calendar) [kW]:')
        for val in results.demands:
            print(val)

        print('\ndemand errors [kW]:')
        for i in range(len(results.demands)):
            print (results.demands[i] - demand_targets.monthly[i])

    # print checksum
    if debug:
        print('\nchecksum: {:.3f}'.format(np.sum(bat.P_kw_nf)))

#
# Print help
#

def help_printout():
    print('\nmamba.py\nmuGrid Analytics LLC\nMichael Wood\nmichael.wood@mugrid.com')
    print('')
    print('Arguments are best issued in this order, with the values following directly after keys')
    print('e.g. % python mamba.py -s fish -bp 20 be 40 .. (etc)')
    print('')
    print('Typical Command Line Arguments')
    print(' Simulation type:        -sim [r=res | ua = utility arbitrage | up = utility peak shaving]    e.g. -sim r')
    print(' Site name:              -s  [sitename]              e.g. -s hradult')
    print(' Battery power:          -bp [power kW]              e.g. -bp 60')
    print(' Battery energy:         -be [energy kWh]            e.g. -be 120')
    print(' Generator power:        -gp [power kW]              e.g. -gp 60')
    print(' Generator tank:         -gt [size gal]              e.g. -gt 200')
    print('')
    print('Optional Command Line Arguments')
    print(' Run "n" simulations:    -r [n]                      e.g. -r 1   default=2920')
    print(' Dispatch vectors ON:    -v                          e.g. -v     default=OFF')
    print(' Battery vector ON:      -vb                           e.g. -vb    default=OFF')
    print(' Load stats ON:          --loadstats                 e.g. --loadstats     default=OFF')
    print(' Skip ahead "h" hours:   -sk [h]                     e.g. -sk 24 default=OFF')
    print(' Superloop enable:       -sl                         e.g. -sl    default=OFF')
    print(' Gen fuel is propane:    -gfp                        e.g. -gfp   default=OFF')
    print(' Days to simulate:       --days [days]               e.g. --days 3')
    print('     --days must come after --sim because it re-modifies some variables')
    print(' Battery depth of dischg:-bd [dod]                   e.g. -bd 0.95')
    print('     -bd must come after -be because it modifies battery energy')
    print('     careful: -bd just changes the battery energy, so soc will still be 0-100%')
    print(' Plots ON (option to plot normal or utility first):    --plots [ | u]                     e.g. --plots     default=OFF')
    print('')



################################################################################
#
# Main
#
################################################################################


t_script_begin_dt = dt.datetime.now()

err =   FaultClass()


#
# Run options
#

# data
site = 'badriverhc'                    # fish, hradult, (hrfire not working)
pv_scaling_factor = 1
load_scaling_factor = 1
filename_param = ''

# simulation
simulation_interval = 3            # hours between simulated outages
runs = 365*24//simulation_interval   # number of iterations
skip_ahead = 0                      # number of hours to skip ahead
solar_data_inverval_15min = 1
superloop_enabled = 0               # run matrix of battery energy and PV scale (oversize) factors
grid_online = 0
peak_shaving = 0
grid_charging = 0
demand_target = 0
days = 14                          # length of grid outage
L = days*24*4                     # length of simulation in timesteps
vary_soc = 0
weekend_arb_off = 0
demand_throttle = 1

# physical capacities
batt_power = 200.         # kw
batt_hrs = 0.           # kWh/kW
batt_energy = 580.       # kwh
dod = 1.                 # depth of discharge
gen_power = 0.           # kw
gen_tank = 100.          # gal
gen_fuel_propane = 0    # 1 = propane, 0 = diesel
batt_power_varies = 0  # batt power such that capacity = 1h

# outputs on/off
output_sim_meta = 1        # only this should be left on normally
output_resilience = 0
output_vectors = 0
plots_on = 0
load_stats = 0
debug = 0
debug_energy = 0
debug_demand = 0
debug_indexing = 0
debug_res = 0
batt_vector_print = 0


# command line run options override defaults
if len(sys.argv) > 1:

    for i in range(1, len(sys.argv)):

        if sys.argv[i] == '-s':
            site = str(sys.argv[i+1])

        elif sys.argv[i] == '-p':
            pv_scaling_factor = float(sys.argv[i+1])

        elif sys.argv[i] == '-ls':
            load_scaling_factor = float(sys.argv[i+1])

        elif sys.argv[i] == '-r':
            runs = int(sys.argv[i+1])

        elif sys.argv[i] == '-bp':
            batt_power = float(sys.argv[i+1])
            batt_power_varies = 0

        elif sys.argv[i] == '-bh':
            batt_hrs = float(sys.argv[i+1])
            batt_power_varies = 1

        elif sys.argv[i] == '-be':
            batt_energy = float(sys.argv[i+1])

        elif sys.argv[i] == '-bd':
            dod = float(sys.argv[i+1])
            batt_energy = batt_energy * dod

        elif sys.argv[i] == '-dt':
            demand_throttle = float(sys.argv[i+1])

        elif sys.argv[i] == '-gp':
            gen_power = float(sys.argv[i+1])

        elif sys.argv[i] == '-sim':
            sim = str(sys.argv[i+1])
            if sim == 'r':
                grid_online = 0
                simulation_interval = 3
                runs = 365*24//simulation_interval   # number of iterations
                days = 14                          # length of grid outage
                L = days*24*4                     # length of simulation in timesteps
                output_resilience = 1

            elif sim == 'ua':
                grid_online = 1
                gen_power = 0
                simulation_interval = 8760
                runs = 1
                days = 365
                L = days*24*4                     # length of simulation in timesteps
                output_vectors = 1

            elif sim == 'up':
                grid_online = 1
                grid_charging = 1
                gen_power = 0
                simulation_interval = 8760
                runs = 1
                days = 365
                L = days*24*4                     # length of simulation in timesteps
                peak_shaving = 1
                output_vectors = 1

            else:
                print('\033[1;31;1m fatality: no simulation type selected')
                quit()

        elif sys.argv[i] == '-gt':
            gen_tank = float(sys.argv[i+1])

        elif sys.argv[i] == '-v':
            output_vectors = 1

        elif sys.argv[i] == '-vb':
            batt_vector_print = 1

        elif sys.argv[i] == '--plots':
            plots_on = 1
            plot_pv_first = 1
            plot_grid_first = 0

            if i+1 < len(sys.argv):         # make sure there are more args
                type = str(sys.argv[i+1])

                if type == 'u':
                    plot_pv_first = 0
                    plot_grid_first = 1

        elif sys.argv[i] == '--loadstats':
            load_stats = 1

        elif sys.argv[i] == '--debug':

            debug = 1                       # basic debug only
            __version__ = __version__ + '_debug'

            if i+1 < len(sys.argv):         # make sure there are more args
                bug = str(sys.argv[i+1])

                if bug == 'energy':
                    debug_energy = 1

                elif bug == 'demand':
                    debug_demand = 1

                elif bug == 'indexing':
                    debug_indexing = 1

                elif bug == 'res':
                    debug_res = 1

        elif sys.argv[i] == '--days':
            days = int(sys.argv[i+1])
            L = days*24*4                     # length of simulation in timesteps

        elif sys.argv[i] == '-sk':
            skip_ahead = int(sys.argv[i+1])

        elif sys.argv[i] == '-sl':
            superloop_enabled = 1

        elif sys.argv[i] == '-slp':
            j = i + 6 + 1                # excpect 6 more args
            superloop_enabled = 1
            strs = sys.argv[i+1:j]
            pv_scale_vector = [float(i) for i in strs]

        elif sys.argv[i] == '-slb':
            j = i + 4 + 1               # excpect 4 more args
            superloop_enabled = 1
            strs = sys.argv[i+1:j]
            batt_energy_vector = [int(i) for i in strs]

        elif sys.argv[i] == '-gfp':
            gen_fuel_propane = 1

        elif sys.argv[i] == '-sv':
            filename_param = str(sys.argv[i+1])

        elif sys.argv[i] == '--help' :
            help_printout()
            quit()

#
# Create output directory
#
now = dt.datetime.now().strftime("%Y_%m_%d__%H_%M_%S")
output_dir = './Data/Output/' + site + '_' + str(now)
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
else:
    print('error: output directory already exists')

#
# Superloop
#

if superloop_enabled:
    pv_scale_vector = [1, 1.25, 1.5, 1.75, 2.0]
    load_scale_vector = [1]
    batt_power_vector = [30,60,90]
    batt_hrs_vector = [1]
    gen_power_vector = [0]
else:
    pv_scale_vector = [pv_scaling_factor]
    batt_energy_vector = [batt_energy]
    load_scale_vector = [load_scaling_factor]
    batt_hrs_vector = [batt_hrs]
    batt_power_vector = [batt_power]
    gen_power_vector = [gen_power]


pv_scale_vals= []
load_scale_vals= []
batt_energy_vals = []
batt_hrs_vals = []
batt_power_vals = []
gen_power_vals = []

max_ttff = []
avg_ttff = []
min_ttff = []

conf_72h = []
conf_336h = []

for load_scaling_factor in load_scale_vector:
    for pv_scaling_factor in pv_scale_vector:
        for batt_power in batt_power_vector:
            for batt_hrs in batt_hrs_vector:
                for gen_power in gen_power_vector:
                    [gen_fuelA, gen_fuelB] = lookup_fuel_curve_coeffs(gen_power, gen_fuel_propane)

                    if batt_power_varies:
                        batt_energy = batt_power * batt_hrs

                    #
                    # Data source
                    #

                    # synthetic data
                    #L = 72
                    #pv = DataClass(L)
                    #load = DataClass(L)
                    #create_synthetic_data()

                    #hr fire load_all =  DataClass(15.*60., 50788)  # timestep[s], hood river fire size
                    load_all =  DataClass(15.*60., 2*46333)  # timestep[s]
                    pv_all =    DataClass(60.*60., 2*8760)     # timestep[s]
                    results =   DataClass(3.*60.*60., runs)
                    dispatch_previous =  DataClass(15.*60., 2*46333)  # timestep[s]
                    demand_targets =  DataClass(15.*60., 2*46333)  # timestep[s]
                    import_load_data(site, load_stats)
                    import_pv_data(site)

                    # look for dispatch file containing soc 35040 - use if available
                    soc_filename = './Data/Dispatch/soc_' + site + '_35040.csv'
                    if not grid_online and os.path.isfile(soc_filename):
                        vary_soc = 1
                        import_soc_35040(site)

                    if peak_shaving:
                        import_demand_targets(site)


                    #
                    # some data prep
                    #

                    #
                    # Simulation Loop
                    #

                    for i in range(runs):
                        h = simulation_interval * i + skip_ahead
                        t0=h*4       # where to start simulation in "all load" vector

                        # start with fresh variables
                        # wish we could pre-allocate these, but it was causing a bug
                        #   even after calling .clear() on everything
                        load =  DataClass(  15.*60.,L)                 # timestep[s]
                        pv =    DataClass(  60.*60.,L)                   # timestep[s]
                        gen =   GenClass(   gen_power,gen_fuelA,gen_fuelB,gen_tank,15.*60.,L)   # kW, fuel A, fuel B, tank[gal], tstep[s]
                        bat =   BattClass(  batt_power,batt_energy,1,15*60.,L)      # kW, kWh, soc0 tstep[s]
                        grid =  GridClass(  1000.,L)                    # kW
                        microgrid = MicrogridClass()

                        # timestamp this data point
                        results.datetime.append(dt.datetime(2019,1,1) + dt.timedelta(hours=h))

                        #
                        # run the dispatch simulation
                        #
                        if not grid_online:
                            results.time_to_grid_import_h_nf[i] = simulate_resilience(t0,L)
                            # calculate one last result
                            results.onlineTime_h_ni[i] = grid.offlineCounter/4.
                        elif grid_online:
                            simulate_utility_on(t0,L)
                        else:
                            print('\033[1;31;1m fatality: no dispatch type selected')
                            quit()

                    # keep track of asset sizes and pv scaling
                    batt_energy_vals.append(batt_energy)
                    pv_scale_vals.append(pv_scaling_factor)
                    load_scale_vals.append(load_scaling_factor)
                    batt_hrs_vals.append(batt_hrs)
                    batt_power_vals.append(batt_power)
                    gen_power_vals.append(gen_power)

                    # calculate some resilience metrics
                    #   80% @ 1 wk, 50% @ 2 wks
                    ttff  = results.time_to_grid_import_h_nf

                    conf_72h.append(len(ttff[ttff >= 72])/runs)
                    conf_336h.append(len(ttff[ttff >= 336])/runs)

                    max_ttff.append(np.max(ttff))
                    min_ttff.append(np.min(ttff))
                    avg_ttff.append(np.average(ttff))

                    if debug_res:
                        print('TTFF [h]: max={} avg={} min={}'.format(max_ttff,avg_ttff,min_ttff))
                        print('Confidence: 72h={} 336h={}'.format(conf_72h,conf_336h))

                    if output_sim_meta:
                        if superloop_enabled:
                              filename = output_dir + '/sim_meta_{:}_{:.1f}_{:.3f}_{:.0f}_{:.1f}_{:.0f}.csv'.format(site, load_scaling_factor, pv_scaling_factor, batt_power, batt_hrs, gen_power)
                        else:
                            filename = output_dir + '/sim_meta_{}_{}.csv'.format(site, filename_param)
                        with open(filename, 'w', newline='') as file:
                            output = csv.writer(file)
                            output.writerow(['Site',site])
                            output.writerow(['Mamba.py v',__version__])
                            output.writerow(['Datetime',dt.datetime.now()])
                            output.writerow(['Simulated outage duration [days]',days])
                            output.writerow(['Outages simulated',runs])
                            output.writerow([])
                            output.writerow(['PV scaling factor', pv_scaling_factor])
                            output.writerow(['Battery power [kW]',batt_power])
                            output.writerow(['Battery energy [kWh]',batt_energy])
                            output.writerow(['Battery hours [kWh]',batt_hrs])
                            output.writerow(['Generator power [kW]',gen_power])
                            output.writerow(['Generator tank [gal]',gen_tank])
                            output.writerow(['Fuel curve A coefficient [gal/h/kW]',gen_fuelA])
                            output.writerow(['Fuel curve B coefficient [gal/h]',gen_fuelB])
                            output.writerow([])
                            output.writerow(['Program call and args',' '.join(sys.argv)])

                    if output_resilience:
                        if superloop_enabled:
                              filename = output_dir + '/resilience_{:}_{:.1f}_{:.3f}_{:.0f}_{:.1f}_{:.0f}.csv'.format(site, load_scaling_factor, pv_scaling_factor, batt_power, batt_hrs, gen_power)
                        else:
                            filename = output_dir + '/resilience_{}_{}.csv'.format(site, filename_param)
                        with open(filename, 'w', newline='') as file:
                            output = csv.writer(file)
                            output.writerow(['Site',site])
                            output.writerow(['Mamba.py v',__version__])
                            output.writerow(['Datetime',dt.datetime.now()])
                            output.writerow(['Simulated outage duration [days]',days])
                            output.writerow(['Outages simulated',runs])
                            output.writerow([])
                            output.writerow(['PV scaling factor', pv_scaling_factor])
                            output.writerow(['Battery power [kW]',batt_power])
                            output.writerow(['Battery energy [kWh]',batt_energy])
                            output.writerow(['Battery hours [kWh]',batt_hrs])
                            output.writerow(['Generator power [kW]',gen_power])
                            output.writerow(['Generator tank [gal]',gen_tank])
                            output.writerow(['Fuel curve A coefficient [gal/h/kW]',gen_fuelA])
                            output.writerow(['Fuel curve B coefficient [gal/h]',gen_fuelB])
                            output.writerow(['Pull soc0 from previous dispatch?',vary_soc])
                            output.writerow([])
                            output.writerow(['Confidence 72 h',conf_72h])
                            output.writerow(['Confidence 336 h',conf_336h])
                            output.writerow(['Max TTFF [h]', max_ttff])
                            output.writerow(['Avg TTFF [h]', avg_ttff])
                            output.writerow(['Min TTFF [h]', min_ttff])
                            output.writerow([])
                            output.writerow(['Outage','Outage Start', 'Time to First Failure [h]', 'Cumulative Operating Time [h]'])
                            for i in range(runs):
                                output.writerow([i+1,results.datetime[i],results.time_to_grid_import_h_nf[i],results.onlineTime_h_ni[i]])






t_script_finished_dt = dt.datetime.now()
t_elapsed_dt = t_script_finished_dt - t_script_begin_dt
results.code_runtime_s = t_elapsed_dt.total_seconds()



#
# Outputs
#

if superloop_enabled:
    filename = output_dir + '/resilience_superloop.csv'
    with open(filename, 'w', newline='') as file:
        output = csv.writer(file)
        output.writerow(['Site',site])
        output.writerow(['Mamba.py v',__version__])
        output.writerow(['Datetime',dt.datetime.now()])
        output.writerow(['Runtime [s]',results.code_runtime_s])
        if batt_power_varies:
            output.writerow(['Battery power varies such that capacity is ', batt_hrs, 'h'])
        else:
            output.writerow(['Battery power [kW]',batt_power])
        output.writerow(['Gen power [kW]', gen_power])
        output.writerow(['Gen tank size [gal]', gen_tank])
        output.writerow([])
        output.writerow(['Load Scaling factor', 'PV scaling factor','Batt energy [kWh]', 'Batt Power [kW]', 'Batt hrs [h]', 'Generator power [kVA]','Confidence 168h','Confidence 336h','Min TTFF [h]','Avg TTFF [h]', 'Max TTFF [h]'])
        for i in range(len(max_ttff)):
            output.writerow([load_scale_vals[i],pv_scale_vals[i],batt_energy_vals[i], batt_power_vals[i], batt_hrs_vals[i], gen_power_vals[i], conf_72h[i],conf_336h[i],min_ttff[i],avg_ttff[i],max_ttff[i]])

# plots
if plots_on:

    # plot type
    fig, ax1 = plt.subplots(figsize=(20,8.5))
    ax2 = ax1.twinx()

    # main vectors
    t = np.arange(1., L+1., 1.)
    p = pv.P_kw_nf
    l = load.P_kw_nf
    g = gen.P_kw_nf
    d = np.clip(bat.P_kw_nf,0,10000)
    c = -1*np.clip(bat.P_kw_nf,-10000,0)
    u = np.clip(grid.P_kw_nf,0,10000)

    # consider putting in a horizontal demand target line
    if peak_shaving and debug_demand:
        m = load.datetime[i].month
        ax1.plot([1, L], [demand_targets.monthly[m-1], demand_targets.monthly[m-1]], 'r', label='demand target', linewidth=.5)

    # plot load and soc
    ax1.plot(t, l,           'k',    label='load', linewidth=0.5)
    ax1.plot(t, l+c,         'k--',  label='load + batt charge', linewidth=0.5)
    ax2.plot(t, bat.soc_nf,  'b:',   label='soc', linewidth=0.5)

    # area plots
    if plot_pv_first and grid_online:
        ax1.plot(t, u,       'r--',  label='grid import',linewidth=1.)
        ax1.fill_between(t, 0,      p,     facecolor='lightyellow', label='pv')#, interpolate=True)
        ax1.fill_between(t, p,      p+d,   facecolor='lightblue',   label='batt discharge')#, interpolate=True)
        ax1.fill_between(t, p+d,    p+d+u, facecolor='pink',        label='grid')#, interpolate=True)

    # area plots
    elif plot_pv_first and not grid_online:
        ax1.fill_between(t, 0,      p,          facecolor='lightyellow', label='pv')#, interpolate=True)
        ax1.fill_between(t, p,      p+d,        facecolor='lightblue',  label='batt discharge')#, interpolate=True)
        ax1.fill_between(t, p+d,    p+d+g,      facecolor='lightgreen', label='gen')#, interpolate=True)

    # area plots
    elif plot_grid_first and grid_online:
        ax1.fill_between(t, 0,       u,     facecolor='pink', label='grid')#, interpolate=True)
        ax1.fill_between(t, u,       u+p,   facecolor='lightyellow', label='pv')#, interpolate=True)
        ax1.fill_between(t, u+p,     u+p+d, facecolor='lightblue', label='batt discharge')#, interpolate=True)

    ax1.set_xlabel('timestep')
    ax1.set_ylabel('power [kW]')
    ax2.set_ylabel('SOC')

    ax1.legend()
    ax2.legend()
    ax1.set_xlim(1,L)
    ax1.set_ylim(1,1.1*np.max(pv.P_kw_nf))
    ax2.set_ylim(0,1)

    #fig.tight_layout()
    plt.show()

    # if debug_demand:
    #     plt.plot(load.datetime,grid.P_kw_nf)
    #     plt.show()


# print errors
err.print_faults()

if not grid_online and not plots_on:

    if platform.system() == 'Darwin':   # macos
        os.system('say "Ding! (resilient) fries are done"')
    else:
        sys.stdout.write('\a') # beep
