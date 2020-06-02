# mambadis.py
# mjw 2020.3.26
# muGrid Analytics
# python 3.7

# Fast, cycle charging dispatch of PV-battery-generator migrogrids

#
# Versions
#

# 0.0 - first attempt with synthetic data (no rev number in file)
#   0.1 - load in real PV and load data
# 1.0 - separate synthetic and real data
#   1.1 - match up load data to PV vector, PV data in kW
#   1.2 - "all" load/pv data as well as L-length data, tested on different t0's
#   1.3 - no more "if LSimbalance>0" logic (not needed)
# 2.0 - "test_vX.py" changed to "lilcc_vX.py"
#   2.1 - include gen fuel calc and tank size (not working now)
#   2.2 - simulate_outage() works for many loops but "indexing" error, can't pre-allocate load. pv. bat. etc
#   2.3 - switch to python 3.7 for dev, note that input data tested with HR Fire
# 3.0 - manual changes to work with FISH data, PV timeseries wrap around 12/31
#   3.1 - ** SOC_0 = 1.0 ***, load stats, tested on HR Adult Center (manual changes)
#   3.2 - now actually tested on python 3, pass sitename to load function, HR fire broken
#   3.3 - rename to mambadis.py ("fast snake" dispatch)
#   3.4 - ** run from command line **, output to csv filename
#   3.5 - measure runtime, command line run options overwrite defaults
#   3.6 - simplify offsets, fuel curve lookup table, some small things
#   3.7 - add _nf to numpy float and _ni to numpy int array variable names, _dt to datetime variables, remove some implicit casts
#   3.8 - can accept 15 min interval PV data (not fully tested)
#   3.9 - ** PV 15 min interval data tested **, new filename
#   3.10 - simulation dispatch vector now output to csv
#   3.11 - fix bug where simulation dispatch vector showed wrong P_pv
# 4.0 - real CC dispatch including min on/off time, nonzero batt efficiency, need to clean up algo & "grid" aspect
#   4.1 - change dispatch when gen tank is empty
#   4.2 - min on/off time was causing early microgrid failures
#   4.3 - help
#   4.4 - help is broken (removed), add pv scaling factor 
#   4.5 - fixed help, -sk skip ahead arg, new folder organization, new filenames for clarity
#   4.6 - keep mambadis.py in main dir, output resilience conf and ttff
# 5.0 - superloop runs a matrix of battery energy and PV scale (oversize) factors, output and superloop files 
#   5.1 - bug where superloop output fails
#   5.2 - for superloop batt power varies to always be 1h cap
#   5.3 - varying batt power doesn't do much, leave out for now
#   5.4 - bug fix: for super small load-gen imbalance code would declare failure, microgrid class

################################################################################
#
# Modules
#
################################################################################

import sys
import numpy as np
#import matplotlib.pyplot as plt
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
        me.rt_eff = 0.9

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
        me.soc_nf[index] = me.soc_update(p_final)
        me.P_kw_nf[index] = p_final
        return p_final

    def soc_update(me,powerkw):
        soc_new = me.soc_prev - (powerkw * me.timestep/3600.0 * me.rt_eff)/me.En_kwh
        me.soc_prev = soc_new
        return soc_new

    def P_max_soc(me):
        return me.soc_prev * me.En_kwh * (3600.0/me.timestep)

    def P_min_soc(me):
        return (me.soc_prev - 1.) * me.En_kwh * (3600.0/me.timestep)
#
# Faults
#

class FaultClass:
    def __init__(me):
        me.mainloop = 0
        me.energybalance = 0
        me.index = 0
        me.fuelCurveCoeffs = 0

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

    def main_loop(me):
        me.mainloop += 1

    def energy_balance(me):
        me.energybalance += 1

    def indexing(me):
        me.index += 1

    def gen_fuel_coeffs(me):
        me.fuelCurveCoeffs = 1


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

def lookup_fuel_curve_coeffs(power):
    # coeffs = [fuel_curve_coeff_A, fuel_curve_coeff_B]
    if power < 25:      coeffs = [0.06 , 0.3]
    elif power < 35:    coeffs = [0.067, 1.23]
    elif power < 50:    coeffs = [0.8, 1.52]
    elif power < 67:    coeffs = [0.067, 1.73]
    elif power < 87:    coeffs = [0.071, 2.33]
    elif power < 120:   coeffs = [0.064, 2.54]
    else:
        coeffs = [0.064, 2.54]
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
    md = my_data[1:,1]
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
        filename = '../Data/Solar/' + site + '_solar.csv'

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
# Simulate outage
#

def simulate_outage(t_0,L):

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

    if debug:
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

    if debug:
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

            if bat.soc_prev < 0.001:
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

            if (bat.soc_prev == 1) or ((bat.soc_prev > 0.5) and (pv.P_kw_nf[i_pv] > 0)):
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
    if vectors_on:
        with open('./Data/Output/vectors.csv', 'w', newline='') as file:
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
                output.writerow([i+1,l,p,b,s,g,G,d])

    if debug:
        print('checksum: {:.1f}'.format(np.sum(bat.P_kw_nf)))

    return time_to_grid_import

#
# Print help 
#

def help_printout():
    print('\nmambadis.py\nmuGrid Analytics LLC\nMichael Wood\nmichael.wood@mugrid.com')
    print('')
    print('Arguments can be included in any order, but values must follow directly after keys')
    print('e.g. % python mambadis.py -s fish -bp 20 be 40 .. (etc)')
    print('')
    print('Required Command Line Arguments')
    print(' Site name:              -s [sitename]       e.g. -s hradult')
    print(' Battery power:          -bp [power kW]      e.g. -bp 60')
    print(' Battery energy:         -be [energy kWh]    e.g. -be 120') 
    print(' Generator power:        -gp [power kW]      e.g. -gp 60')
    print(' Generator tank:         -gt [size gal]      e.g. -gt 200')
    print('')
    print('Optional Command Line Arguments')
    print(' Run "n" simulations:     -r [n]             e.g. -r 1   default=2920')
    print(' Dispatch vectors ON:     -v                 e.g. -v     default=OFF')
    print(' Load stats ON:           -l                 e.g. -l     default=OFF')
    print(' Debug ON:                -d                 e.g. -d     default=OFF')
    print(' Skip ahead "h" hours:    -sk [h]            e.g. -sk 24 default=OFF')
    print(' Superloop enable:        -sl                e.g. -sl    default=DISABLED')
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

runs = 365*8                  # number of iterations
skip_ahead = 0                  # number of hours to skip ahead
site = 'hradult'                   # fish, hradult, (hrfire not working)
solar_data_inverval_15min = 1
superloop_enabled = 0           # run matrix of battery energy and PV scale (oversize) factors

# physical capacities
batt_power = 60.         # kw
batt_energy = 60.       # kwh
gen_power = 60.           # kw
gen_tank = 250.          # gal

# outputs on/off
superloop_file_on = 0
output_file_on = 1
vectors_on = 0
plots_on = 0
load_stats = 0
debug = 0

# pv scaling
pv_scaling_factor = 1

# command line run options
if len(sys.argv) > 1:

    for i in range(1, len(sys.argv)):

        if sys.argv[i] == '-s':
            site = str(sys.argv[i+1])

        elif sys.argv[i] == '-p':
            pv_scaling_factor = float(sys.argv[i+1])

        elif sys.argv[i] == '-r':
            runs = int(sys.argv[i+1])

        elif sys.argv[i] == '-bp':
            batt_power = float(sys.argv[i+1])

        elif sys.argv[i] == '-be':
            batt_energy = float(sys.argv[i+1])

        elif sys.argv[i] == '-gp':
            gen_power = float(sys.argv[i+1])

        elif sys.argv[i] == '-gt':
            gen_tank = float(sys.argv[i+1])

        elif sys.argv[i] == '-v':
            vectors_on = 1

        elif sys.argv[i] == '-l':
            load_stats = 1

        elif sys.argv[i] == '-d':
            debug = 1

        elif sys.argv[i] == '-sk':
            skip_ahead = int(sys.argv[i+1])
        
        elif sys.argv[i] == '-sl':
            superloop_enabled = 1

        elif sys.argv[i] == '--help' :
            help_printout()
            quit()

if superloop_enabled:
    pv_scale_vector = [1,1.5,2,2.5,3,3.5]
    batt_energy_vector = [30,60,125,250]
else:
    pv_scale_vector = [pv_scaling_factor]
    batt_energy_vector = [batt_energy]


pv_scale_vals= []
batt_energy_vals = []

max_ttff = []
avg_ttff = []
min_ttff = []

conf_168h = []
conf_336h = []

for pv_scaling_factor in pv_scale_vector:
    for batt_energy in batt_energy_vector: 
        [gen_fuelA, gen_fuelB] = lookup_fuel_curve_coeffs(gen_power)

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
        import_load_data(site, load_stats)
        import_pv_data(site)

        #
        # some data prep
        #

        # window start and size
        days = 14
        L = days*24*4                     # length of simulation in timesteps

        #
        # Simulation Loop
        #

        for i in range(runs):
            h=3*i + skip_ahead
            t0=h*4       # where to start simulation in "all load" vector
            
            # start with fresh variables
            # wish we could pre-allocate these, but it was causing a bug
            #   even after calling .clear() on everything
            load =  DataClass(  15.*60.,L)                 # timestep[s]
            pv =    DataClass(  60.*60.,L)                   # timestep[s]
            gen =   GenClass(   gen_power,gen_fuelA,gen_fuelB,gen_tank,15.*60.,L)   # kW, fuel A, fuel B, tank[gal], tstep[s]
            bat =   BattClass(  batt_power,batt_energy,1.0,15*60.,L)      # kW, kWh, soc0 tstep[s]
            grid =  GridClass(  100.,L)                    # kW
            microgrid = MicrogridClass()

            # timestamp this data point
            results.datetime.append(dt.datetime(2019,1,1) + dt.timedelta(hours=h))
            
            #
            # run the dispatch simulation
            #
            results.time_to_grid_import_h_nf[i] = simulate_outage(t0,L)

            # calculate one last result
            results.onlineTime_h_ni[i] = grid.offlineCounter/4.

        # keep track of battery power and pv scaling
        batt_energy_vals.append(batt_energy)
        pv_scale_vals.append(pv_scaling_factor)

        # calculate some resilience metrics
        #   80% @ 1 wk, 50% @ 2 wks
        ttff  = results.time_to_grid_import_h_nf

        conf_168h.append(len(ttff[ttff >= 168])/runs) 
        conf_336h.append(len(ttff[ttff >= 336])/runs) 

        max_ttff.append(np.max(ttff))
        min_ttff.append(np.min(ttff))
        avg_ttff.append(np.average(ttff))

        print(str(pv_scaling_factor) +' ' + str(batt_power) + ' ' + str(batt_energy) + ' ' + str(len(ttff[ttff>=168])/runs))

t_script_finished_dt = dt.datetime.now()
t_elapsed_dt = t_script_finished_dt - t_script_begin_dt
results.code_runtime_s = t_elapsed_dt.total_seconds()


#
# Outputs
#

if superloop_enabled:
    with open('./Data/Output/superloop.csv', 'w', newline='') as file:
        output = csv.writer(file)
        output.writerow(['Datetime',dt.datetime.now()])
        output.writerow(['Runtime [s]',results.code_runtime_s])
        output.writerow(['Site',site])
        output.writerow(['Batt power [kW]',batt_power])
        output.writerow(['Gen power [kW]', gen_power])
        output.writerow(['Gen tank size [gal]', gen_tank])
        output.writerow([])
        output.writerow(['PV scaling factor','Batt energy [kWh]','Confidence 168h','Confidence 336h','Min TTFF [h]','Avg TTFF [h]', 'Max TTFF [h]'])
        for i in range(len(max_ttff)):
            output.writerow([pv_scale_vals[i],batt_energy_vals[i],conf_168h[i],conf_336h[i],min_ttff[i],avg_ttff[i],max_ttff[i]])

if output_file_on:
    with open('./Data/Output/output.csv', 'w', newline='') as file:
        output = csv.writer(file)
        output.writerow(['Site',site])
        output.writerow(['Datetime',dt.datetime.now()])
        output.writerow(['Runtime [s]',results.code_runtime_s])
        output.writerow(['Simulated outage duration [days]',days])
        output.writerow(['Outages simulated',runs])
        output.writerow([])
        output.writerow(['PV scaling factor', pv_scaling_factor])
        output.writerow(['Battery power [kW]',batt_power])
        output.writerow(['Battery energy [kWh]',batt_energy])
        output.writerow(['Generator power [kW]',gen_power])
        output.writerow(['Generator tank [gal]',gen_tank])
        output.writerow(['Fuel curve A coefficient [gal/h/kW]',gen_fuelA])
        output.writerow(['Fuel curve B coefficient [gal/h]',gen_fuelB])
        output.writerow([])
        output.writerow(['Confidence 168 h',conf_168h])
        output.writerow(['Confidence 336 h',conf_336h])
        output.writerow(['Max TTFF [h]', max_ttff])
        output.writerow(['Avg TTFF [h]', avg_ttff])
        output.writerow(['Min TTFF [h]', min_ttff])
        output.writerow([])
        output.writerow(['Outage','Outage Start', 'Time to First Failure [h]', 'Cumulative Operating Time [h]'])
        for i in range(runs):
            output.writerow([i+1,results.datetime[i],results.time_to_grid_import_h_nf[i],results.onlineTime_h_ni[i]])

# plots
if plots_on:
    t_ni = np.arange(1., L+1., 1.)

    fig, ax1 = plt.subplots(figsize=(20,10))


    bat_p_kw_pos_nf = np.clip(bat.P_kw_nf,0,10000)
    bat_p_kw_neg_nf = -1*np.clip(bat.P_kw_nf,-10000,0)

    ax2 = ax1.twinx()
    ax1.plot(t_ni, pv.P_kw_nf,                             'y',    label='PV')
    ax1.plot(t_ni, pv.P_kw_nf + bat_p_kw_pos_nf,              'b',    label='battdis+PV')              # implicit cast? x3
    ax1.plot(t_ni, pv.P_kw_nf + bat_p_kw_pos_nf + gen.P_kw_nf,   'g',    label='gen+battdis+PV')
    ax1.plot(t_ni, load.P_kw_nf,                           'k',    label='load')
    ax1.plot(t_ni, load.P_kw_nf + bat_p_kw_neg_nf,            'k-*',  label='load+battchg')
    ax2.plot(t_ni, bat.soc_nf,                             'm-',   label='soc', marker='')

    ax1.fill_between(t_ni, 0, pv.P_kw_nf, facecolor='yellow')#, interpolate=True)
    ax1.fill_between(t_ni, pv.P_kw_nf, pv.P_kw_nf+bat_p_kw_pos_nf, facecolor='blue')#, interpolate=True)
    ax1.fill_between(t_ni, pv.P_kw_nf+bat_p_kw_pos_nf, pv.P_kw_nf + bat_p_kw_pos_nf + gen.P_kw_nf, facecolor='green')#, interpolate=True)
    ax1.set_ylabel('between y1 and 0')

    ax1.set_xlabel('h')
    ax1.set_ylabel('kW')
    ax2.set_ylabel('SOC')

    #ax1.legend(loc='upper left')
    ax2.legend(loc='lower left')
    ax1.set_xlim(1,L)

    #fig.tight_layout()
    plt.show()


# errors
err.print_faults()

# ding fries are done
print('\a')
