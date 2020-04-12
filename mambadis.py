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


################################################################################
#
# Modules
#
################################################################################

import sys
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
        me.P_kw = np.zeros((length,), dtype=float)
        me.onlineTime = np.zeros((length,),dtype=int)
        me.time_to_grid_import_h = np.zeros((length,), dtype=float) # [h]
        me.code_runtime_s = 0


    def clear(me):
        me.datetime = []
        me.P_kw.fill(0)
        me.time_to_grid_import_h.fill(0)

#
# Generator
#

class GenClass:
    def __init__(me, Pn_kw, a, b, tankCap, tstep, length):
        me.timestep = tstep                             # units of seconds
        me.dpph = 3600./tstep
        me.tstepHr = tstep/3600.
        me.Pn_kw = Pn_kw
        me.fuelcurve_Acoeff = a                         # [gal/h/kw]
        me.fuelcurve_Bcoeff = b                         # [gal/h]
        me.fuelcurve_Ainv = 1/a
        me.fuelTankCapacity = tankCap                   # [gal]
        me.P_kw = np.zeros((length,), dtype=float)
        me.fuelConsumed = 0                              # [gal]

    def clear(me):
        me.P_kw.fill(0)
        me.fuelConsumed = 0


    def power_request(me, i, p_req):
        if p_req > 0:
            p_final = min(p_req,  me.Pn_kw, me.Pmax_tank())
        elif p_req < 0:
            p_final = 0.
        else:
            p_final = 0.
        me.P_kw[i] = np.array([p_final])
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
#
# Grid
#

class GridClass:
    def __init__(me, Pn_kw, length):
        me.Pn_kw = Pn_kw
        me.P_kw = np.zeros((length,), dtype=float)
        me.time_to_import = 0
        me.import_on = 0
        me.offlineCounter = 0

    def clear(me):
        me.P_kw.fill(0)
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
        me.P_kw[index] = p_final
        return p_final

    def timer_tick(me):
        if not me.import_on:
            me.time_to_import += 1

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
        me.P_kw = np.zeros((length,), dtype=float)
        me.soc = np.zeros((length,), dtype=float)
        me.soc_flag = 0                             # binary

    def clear(me):
        me.soc_prev = me.soc0
        me.P_kw.fill(0)
        me.soc.fill(0)
        me.soc_flag = 0

    def power_request(me, index, p_req):
        if p_req > 0:
            p_final = np.amin([p_req, me.Pn_kw, me.P_max_soc()])
        elif p_req < 0:
            p_final = np.amax([p_req, -me.Pn_kw, me.P_min_soc()])
        else:
            p_final = 0.
        me.soc[index] = me.soc_update(p_final)
        me.P_kw[index] = p_final
        return p_final

    def soc_update(me,powerkw):
        soc_new = me.soc_prev - (powerkw * me.timestep/3600.0)/me.En_kwh
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
    if power < 25:
        fuelA = 0.06
        fuelB = 0.3
    elif power < 35:
        fuelA = 0.067
        fuelB = 1.23
    elif power < 50:
        fuelA = 0.8
        fuelB = 1.52
    elif power < 67:
        fuelA = 0.067
        fuelB = 1.73
    elif power < 87:
        fuelA = 0.071
        fuelB = 2.33
    elif power < 120:
        fuelA = 0.064
        fuelB = 2.54
    else:
        fuelA = 0.064
        fuelB = 2.54
        err.gen_fuel_coeffs()
    return [fuelA, fuelB]


#
# Synthetic data creation
#

def create_synthetic_data():
    pv_oneday = 100*np.array([0., 0., 0., 0., 0., 0., 1., 2., 3., 4., 5., 6., 7., 6., 5., 4., 3., 2., 1., 0., 0., 0., 0., 0.])
    load_oneday = 100*np.array([1., 1., 1., 1., 1., 1., 5., 3., 1., 1., 4., 2., 1., 1., 1., 1., 6., 8., 12., 8., 5., 3., 2., 1.])
    pv.P_kw = np.concatenate([pv_oneday,pv_oneday,pv_oneday])
    load.P_kw = np.concatenate([load_oneday,load_oneday,load_oneday])

#
# Import load data
#

def import_load_data(site, load_stats):
    filename = site + 'load.csv'
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
    load_all.P_kw = np.concatenate((md,md),axis=0)

    if load_stats:
        print('max load [kw] = {:.1f}'.format(np.amax(load_all.P_kw)))
        print('avg load [kw] = {:.1f}'.format(np.average(load_all.P_kw)))
        print('min load [kw] = {:.1f}'.format(np.amin(load_all.P_kw)))

#
# Import solar vector
#

def import_pv_data(site):
    filename = site + 'solar.csv'
    with open(filename,'r') as f:
        datacsv = list(csv.reader(f, delimiter=","))
        del datacsv[0]
        t = []
        p = []
        for line in datacsv:
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
    pv_all.P_kw = np.concatenate((pv,pv),axis=0)

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

    load.P_kw = load_all.P_kw[m_0:m_end]
    load.datetime = load_all.datetime[m_0:m_end]

    if debug:
        print('LOAD')
        print('m_0={:d}'.format(m_0) + ' m_end={:d}'.format(m_end))
        print(load.datetime[0])
        print(load.datetime[-1])
        print(load.P_kw.size)
        print(len(load.datetime))
        print('')

    # temporarily store small chunk "all pv" vector in "pv"

    # where in "all PV" data this run will begin
    n_0 = t_0//4 + pv_all.offset      # offset usually 0

    # where in "all load" data this run will end
    n_end = n_0 + L//4                # offset usually 0

    pv.P_kw = pv_all.P_kw[n_0:n_end]
    pv.datetime = pv_all.datetime[n_0:n_end]

    if debug:
        print('PV')
        print('n_0={:d}'.format(n_0) + ' n_end={:d}'.format(n_end))
        print(pv.datetime[0])
        print(pv.datetime[-1])
        print(pv.P_kw.size)
        print(len(pv.datetime))
        print('')

    # check indexing
    # beginning and ending date and hour should match between load and pv
    if (load.datetime[0].day-pv.datetime[0].day):
        err.indexing()
    if (load.datetime[-1].day-pv.datetime[-1].day):
        err.indexing()
    if (load.datetime[0].hour-pv.datetime[0].hour):
        err.indexing()
    if (load.datetime[-1].hour-pv.datetime[-1].hour):
        err.indexing()


#
# Algorithm
#

    for i in range(L):

        # only increment i_pv every 4 i-increments
        i_pv = i//4

        LSimbalance = load.P_kw[i]      -   pv.P_kw[i_pv]   # load-solar imbalance

        battpower = bat.power_request(i,LSimbalance)

        LSBimbalance =  LSimbalance     -   battpower       # load-solar-batt imbalance

        genpower = gen.power_request(i,LSBimbalance)

        LSBGimbalance = LSBimbalance    -   genpower        # etc

        gridpower = grid.power_request(i,LSBGimbalance)

        if gridpower <= 0:
            grid.offlineCounter += 1                        # time that microgrid services load

        # check energy balance
        if np.absolute((LSimbalance - bat.P_kw[i] - gen.P_kw[i] - grid.P_kw[i])) > 0.001:
            err.energy_balance()

    time_to_grid_import = grid.time_to_import/(3600./load.timestep)

    # vectors
    if vectors_on:
        print('{:},'.format('time') + '{:},'.format('load') + '{:},'.format('pv') + '{:},'.format('b_kw') + '{:},'.format('b_soc') + '{:},'.format('gen') + '{:},'.format('grid') + '{:}'.format('diff'))
        for i in range(L):
            l=load.P_kw.item(i)
            p=pv.P_kw.item(i//4)
            b=bat.P_kw.item(i)
            s=bat.soc.item(i)
            g=gen.P_kw.item(i)
            G=grid.P_kw.item(i)
            d=l-p-b-g-G
            print('{:d},'.format(i+1) + '{:.1f},'.format(l) + '{:.1f},'.format(p) + '{:.1f},'.format(b) + '{:.2f},'.format(s) + '{:.1f},'.format(g) + '{:.1f},'.format(G) + '{:.1f}'.format(d))
        print('')

    if debug:
        print('raw sum of batt power [kw]: {:.1f}'.format(np.sum(bat.P_kw).item(0)))

    return time_to_grid_import





################################################################################
#
# "main"
#
################################################################################

t_begin = dt.datetime.now()

err =   FaultClass()

#
# Run options
#

runs = 1#365*8                # number of iterations
skip_ahead = 0               # number of hours to skip ahead
site = 'fish'             # fish, hradult, (hrfire not working)

# physical capacities
batt_power = 25.         # kw
batt_energy = 950.       # kwh
gen_power = 0.           # kw
gen_tank = 250.         # gal

# outputs on/off
output_file_on = 1
vectors_on = 0
plots_on = 0
load_stats = 0
debug = 0

# command line run options
if len(sys.argv) > 1:

    for i in range(1, len(sys.argv)):

        if sys.argv[i] == '-s':
            site = str(sys.argv[i+1])

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

    results.datetime.append(dt.datetime(2019,1,1) + dt.timedelta(hours=h))

    results.time_to_grid_import_h[i] = simulate_outage(t0,L)

    results.onlineTime[i] = grid.offlineCounter/4.

t_end = dt.datetime.now()
t_elapsed = t_end - t_begin
results.code_runtime_s = t_elapsed.total_seconds()



#
# Outputs
#

if output_file_on:
    with open('output.csv', 'w', newline='') as file:
        output = csv.writer(file)
        output.writerow(['Datetime',dt.datetime.now()])
        output.writerow(['Runtime [s]',results.code_runtime_s])
        output.writerow(['Site',site])
        output.writerow(['Runs',runs])
        output.writerow(['Simulation period [days]',days])
        output.writerow(['Battery power [kW]',batt_power])
        output.writerow(['Battery energy [kWh]',batt_energy])
        output.writerow(['Generator power [kW]',gen_power])
        output.writerow(['Generator tank [gal]',gen_tank])
        output.writerow(['Fuel curve A coefficient [gal/h/kW]',gen_fuelA])
        output.writerow(['Fuel curve B coefficient [gal/h]',gen_fuelB])
        output.writerow([])
        output.writerow(['Run','Outage Start', 'Time to First Failure [h]', 'Cumulative Operating Time [h]'])
        for i in range(runs):
            output.writerow([i+1,results.datetime[i],results.time_to_grid_import_h[i],results.onlineTime[i]])

# plots
if plots_on:
    t = np.arange(1., L+1., 1.)

    fig, ax1 = plt.subplots(figsize=(20,10))


    bat_p_kw_pos = np.clip(bat.P_kw,0,10000)
    bat_p_kw_neg = -1*np.clip(bat.P_kw,-10000,0)

    ax2 = ax1.twinx()
    ax1.plot(t,pv.P_kw,                             'y',    label='PV')
    ax1.plot(t,pv.P_kw + bat_p_kw_pos,              'b',    label='battdis+PV')
    ax1.plot(t,pv.P_kw + bat_p_kw_pos + gen.P_kw,   'g',    label='gen+battdis+PV')
    ax1.plot(t,load.P_kw,                           'k',    label='load')
    ax1.plot(t,load.P_kw + bat_p_kw_neg,            'k-*',  label='load+battchg')
    ax2.plot(t,bat.soc,                             'm-',   label='soc', marker='')

    ax1.fill_between(t, 0, pv.P_kw, facecolor='yellow')#, interpolate=True)
    ax1.fill_between(t, pv.P_kw, pv.P_kw+bat_p_kw_pos, facecolor='blue')#, interpolate=True)
    ax1.fill_between(t, pv.P_kw+bat_p_kw_pos, pv.P_kw+bat_p_kw_pos+gen.P_kw, facecolor='green')#, interpolate=True)
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
