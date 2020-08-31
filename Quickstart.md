# Quickstart

mamba.py v6.1 | Michael Wood | 2020.8.31




## Overview
The 'mamba dispatch' software uses simplified dispatch strategies for very fast simulation of utility connected or islanded microgrids. Peak shaving and arbitrage are the two main strategies. Mamba is written and tested in Python 3.7.

1. **Single Simulation, or 'run'.**
  - For utility connected simulation we simply calculate the residual load (actual load less PV production) and make a decision to dispatch the battery based on the selected strategy - any available generator is probably not used. This simulation can be anywhere from a few timesteps to a year long.
  - For resilience simulations we initiate a grid outage at time t, and then simulate the PV-battery-generator microgrid for a reasonable time period, say two weeks. At the end we are given basic information about how long until the microgrid failed to serve the load (Time to First Failure) and of the whole time period how much time the load was served.
2. **Year of Resilience Simulations.** A single outage is interesting, but we really want to know how our microgrid will perform at any time of the year. So we simulate one outage at Jan 1 0:00, then another at  Jan 1 3:00, then Jan 1 6:00, all the way to Dec 31 23:00. And by default we jump in 3 hr increments. A year takes about one minute to simulate. This does not apply to utility connected cases.
3. **Matrix of Simulated Years, or 'superloop'.** The year of outages was only for a single given battery size, PV size, and generator size (could be 0). If you want to know the results for an NxM matrix of battery and PV sizes, run a superloop. In this version a generator axis is not implemented.



## Setup

### Installing Python

The internet is not-literally full of great ways to get python 3 on your machine, but if you have any doubt I recommend the ~~PyCharm~~ Atom editor and Anaconda (often just *conda*) library manager. However pip works well too.

### Dependencies

The only non-standard library ('package') should be numpy, which Anaconda can help install. But on a clean Python 3 install the pip could be quicker.

`pip install numpy`

All the commands in this doc are from a normal terminal shell.

### Input Data

Input load and solar data should be 15-minute interval, but if that’s not available we can flip a switch. In general site names should match the load and solar datafiles in the I-hope-obvious way.



## Running

### Inputs

Besides input data, the most important considerations are choosing the simulation type and defining the simulation parameters. You could always edit these in the script directly. But giving command line arguments should be faster, clearer, and avoid mistakes.

* Arguments, or 'flags' are explained if you search on “help” in the script or do:

  `python mambadis.py --help`

* Basic example argument list. In order: simulation type, site, battery power, battery energy, generator power, generator tank size, generator fuel is propane.

  `python mambadis.py -sim r -s fish -bp 30 -be 60 -gp 30 -gt 120 -gfp`

* For sites without a generator, ~~set~~ the generator power **is set by default** to 0. For no battery set the battery power to 0.

* Battery **power** defaults to the same numerical value (different units) as the battery **energy**, which is just to say that we always simulate a '1 h' battery. This can be overridden with the -bp [kW] argument.

* To run a superloop you will need to edit the code to define the vector of battery sizes and PV scaling factors. Search on 'if superloop enabled' and remember to include the superloop flag for a successful execution:

  `python mambadis.py -sim r -s fish -sl`

* When running a superloop, consider first testing your matrix of pv-battery sizes with a single run per Simulated 'Year'. If you messed up, you'll know a lot sooner. Example:

  `python mambadis.py -sim r -s fish -sl -r 1`

### Errors

You may see an error output about mismatched indices, this is a problem but I don't think it affects the final numbers in a big way. I have it on a to-fix list.



## Analysis


### Outputs
Harkening back to the three nested loops, we have:

1. **Single Simulation, or 'run.'** The file vectors.csv (off by default) gives you all the dispatch info for a single two-week outage. Warning: only turn this on for a single simulation, or the code will run a lot slower. This is because normally we simulate 2920 (a whole year, in 3 hr increments) at a time, so we throw out the dispatch data for speed. Try:

   `	python mambadis.py -s fish -be 60 -r 1 -v`

With the above arguments you'll get only the dispatch starting at 12am Jan. To skip ahead and get any 14-day window of the year do (e.g. 24 hrs):

​		`	python3 mambadis.py -s fish -be 60 -r 1 -v -sk 24`


2. **Year of Resilience Simulations (default).** Always produced by default, output.csv contains a year of time to first failure and cumulative operating time data for your chosen pv-battery-generator size. Each row of data is for one Single Outage.

3. **Matrix of Simulated Years, or 'superloop.'** Each superloop.csv row represents a Simulated Year of outages, summarized with confidence and TTFF statistics. Clearly this file is only output with the -sl flag.

### More Analysis

I included an Analysis folder which has some dumb spreadsheets for Apple Numbers, if you don’t have that I can export to excel but it’s a little messy. This was a time-crunch kludge — clearly the code should produce all the plots.
