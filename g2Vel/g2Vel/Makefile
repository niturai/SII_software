# create directories to save files
folder:
	@mkdir -p data/orbit data/intfer/465                                            # orbit data and signal for observation wavelength
	@mkdir -p fig/orbit/pos fig/orbit/orb fig/intfer/465/sig fig/intfer/465/cntr    # fig for position of objects and for signal
	@mkdir -p gif/orbit gif/sig465                                                  # gif movie for orbit and signal on observational day

# the orbit simulation
orbit:
	@python3 parameter.py                                                            # mass, period and 6 orbital component
	@python3 initial.py                                                              # initial position and velocity of objects
	@python3 orbit.py                                                                # integrated orbit for one period
	
# intferometric simulation
intfer:
	@python3 julian.py                                                               # obbservational time for baselines
	@python3 data.py                                                                 # simulated data for observational wavelength
	@sbatch simu.sh                                                                  # interferometric data

# plot all the data also gif movie	
plot:
	@python3 orbitplot.py                                                            # plot the orbit
	@python3 dataplot.py                                                             # Interferometric signal 
	@python3 paraplot.py                                                             # estimated parameters with corner
	@python3 makegif.py                                                              # movies of data for 17 day of observation

# cleaning	
clean:
	@rm -rf __pycache__
