import numpy as np

# path to save file
drk = 'data/intfer/'

# Average time of observation
Nd = 0.02083333                          # in day (30 minutes)
dt = 1800                                # in seconds (30 minutes)

# number of steps for integration
nstep = 3770                             # so that one orbit can be covered, as period is 78.53 days

# Number of observation night
Night = 17

# starting days of observation for SII (from April 26, 2023 to May 12, 2023)
start = ['2460061.23264', '2460062.23264', '2460063.23264', '2460064.23264', '2460065.23264', '2460066.23264', '2460067.23264', '2460068.23264', '2460069.23264', '2460070.23264', '2460071.23264', '2460072.24097', '2460073.27986', '2460074.32292', '2460075.36875', '2460076.41458', '2460077.45972']

# ending days of observation for SII (from April 26, 2023 to May 12, 2023)
end = ['2460061.40208', '2460062.44028', '2460063.47708', '2460064.51389', '2460065.55000', '2460066.58611', '2460067.62222', '2460068.65972', '2460069.67292', '2460070.67292', '2460071.67292', '2460072.67292', '2460073.67292', '2460074.67292', '2460075.67292', '2460076.67292', '2460077.67292']

# number of steps each day of observatin according to average time dt
def steps():
    steps = []
    for i in range(0, Night, 1):
        dobs = float(end[i]) - float(start[i])
        dobs /= Nd
        steps.append(int(dobs + 1 ))
    return np.array(steps,dtype=str)

step = list(steps())

# open a txt file to save the starting time, ending time and steps to be taken for observation
if __name__ == "__main__":   
   f = open(drk + "julian.txt", "w")

   for i in start:
       f.write(i + str('\t\t'))
   f.write('\n')

   for i in end:
       f.write(i + str('\t\t'))
   f.write('\n')

   for i in step:
       f.write(i + str('\t\t\t'))
   f.write('\n')
   f.close()

# number of starting, observing and ending step (in nstep) each day. 
#This will help to ignore the unnecessary step of observation because of day time and gamma observation time.
ds = []
for i in range(0, Night, 1):
    a = float(start[i]) - float(start[0])
    a /= Nd
    ds.append(int(a))

dst = np.array(ds)                      # starting steps of each day
step = np.array(step,dtype=int)         # steps to take
den = dst + step                        # ending steps after step

