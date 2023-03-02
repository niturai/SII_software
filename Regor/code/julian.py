import numpy as np

# Average time of observation
Nd = 0.02083333                          # in day
dt = 1800                                # in seconds (30 seconds)


# number of steps for integration of one orbit 
nstep = 3770                          # number of time step, so that one orbit can cover as period is 78.53 days


# Number of observation night (observation length)
Night = 17

# starting days of observation for SII
start = ['2460061.23264', '2460062.23264', '2460063.23264', '2460064.23264', '2460065.23264', '2460066.23264', '2460067.23264', '2460068.23264', '2460069.23264', '2460070.23264', '2460071.23264', '2460072.24097', '2460073.27986', '2460074.32292', '2460075.36875', '2460076.41458', '2460077.45972']

# ending time of observation for SII
end = ['2460061.40208', '2460062.44028', '2460063.47708', '2460064.51389', '2460065.55000', '2460066.58611', '2460067.62222', '2460068.65972', '2460069.67292', '2460070.67292', '2460071.67292', '2460072.67292', '2460073.67292', '2460074.67292', '2460075.67292', '2460076.67292', '2460077.67292']


# number of steps for average time Nd
def steps():
    steps = []
    for i in range(0, Night, 1):
        dobs = float(end[i]) - float(start[i])
        dobs /= Nd
        steps.append(int(dobs + 1 ))
    return np.array(steps,dtype=str)

step = list(steps())

if __name__ == "__main__":
   # open a txt file to save the julian day time of observation
   f = open("data/julian.txt", "w")

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


# number of steps of orbit to ignore because of day time and gamma observation. 
#It is usefull during the plot of orbit on the observation time
ds = []
for i in range(0, Night, 1):
    a = float(start[i]) - float(start[0])
    a /= Nd
    ds.append(int(a))

dst = np.array(ds)                      # starting steps in nstep
step = np.array(step,dtype=int)         # steps to take
den = dst + step                        # ending steps after step


