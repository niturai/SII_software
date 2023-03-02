## parameters which are needed for orbit integration
par =  "M (Msun),      P (days),        a (au),         e,            I (deg),       Omega (deg),        epoch (days),          omega (deg)"



# number of objects in the system 
obj1 = ['28.5']                                                                      # Blue Supergiant mass (the origin point)
obj2 = ['9.0', '78.53', '1.2', '0.334', '65.5', '247.7', '2450120.9', '67.4']        # Wolf-Rayet parameters (with respect to Blue Supergiant)
obj3 = ['14', '1.6e8', '13791', '0.326', '65', '120', '2450120.9', '6']       # Gamma 1 Velorum situated far away from (BS+WR)


# open a txt file to save the parameter
f = open("data/para.txt", "w")                         # w stands for writing

# write the parameter
f.write(par)
f.write('\n')                                          # new line                                       

# write the objects
for i in obj1:
    f.write(i + str('\t\t\t'))
f.write('\n')

for i in obj2:
    f.write(i + str('\t\t'))
f.write('\n')

for i in obj3:
    f.write(i + str('\t\t'))

# close the txt file
f.close()
