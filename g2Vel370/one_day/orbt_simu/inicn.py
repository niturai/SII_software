par =  "#M (Msun),      P (days),         a (au),         e,            I (deg),       Omega (deg),        epoch (days),          omega (deg)"


obj1 = ['28.5']                                                                                      # parameters for Gamma Velorum 2O

obj2 = ['9.0', '78.53', '1.2', '0.334', '65.5', '247.7', '2450120.9', '67.4']                        # parameters for Gamma Velorum 2W

obj3 = ['14', '157639398.94', '13791', '0.326', '65', '120', '2450120.9', '6']                                  # parameters for Gamma Velorum 1

f = open("data/gamma_ini.txt", "w")


f.write(par)
f.write('\n')

for i in obj1:
    f.write(i + str('\t\t'))

f.write('\n')

for i in obj2:
    f.write(i + str('\t\t'))
   
f.write('\n')

for i in obj3:
    f.write(i + str('\t\t'))

f.close()
