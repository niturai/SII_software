par =  "#M (Msun),        P (yr),               a (au),         e,            I (deg),       Omega (deg),     epoch (yr),          omega (deg)"


obj1 = ['28.5']                                                                                      # parameters for Gamma Velorum 2O

obj2 = ['9.0', '0.215', '1.2', '0.326', '65', '232.7', '2450120.5', '68']                      # parameters for Gamma Velorum 2W

obj3 = ['14', '271394', '15537', '0.326', '65', '120', '24e9', '6']                      # parameters for Gamma Velorum 1

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
