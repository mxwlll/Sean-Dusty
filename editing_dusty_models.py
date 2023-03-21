# -*- coding: utf-8 -*-
"""
Created on Thu May  5 20:35:58 2022

@author: seand
"""


import os

#This directs your code to your Model Library.
library_path = "C:/UTSA/Research/DUSTY/DUSTY/Sean_mod/Model_library/"

# Update this for every new model you make. It will create a new folder in your 
# library with a 4-digit model number. 

directory = '9999'

pathy = os.path.join(library_path, directory)


if os.path.exists(pathy)==False:
    os.mkdir(pathy)



'''
Here is where things get messy. I wrote this to open up a former model I made, 
copy it into a new file, and then edit the new file. You'll need to change the 
prev_mod into whatever model number is your base file. 

'''
prev_mod = '0000'

prev_file = open(library_path+prev_mod+'/sphere-1_{0}.inp'.format(prev_mod), 'r')

new_inp = open(pathy+'/sphere-1_{0}.inp'.format(directory), 'w')

for line in prev_file:
    new_inp.write(line)

'''
Here's the easy part: here's where you can adjust individual parameters in your 
dusty input! 
'''


temperature = 3000
Td = 1200
abundance = 0.002, 0.013, 0.50, 0.03, 0.20,  0.007, 0.0
shell_thickness = 1000.0
op_prop_index = 3       #used for C_abs
power = 2.5
tau_min = 0.50
tau_max = 1.0
num_mods = int(3)
sil_ow = 0.0
sil_oc = 0.0
sil_dl = 0.0
grf_dl = 0.0
amc_hn = 0.0
sic_pg = 0.0
num_pwrs = 1
comp1 = 'Corundum'
comp2 = 'Mg2Fe8O'
comp3 = 'FeHofm'
comp4 = 'Al2O3-comp'
comp5 = 'cosmic_glass_20120815_DUSTY'
comp6 = 'Mg1Fe9O'
comp7 = 'FeO'
num_comp = '7'
spec_shape = 5
sio_dep = 10
f_lam_dat = 'phnx0'
'''
spec_shape = 1 ----> Blackbody
spec_shape = 2 ----> Engelke-Marengo 
spec_shape = 3 ----> Broken Power law
spec_shape = 4 ----> lamda_F_lamda file
spec_shape = 5 ----> F_lamda (i.e. kurucz10.dat, Mstar, Phnx3000)
spec_shape = 6 ----> F_nu
'''


##############################################################################


with open(pathy+'/sphere-1_{0}.inp'.format(directory), 'r') as file:
    data = file.readlines()
'''
This might mess some things up. You'll need to adjust the first string after 
the parenthese after replace with whatever the base file is. 
'''
temp_loc = data[44]
data[44] = temp_loc.replace('2000', str(temperature))

td_loc = data[47]
data[47] = td_loc.replace('1000', str(Td))

sil_loc = data[56]








sil_list = sil_loc.split()
sil_list[0] = '     x'
if op_prop_index == 3:
    sil_list[1] = ' ' 
sil_list[2] = str(sil_ow)
sil_list[3] = str(sil_oc)
sil_list[4] = str(sil_dl)
sil_list[5] = str(grf_dl)
sil_list[6] = str(amc_hn)
sil_list[7] = str(sic_pg)+'\n'


data[56] = '     '.join(sil_list)





# data[56] = sil_loc.replace('0.50', str(sil_ow))

num_comp_loc = data[57]
data[57] = num_comp_loc.replace('4', num_comp)

comp1_loc = data[58]
data[58] = comp1_loc.replace('Corundum', comp1)

comp2_loc = data[59]
data[59] = comp2_loc.replace('FeO', comp2)

comp3_loc = data[60]
data[60] = comp3_loc.replace('FeHofm', comp3)

comp4_loc = data[61]
data[61] = comp4_loc.replace('Al2O3-comp', comp4)

data.insert(62, '	Lib_nk\{}.nk\n'.format(comp5))

abun_loc = data[63]
data[63] = abun_loc.replace('0.005,0.1,0.50,0.13', str(abundance))

num_pwr_loc = data[72]
data[72] = num_pwr_loc.replace('1', str(num_pwrs))

shell_thic_loc = data[73]
data[73] = shell_thic_loc.replace('50.0', str(shell_thickness))

power_loc = data[74]
data[74] = power_loc.replace('2.7', str(power))

tau_loc1 = data[81]
data[81] = tau_loc1.replace('1.1', str(tau_min))

tau_loc2 = data[81]
data[81] = tau_loc2.replace('2.0', str(tau_max))

num_mod_loc = data[82]
data[82] = num_mod_loc.replace('3', str(num_mods))

#anything beyond the first 5 components can get posted here
#abundences can be posted at the end of the tuple above.
data.insert(63, '	Lib_nk\{}.nk\n'.format(comp6))
data.insert(63, '	Lib_nk\{}.nk\n'.format(comp7))


spec_shapeloc = data[42]
data[42] = spec_shapeloc.replace('1', str(spec_shape))
bb_loc = data[43]
bb_loc1 = data[44]

if int(spec_shape)==5:
    data[43] = bb_loc.replace('=', ' ')
    data[44] = bb_loc1.replace('=', ' ')
    data.insert(43, '                 {}.dat\n'.format(f_lam_dat))

if int(spec_shape)==2:
    data[43] = bb_loc.replace('=', ' ')
    data.insert(45, '                 SiO absorption depth = {} percents\n'.format(sio_dep))





with open(pathy+'/sphere-1_{0}.inp'.format(directory), 'w') as file:
    file.writelines(data)



# This is the string that gets updated to our dusty.inp file. The dusty_string
# is what gets added to dusty.inp, and dusty_dir is where dusty.inp file is  
dusty_string = 'Sean_mod/Model_Library/{0}/sphere-1_{0}\n'.format(directory)

dusty_dir = "C:/UTSA/Research/DUSTY/DUSTY/dusty.inp"

with open(dusty_dir, 'r') as dusty:
    data2 = dusty.readlines()
    data2[9] = dusty_string
with open(dusty_dir, 'w') as dusty:
    dusty.writelines(data2)



# os.system('python pandas_excel_updater.py')




