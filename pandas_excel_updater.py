# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 14:40:02 2022

@author: seand
"""


import pandas as pd
from editing_dusty_models import (directory, temperature, Td, power, tau_min,
tau_max, shell_thickness, abundance, comp1, comp2, comp3, comp4, comp5, comp6,
comp7, data, sil_ow, sil_oc, sil_dl, grf_dl, amc_hn, sic_pg, spec_shape, f_lam_dat,
sio_dep)




excel_sheet = pd.read_excel('./Model Parameters.xlsx')
rowname = int(directory) - 5

excel_sheet.iloc[rowname, 0] = int(directory)



if spec_shape==5:
    excel_sheet.iloc[rowname, 1] = f_lam_dat
elif spec_shape==2:
    excel_sheet.iloc[rowname, 1] = 'EM '+ str(temperature)+ ' SiO % = ' + str(sio_dep)
elif spec_shape==3:
    excel_sheet.iloc[rowname, 1] = 'brok-Pwr-law'
else:
    excel_sheet.loc[rowname, 1] = temperature
        

excel_sheet.iloc[rowname, 2] = Td
if type(power)==tuple:
    excel_sheet.iloc[rowname, 3] = '{0}, {1}'.format(power[0], power[1])
else:
    excel_sheet.iloc[rowname, 3] = power
excel_sheet.iloc[rowname, 4] = tau_min
excel_sheet.iloc[rowname, 5] = tau_max
if type(shell_thickness)==tuple:
    excel_sheet.iloc[rowname, 6] = '{0}, {1}'.format(shell_thickness[0], shell_thickness[1])
else:
    excel_sheet.iloc[rowname, 6] = shell_thickness

components= [comp1, comp2, comp3, comp4, comp5, comp6, comp7]
dust_recipe= {components[i]: abundance[i] for i in range(len(abundance)) if abundance[i]>0.0}

silicates = {'sil ow':sil_ow, 'sil oc':sil_oc, 'sil dl': sil_dl,
             'grf dl': grf_dl, 'amc hn': amc_hn, 'sic pg': sic_pg}

silicate_recipe = {i: silicates[i] for i in silicates if silicates[i]>0.0}

dust_total = {**silicate_recipe, **dust_recipe}
li = list(dust_total.items())
nms = [x[0] for x in li]
abn = [x[1] for x in li]
abnd = [round(g/sum(abn),3) for g in abn]

j = 7

for i in range(len(dust_total) ):
    excel_sheet.iloc[rowname, j] = nms[i]
    excel_sheet.iloc[rowname, j+1] = abnd[i]
    j += 2

excel_sheet.to_excel('Model Parameters.xlsx', index = False, sheet_name='Sheet3')



















































print('hello world')
