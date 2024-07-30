# %% [markdown]
# # Calculating the absorption and scattering cross section using a probability distribution

# Max TO-DO
# Integrate with new Sean code DONE
# Individual radius sliders for each dust
# Allow selecting less than max dust DONE
# Checkboxes for each dist DONE
# More saving options
# Fix CDE-tCDE
# Save logs to root folder
# %%
#import time
# AA = time.time()
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spit
import os
import tkinter.filedialog
import tkinter as tk
from tkinter import ttk
from ttkthemes import ThemedTk
import configparser

print('Imports done')


# %% [markdown]
# ### This defines the probability distribution


config = configparser.ConfigParser()

def save_config():
    config['GUI SETTINGS'] = {
        'Dust Directory': dust_dir.get(),
    }
    with open('gui_config.ini', 'w') as configfile:
        config.write(configfile)

def load_config():
    if os.path.exists('gui_config.ini'):
        config.read('gui_config.ini')
        dust_dir.set(config['GUI SETTINGS'].get('Dust Directory', ''))

def save_and_exit():
    save_config()
    root.destroy()

def select_directory():
    dir_name = tkinter.filedialog.askdirectory()
    dust_dir.set(dir_name)
    print(dust_dir.get())


def select_dust(index):
    file_path = tkinter.filedialog.askopenfilename(initialdir=dust_dir.get(), filetypes=[("nk files", "*.nk")])
    file_name = os.path.basename(file_path)
    dusts[index].set(file_name)
    print(dusts[index].get())


def add_dust():
    global dust_count
    dust_count += 1
    dust_var = tk.StringVar()
    weight_var = tk.DoubleVar()
    dist_var = tk.StringVar(value="spheres")
    minr_var = tk.DoubleVar(value=0.005)
    maxr_var = tk.DoubleVar(value=0.25)
    q_var = tk.DoubleVar(value=3.5)
    dusts.append(dust_var)
    weights.append(weight_var)
    distributions.append(dist_var)
    minrs.append(minr_var)
    maxrs.append(maxr_var)
    qs.append(q_var)
    row = dust_count + 1

    ttk.Label(root, text=f"Dust: {dust_count}").grid(row=row, column=0, padx=10, pady=5)
    ttk.Entry(root, textvariable=dust_var, state=tk.DISABLED).grid(row=row, column=1, padx=10, pady=5)
    ttk.Button(root, text="Select", command=lambda index=dust_count-1: select_dust(index)).grid(row=row, column=2, padx=10, pady=5)
    ttk.Entry(root, textvariable=weight_var, width=5).grid(row=row, column=3, padx=10, pady=5)
    ttk.OptionMenu(root, dist_var, "spheres", "spheres", "CDE", "CDE2", "ERCDE", "tCDE").grid(row=row, column=4, padx=10, pady=5)
    ttk.Entry(root, textvariable=minr_var, width=7).grid(row=row, column=5, padx=10, pady=5)
    ttk.Entry(root, textvariable=maxr_var, width=7).grid(row=row, column=6, padx=10, pady=5)
    ttk.Entry(root, textvariable=q_var, width=5).grid(row=row, column=7, padx=10, pady=5)

    add_button.grid(row=dust_count + 2, column=0, padx=10, pady=10)
    remove_button.grid(row=dust_count + 2, column=1, padx=10, pady=10)

def remove_dust():
    global dust_count
    if dust_count > 1:

        for widget in root.grid_slaves(row=dust_count + 1):
            widget.grid_forget()

        dusts.pop()
        weights.pop()
        distributions.pop()
        minrs.pop()
        maxrs.pop()
        qs.pop()
        dust_count -= 1

        add_button.grid(row=dust_count + 2, column=0, padx=10, pady=10)
        remove_button.grid(row=dust_count + 2, column=1, padx=10, pady=10)
    
# %%
def run_program():
    print('Program started')
    def probability(dis_name, l1, l2, lmin=0.05, m1=0, m2=0, d=0):
        '''
        This is the probability distribution as a function of L1 and L2, the 
        geometric parameters. This parameter gets inserted into the integral that 
        calculates the average polarizability per unit volume

        Parameters
        ----------
        dis_name : String
            This specifies the distribution we will be using. 
            'CDE' = Continuous Distribution of Ellipsoids
            'ERCDE' = Externally Restricted CDE, returns CDE if lmin=0
            'tCDE' = truncated CDE, REQUIRES MORE WORK
        l1 : Float
            Largest Geometric Constant, lmin<l1<1.0
        l2 : Float
            Second Largest Geometric Constant, lmin<l2<=l1
        lmin : Float, optional
            Minimum allowed geometric constant.  The default is 0.

        Returns
        -------
        Float
            Function dependent on l1 and l2

        '''
        l3 = 1 - l1 - l2
        if dis_name == 'CDE':
            return 2
        elif dis_name == 'CDE2':
            return 120 * l1 * l2 * l3
        elif dis_name == 'ERCDE':
            return 2/((1 - (3*lmin))**2)
        elif dis_name == 'tCDE':
            return 1/((1-d-m2)*(1-m1-m2-d) - 0.5*((1-d-m2)**2) - m1**2)
        else:
            return True
            
        
    print('probability done')

    # %% [markdown]
    # ### I think it would be easier to make the tCDE its own function

    # %%
    def prob_tCDE(m1,m2,d):
        return 1/((1-d-m2)*(1-m1-m2-d) - 0.5*((1-d-m2)**2) - m1**2)

    # %% [markdown]
    # ### Now we define our volume function. v_avg is used to calculate sigma below

    # %%
    def volume_integrand_mrn(r, q):
        v = r**(-q)
        return v

    # UNITS ARE IN MICRONS
   # rmin = 0.005
   # rmax = 0.25
   # q = 3.5





    # %% [markdown]
    # ### creates function that calculates Sigma, defined in eq 13 in Min et al 2003

    # %%
    def sigma(m, lamda, v):
        sig = []
        for i in range(len(lamda)):
            k = (2.0 * np.pi)/lamda[i]
            term1 = (6.0*np.pi) / (v * (k**3))
            term2 = np.imag((m[i]**2))
            term3 = 1.0 / abs(m[i]**2 - 1)**2
            sig.append(term1 * term2 * term3)
        return sig
    print(minrs, maxrs, qs)

    # %% [markdown]
    # ### creates our bounds for our geometric factors. The bounds are the sides of a triangle in (l1, l2) space

    # %%
    def bounds_l1():
        return [0,1]

    def bounds_l2(l1):
        return [0,1-l1]
    print('bounds done')


    # %% [markdown]
    # ### This is where we calculate the absorption cross-section (Cabs). It creates an empty list, then calculates Cabs for a given distribution at each wavelength as described in Min 03, eqn 15. It then uses this to find the shape averaged mass absorption coefficient for particles of a given volume, as described in Min 03, eqn 39

    # %%
    def cabs(m, dis_name, bounds_l2, bounds_l1):
        cabs = []
        if dis_name=='spheres':
            for j in range(len(m)):
                cabs.append(np.imag(3*(m[j]**2 - 1)/(m[j]**2 + 2)))
        else:
            for j in range(len(m)):
                def f(l1, l2, n=m[j], dis_name='CDE'):
                    b = 1/(n**2 - 1)
                    term1 = 1/3 * 1/(b + l1)
                    term2 = 1/3 * 1/(b + l2)
                    term3 = 1/3 * 1/(b + 1 - l1 - l2)
                # r = np.real((term1 + term2 + term3)*probability(dis_name, l1, l2))
                    j = np.imag((term1 + term2 + term3)*probability(dis_name, l1, l2))
                    return j
                # return np.real((term1 + term2 + term3)*probability(dis_name, l1, l2)) + np.imag((term1 + term2 + term3)*probability(dis_name, l1, l2))
                cabs.append(spit.nquad(f, [bounds_l2, bounds_l1])[0])
        return cabs

    # j = cabs(m, 'spheres', bounds_l2, bounds_l1)
    # # dust = 'grph1-dl.nk'                  #DUST NAME HERE #grf
  #  rho = 3.33 #grams cm**-3            #density
    # # pathy = os.path.join(nk_path, dust) #pipeline is open
    # # wavelen, n_dust, k_dust = np.loadtxt(pathy, skiprows=7, unpack=True)
    # #                                     #lamda, n, and k values are extracted
    # # m = np.array([complex(n_dust[i], k_dust[i]) for i in range(len(wavelen))])
    # k = cabs(m, 'CDE', bounds_l2, bounds_l1)

    # print(j)
    # # print('')
    # print(k)

    # %% [markdown]
    # 
    # 
    # 
    # 
    # 
    # 
    # # Specify which dust we are using

    # %%
   # dustlist = [(dust1.get(), dist_a.get()), (dust2.get(), dist_b.get()), (dust3.get(), dist_c.get())]

   # weightlist = [wt_a.get(), wt_b.get(), wt_c.get()]
    dustlist = [(dusts[j].get(), distributions[j].get()) for j in range(dust_count)]
    print(dustlist)
    namelist = [dustlist[j][0][:-3]+dustlist[j][1]+'.dat' for j in range(len(dustlist))]
    weightlist = [weights[j].get() for j in range(dust_count)]
    rminlist = [minrs[j].get() for j in range(dust_count)]
    rmaxlist = [maxrs[j].get() for j in range(dust_count)]
    qlist = [qs[j].get() for j in range(dust_count)]
    #dust_dir = ['/home/physics/Research/DUSTY/DUSTY/Lib_nk/', 
    #            "C:/Users/Max/Documents/DUSTY/dusty-master/data/Lib_nk/"]
    # this is the possible locations of where dust can be
    print(rminlist,rmaxlist,qlist)

    #nk_path = dust_dir[1]               #where the dust is 
    nk_path = dust_dir.get()      
    # lam_max, cabs_max, csca_max = np.loadtxt(max(namelist, key=os.path.getsize), unpack=True)

    # output = np.zeros((len(dustlist), len(lam_max), 3))


    for j in range(len(dustlist)):
        pathy = os.path.join(nk_path, dustlist[j][0]) #pipeline is open
        r_integral = spit.quad(volume_integrand_mrn, rminlist[j], rmaxlist[j], args=qlist[j])
        r_average = ((1/(rmaxlist[j] - rminlist[j])) * r_integral[0])**(1/-qlist[j])
        v_avg = (4./3.) * np.pi * r_average**3
        wavelen, n_dust, k_dust = np.loadtxt(pathy, skiprows=8, unpack=True)
        m = np.array([complex(n_dust[i], k_dust[i]) for i in range(len(wavelen))])
        cab = cabs(m, dustlist[j][1], bounds_l2, bounds_l1)
        Cabs_array = np.array((cab))
        Cabs_array *= (2 * np.pi / (wavelen)) * v_avg
        sig = np.array((sigma(m, wavelen, v_avg)))
        Csca_array = Cabs_array/sig
        output = np.transpose((wavelen, Cabs_array, Csca_array))
        f = open(dustlist[j][0][:-3]+dustlist[j][1]+'.dat', 'w')
        for i in range(len(output)):
            f.write(f"{output[i,0]} \t {output[i,1]} \t {output[i,2]}\n")
        f.close()
        print(pathy)
        
        


    lam_final = np.geomspace(0.2, 500, num=500)
    total_array = np.ndarray((len(dustlist), len(lam_final), 3))
    total_array[:,:,0] = lam_final

    for k in range(len(namelist)):
        lam, cabs_tot, csca_tot = np.loadtxt(namelist[k], unpack=True)
        total_array[k,:,1] = np.interp(lam_final, lam, cabs_tot)
        total_array[k,:,2] = np.interp(lam_final, lam, csca_tot)


        
    avg_array = np.ndarray((len(lam_final),3))
    avg_array[:,0] = lam_final
    for j in range(len(lam_final)):
        avg_array[j,1] = np.average(total_array[:,j,1], weights=weightlist)
        avg_array[j,2] = np.average(total_array[:,j,2], weights=weightlist)

        
        
    titlestring=''
    for g in range(len(namelist)):
        titlestring += namelist[g][:3] + str(weightlist[g]).replace('.','')
        
    f = open(titlestring+'.dat','w')
    for i in range(len(lam_final)):
        f.write(f"{avg_array[i,0]} \t {avg_array[i,1]} \t {avg_array[i,2]}\n")
    f.close() 
        
        
    print(titlestring)


    # %%


    # %%
    # lam_new, cab_new, csa_new = np.loadtxt('ism.dat', unpack=True, skiprows=3)
    lam_old, cab_old, csa_old = np.loadtxt(titlestring+'.dat', unpack=True)

    fig, ax = plt.subplots()
    y_axis = yaxis.get()
    if y_axis == "C_abs":
        y_axis = cab_old
        ax.set_ylabel(r'$C_{\mathrm{abs}}$ (cm$^2$/g)', fontsize=14)
        title = 'x=lambda y=c_abs'
    else:
        y_axis = csa_old
        ax.set_ylabel(r'$C_{\mathrm{sca}}$ (cm$^2$/g)', fontsize=14)
        title = 'x=lambda y=c_sca'

    x_min = min(lam_old)
    x_max = max(lam_old)
    y_min = min(y_axis)
    y_max = max(y_axis)
    #fig, ax = plt.subplots()
    #ax.set(xscale='linear', yscale='log')

    #title = 'Lorem ipsum dolor sit amet. '
    #ax.set(xscale='log', yscale='log', xlim=(2,500))
    #ax.set_title(title, fontsize=16)
    #ax.set_xlabel(r'$\lambda (\mu m)$', fontsize=14)
    #ax.set_ylabel(r'$<C_{abs}>$ cm$^{2}$ g$^{-1}$', fontsize=14)
    #ax.plot(wavelen, avg_array, label='please')
    #ax.legend()

    title = 'Lorem ipsum dolor sit amet. '
    ax.set(xscale='log', yscale='log', xlim=(2,500))
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_title(title, fontsize=16)
    ax.set_xlabel(r'$\lambda (\mu m)$', fontsize=14)

    #ax.plot(lam_new, csa_new, label='new')
    ax.plot(lam_old, y_axis, label='old')
    # ax.plot(wavelen, np.array((cabs_ercde)), label='ERCDE')
    ax.legend()
    plt.tight_layout()

    plt.show()
    print("DONE!")

# %%
#tk.Scale(root, from_=0, to=100, orient=tk.HORIZONTAL, variable=wt_b).grid(row=5, column=1)
root = ThemedTk(theme="ubuntu")
root.title("Dust Properties")

dust_dir = tk.StringVar()
dusts = [tk.StringVar()]
weights = [tk.DoubleVar(value=0.0)]
minrs = [tk.DoubleVar(value=0.005)]
maxrs = [tk.DoubleVar(value=0.25)]
qs = [tk.DoubleVar(value=3.5)]
distributions = [tk.StringVar(value="spheres")]
yaxis = tk.StringVar(value="C_abs")
dust_count = len(dusts)

ttk.Label(root, text="Dust Directory:").grid(row=0, column=0, padx=10, pady=5)
ttk.Entry(root, textvariable=dust_dir, state=tk.DISABLED, width=60).grid(row=0, column=1, columnspan=3, padx=10, pady=5)
ttk.Button(root, text="Select", command=select_directory).grid(row=0, column=4, padx=10, pady=5)

ttk.Label(root, text="File").grid(row=1, column=1, padx=10, pady=5)
ttk.Label(root, text="Weight").grid(row=1, column=3, padx=10, pady=5)
ttk.Label(root, text="Distribution").grid(row=1, column=4, padx=10, pady=5)
ttk.Label(root, text="Radius Range (µm)").grid(row=1, column=5, columnspan=2, padx=0, pady=5)
ttk.Label(root, text="Power Law\n Exponent").grid(row=1, column=7, padx=10, pady=5)

ttk.Label(root, text="Dust 1:").grid(row=2, column=0, padx=10, pady=5)
ttk.Entry(root, textvariable=dusts[0], state=tk.DISABLED).grid(row=2, column=1, padx=10, pady=5)
ttk.Button(root, text="Select", command=lambda index=0: select_dust(index)).grid(row=2, column=2, padx=10, pady=5)
ttk.Entry(root, textvariable=weights[0], width=5).grid(row=2, column=3, padx=10, pady=5)
ttk.OptionMenu(root, distributions[0], "spheres", "spheres", "CDE", "CDE2", "ERCDE", "tCDE").grid(row=2, column=4, padx=10, pady=5)
ttk.Entry(root, textvariable=minrs[0], width=7).grid(row=2, column=5, padx=10, pady=5)
ttk.Entry(root, textvariable=maxrs[0], width=7).grid(row=2, column=6, padx=10, pady=5)
ttk.Entry(root, textvariable=qs[0], width=5).grid(row=2, column=7, padx=10, pady=5)

add_button = ttk.Button(root, text="Add Dust", command=add_dust).grid(row=98, column=0, padx=10, pady=10)
remove_button = ttk.Button(root, text="Remove Dust", command=remove_dust).grid(row=98, column=1, padx=10, pady=10)

ttk.Label(root, text="Y-Axis:").grid(row=99, column=0, padx=10, pady=5)
ttk.OptionMenu(root, yaxis, "C_abs", "C_abs", "C_sca").grid(row=99, column=1, padx=10, pady=5)

run_button = ttk.Button(root, text="Run Program", command=run_program, style='Accent.TButton', state=tk.ACTIVE)
run_button.grid(row=100, column=100, padx=10, pady=10)

load_config()
root.protocol("WM_DELETE_WINDOW", save_and_exit)

root.mainloop()

 
# %%
