# %% [markdown]
# # Calculating the absorption and scattering cross section using a probability distribution

# %%
import time
# AA = time.time()
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spit
import os
import tkinter.filedialog
import tkinter as tk
print('hello world')


# %% [markdown]
# ### This defines the probability distribution
# Function to select directory
def select_directory():
    dir_name = tkinter.filedialog.askdirectory()
    dust_dir.set(dir_name)
    print(dust_dir.get())

# Functions to select dust files
def select_dust1():
    file_path = tkinter.filedialog.askopenfilename(initialdir=dust_dir.get(), filetypes=[("nk files", "*.nk")])
    file_name = os.path.basename(file_path)
    dust1.set(file_name)
    print(dust1.get())

def select_dust2():
    file_path = tkinter.filedialog.askopenfilename(initialdir=dust_dir.get(), filetypes=[("nk files", "*.nk")])
    file_name = os.path.basename(file_path)
    dust2.set(file_name)
    print(dust2.get())

def select_dust3():
    file_path = tkinter.filedialog.askopenfilename(initialdir=dust_dir.get(), filetypes=[("nk files", "*.nk")])
    file_name = os.path.basename(file_path)
    dust3.set(file_name)
    print(dust3.get())
# %%
def run_program():
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
            
        
    print('hello world')

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
    rmin = 0.005
    rmax = 0.25
    q = 3.5

    r_integral = spit.quad(volume_integrand_mrn, rmin, rmax, args=q)
    r_average = ((1/(rmax - rmin)) * r_integral[0])**(1/-q)
    v_avg = (4./3.) * np.pi * r_average**3



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
    print('hello world')

    # %% [markdown]
    # ### creates our bounds for our geometric factors. The bounds are the sides of a triangle in (l1, l2) space

    # %%
    def bounds_l1():
        return [0,1]

    def bounds_l2(l1):
        return [0,1-l1]
    print('hello world')


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
    rho = 3.33 #grams cm**-3            #density
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
    dustlist = [(dust1.get(), 'spheres'), (dust2.get(), 'spheres'), (dust3.get(), 'spheres')]
    namelist = [dustlist[j][0][:-3]+dustlist[j][1]+'.dat' for j in range(len(dustlist))]
    weightlist = [wt_a.get(), wt_b.get(), wt_c.get()]

    #dust_dir = ['/home/physics/Research/DUSTY/DUSTY/Lib_nk/', 
    #            "C:/Users/Max/Documents/DUSTY/dusty-master/data/Lib_nk/"]
    # this is the possible locations of where dust can be


    #nk_path = dust_dir[1]               #where the dust is 
    nk_path = dust_dir.get()      
    # lam_max, cabs_max, csca_max = np.loadtxt(max(namelist, key=os.path.getsize), unpack=True)

    # output = np.zeros((len(dustlist), len(lam_max), 3))


    for j in range(len(dustlist)):
        pathy = os.path.join(nk_path, dustlist[j][0]) #pipeline is open
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

        
        


    lam_final = np.geomspace(0.2, 500, num=500)
    total_array = np.ndarray((3,len(lam_final),len(dustlist)))
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

    x_min = min(lam_old)
    x_max = max(lam_old)
    y_min = min(csa_old)
    y_max = max(csa_old)

    #fig, ax = plt.subplots()
    #ax.set(xscale='linear', yscale='log')

    #title = 'Lorem ipsum dolor sit amet. '
    #ax.set(xscale='log', yscale='log', xlim=(2,500))
    #ax.set_title(title, fontsize=16)
    #ax.set_xlabel(r'$\lambda (\mu m)$', fontsize=14)
    #ax.set_ylabel(r'$<C_{abs}>$ cm$^{2}$ g$^{-1}$', fontsize=14)
    #ax.plot(wavelen, avg_array, label='please')
    #ax.legend()
    fig, ax = plt.subplots()
    title = 'Lorem ipsum dolor sit amet. '
    ax.set(xscale='log', yscale='log', xlim=(2,500))
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_title(title, fontsize=16)
    ax.set_xlabel(r'$\lambda (\mu m)$', fontsize=14)
    ax.set_ylabel(r'$<C_{abs}>$ cm$^{2}$ g$^{-1}$', fontsize=14)
    #ax.plot(lam_new, csa_new, label='new')
    ax.plot(lam_old, csa_old, label='old')
    # ax.plot(wavelen, np.array((cabs_ercde)), label='ERCDE')
    ax.legend()
    plt.tight_layout()

    plt.show()



# %%
root = tk.Tk()
root.title("GUI still wip")

dust_dir = tk.StringVar()
dust1 = tk.StringVar()
dust2 = tk.StringVar()
dust3 = tk.StringVar()

tk.Label(root, text="Dust File Directory").grid(row=0, column=0)
tk.Entry(root, textvariable=dust_dir).grid(row=0, column=1)
tk.Button(root, text="Select Directory", command=select_directory).grid(row=0, column=2)

tk.Label(root, text="Dust File 1").grid(row=1, column=0)
tk.Entry(root, textvariable=dust1).grid(row=1, column=1)
tk.Button(root, text="Select Dust 1", command=select_dust1).grid(row=1, column=2)

tk.Label(root, text="Dust File 2").grid(row=2, column=0)
tk.Entry(root, textvariable=dust2).grid(row=2, column=1)
tk.Button(root, text="Select Dust 2", command=select_dust2).grid(row=2, column=2)

tk.Label(root, text="Dust File 3").grid(row=3, column=0)
tk.Entry(root, textvariable=dust3).grid(row=3, column=1)
tk.Button(root, text="Select Dust 3", command=select_dust3).grid(row=3, column=2)

tk.Label(root, text="Set Weights for Each Dust File").grid(row=4, column=0)

wt_a = tk.DoubleVar()
wt_b = tk.DoubleVar()
wt_c = tk.DoubleVar()

tk.Scale(root, from_=0, to=100, orient=tk.HORIZONTAL, variable=wt_a).grid(row=5, column=0)
tk.Scale(root, from_=0, to=100, orient=tk.HORIZONTAL, variable=wt_b).grid(row=5, column=1)
tk.Scale(root, from_=0, to=100, orient=tk.HORIZONTAL, variable=wt_c).grid(row=5, column=2)

tk.Button(root, text="Run Program", command=run_program).grid(row=6, column=1)

root.mainloop()


