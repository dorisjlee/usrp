import shutil
import numpy as np
import h5py
import glob
N_PHI = 256
NUM_MESHBLOCK=264
def compute_fastest_growing_mode(tstep_filename,debug=False):
    if debug: print "Working on {}".format(tstep_filename)
    hdf5 = h5py.File(tstep_filename, 'r+')
    D = np.zeros(N_PHI)
    for M in np.arange(NUM_MESHBLOCK):
        meshblock =hdf5["MeshBlock{}".format(M)]
        density= meshblock["rho"].value#density for each meshblock
        N_r = density.shape[2]
        N_theta=density.shape[1]
        N_phi= density.shape[0]
        x1f = meshblock["x1f"].value
        x2f = meshblock["x2f"].value
        level = meshblock.attrs["Level"][0]
        if level ==0:
            for n in np.arange(N_phi):
                loc = meshblock.attrs["LogicalLocation"]
                phi_i =(n+32*loc[2])*2
                #Compute D_i in each meshblock
                m=0#mass summed over r,theta per density slice
                for ti in np.arange(N_theta):
                    thetai = (x2f[ti]+x2f[ti+1])*0.5
                    dthetai = abs(x2f[ti]-x2f[ti+1])
                    for rii in np.arange(N_r):
                        #for each cell compute the coordinate to calculate the mass in each cell
                        ri = (x1f[rii]+x1f[rii+1])*0.5
                        dri = abs(x1f[rii]-x1f[rii+1])
                        m=m+density[n,ti,rii]*(ri*ri)*np.sin(thetai)*dri*dthetai
                D[phi_i]=D[phi_i]+m*0.5
                D[phi_i+1]=D[phi_i+1]+m*0.5
        else: #level=1
            for n in np.arange(N_phi):
                loc = meshblock.attrs["LogicalLocation"]
                phi_i =n+32*loc[2]
                #Compute D_i in each meshblock
                m=0#mass summed over r,theta per density slice
                for ti in np.arange(N_theta):
                    thetai = (x2f[ti]+x2f[ti+1])*0.5
                    dthetai = abs(x2f[ti]-x2f[ti+1])
                    for rii in np.arange(N_r):
                        #for each cell compute the coordinate to calculate the mass in each cell
                        ri = (x1f[rii]+x1f[rii+1])*0.5
                        dri = abs(x1f[rii]-x1f[rii+1])
                        m=m+density[n,ti,rii]*(ri*ri)*np.sin(thetai)*dri*dthetai
                D[phi_i]=D[phi_i]+m
    if debug: print "D:",D
    return D
for f in np.sort(glob.glob("mri_hires.out1.*.athdf")):
    file = open('amp.txt', 'a')
    density =compute_fastest_growing_mode(f,True)
    dens = np.fft.rfft(density)
    amp = np.real(np.sqrt(dens*np.conj(dens)))
    file.write("{}         {}           {}          {}            {}		{}             {}           {}         {}           {}          \n".format(amp[0],amp[1],amp[2],amp[3],amp[4],amp[5],amp[6],amp[7],amp[8],amp[9]))
    shutil.move(f, "done/")
    file.close()
