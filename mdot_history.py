import numpy as np
import h5py
import glob
NUM_MESHBLOCK=96
def compute_mass_accretion(tstep_filename):
    debug =True
    if debug: print "Working on {}".format(tstep_filename) 
    hdf5 = h5py.File(tstep_filename, 'r+')
    D = np.zeros(128)
    for M in np.arange(NUM_MESHBLOCK):
    #     print "Working on MeshBlock {}".format(M)
        meshblock =hdf5["MeshBlock{}".format(M)]
        density= meshblock["rho"].value#density for each meshblock
        N_r = density.shape[2]
        N_theta=density.shape[1]
        N_phi= density.shape[0]
        x1f = meshblock["x1f"].value
        x2f = meshblock["x2f"].value
        x3f = meshblock["x3f"].value
        vr = meshblock["vel1"].value
        loc = meshblock.attrs["LogicalLocation"]
    #     for rii in np.arange(N_r):
    #         phi_i =n+32*loc[2]
            #Compute D_i in each meshblock
    #         if rii = 0
        if (loc[0]==0):
            rii=0 # we only want the mass flux onto the inner most cell
            val=0 #val=integrand to sum over rho*vr*dS
            ri = (x1f[rii]+x1f[rii+1])*0.5
            dri = abs(x1f[rii]-x1f[rii+1])
            for ti in np.arange(N_theta):
                tii =ti+32*loc[2]
                thetai = (x2f[ti]+x2f[ti+1])*0.5
                dthetai = abs(x2f[ti]-x2f[ti+1])
                for pii in np.arange(N_phi):
                    dphii=abs(x3f[ti]-x3f[ti+1])
                    dS= ri*np.sin(thetai)*dri*dphii
                    val = val +density[pii,ti,rii]*vr[pii,ti,rii]*dS
                D[tii]=D[tii]+val
    return sum(D)
for f in  np.sort(glob.glob("mri_hires_smaller.out1.*.ath")):
    file = open('mdot.txt', 'a')
    mdot = compute_mass_accretion(f)
    print f , mdot
    file.write("{}\n".format(mdot))
    file.close()
