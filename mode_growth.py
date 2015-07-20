import numpy as np
import h5py
import glob
def compute_fastest_growing_mode(tstep_filename,debug=False):
    if debug: print "Working on {}".format(tstep_filename) 
    hdf5 = h5py.File(tstep_filename, 'r+')
    NUM_MESHBLOCK=24
    D = np.zeros(128)
    for M in np.arange(NUM_MESHBLOCK):
        meshblock =hdf5["MeshBlock{}".format(M)]
        density= meshblock["rho"].value#density for each meshblock
        N_r = density.shape[0]
        N_theta=density.shape[1]
        N_phi= density.shape[2]
        for n in np.arange(N_phi):
            loc = meshblock.attrs["LogicalLocation"]
            phi_i =n+32*loc[2]
            #Compute D_i in each meshblock
            x1f = meshblock["x1f"].value
            x2f = meshblock["x2f"].value
            m=[]#mass summed over r,theta per density slice
            for ri in np.arange(N_r):
                for ti in np.arange(N_theta):
                    ri = (x1f[ri]+x1f[ri+1])/2.
                    dri = abs(x1f[ri]-x1f[ri+1])
                    thetai = (x2f[ti]+x2f[ti+1])/2.
                    dthetai = abs(x2f[ti]-x2f[ti+1])
                    m.append(density[ri,ti,n]*(ri**2)*np.sin(thetai)*dri*dthetai)
            D[phi_i]=D[phi_i]+sum(m)
    #if debug: print "D:",D
    amp = float(np.real(max(np.fft.fft(D))))
    if debug: print "amp:", amp
    return amp
mode_amp = []
for f in np.sort(glob.glob("log.out1.000*.ath")):
    file = open('amp.txt', 'a')
    val =  compute_fastest_growing_mode(f,True)
    mode_amp.append(val)
    file.write("{}\n".format(val))
    file.close()
print mode_amp
#plt.plot(mode_amp)
#plt.ylabel("Mode Amplitude",fontsize=13)
#plt.xlabel("Time",fontsize=13)
#plt.show()
