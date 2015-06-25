from ath_read import *
from ath_plots import *
from ath_styles import *
from ath_units import *
import matplotlib.pyplot as plt
import matplotlib as mpl

# Data parameters
#datapath = '/Users/askinner/Work/codes/Hyperion/bin/radplanesrc_kt12-2D_test'
datapath = '/Users/askinner/Downloads'
savepath = datapath
#basename = 'RadPlaneSrc_kt12'
basename = 'RadPlaneSrc'
stepmin = 50
stepmax = 50
skip = 1

# Specify which plots to make
plot_lineout_all       = 1
#----------------------------
plot_lineout_rho       = 1
plot_lineout_f         = 0
#============================
plot_pcolor_all        = 1
#----------------------------
plot_pcolor_rho        = 0
plot_pcolor_Trad       = 0
#----------------------------

# Set plotting/saving parameters
saveflag = 0
savetype = 'png'
dpi = 400
cmap = 'afmhot'
mpl.rcdefaults()
if (saveflag):
    set_style()

# Set code units    
h_star   = 7.874994091733578125e+13
rho_star = 5.938676856256753506e-14
c_star   = 5.396264243999999599e+04
T_star   = 8.222995156475339229e+01
units = Units_LDVK(h_star,rho_star,c_star,T_star)
## Set physical units
#units = Units_LMT()

# Create figure windows
plot_lineout_list = np.array([plot_lineout_rho,
                              plot_lineout_f])
plot_pcolor_list = np.array([plot_pcolor_rho,
                             plot_pcolor_Trad])
plt.close('all')
#plt.ion()
nfig = 0
for i in range(plot_lineout_list.size):
    if plot_lineout_list[i]==1 or plot_lineout_all==1:
        plt.figure(nfig)
        nfig += 1
for i in range(plot_pcolor_list.size):
    if plot_pcolor_list[i]==1 or plot_pcolor_all==1:
        plt.figure(nfig)
        nfig += 1
print 'Creating {0:d} figure windows...'.format(nfig)

for step in range(stepmin,stepmax+1,skip):
    print 'step',step
    ifig = 0
    fignames = []
    
    data = AthenaData(datapath,basename,step,units)
    time = float(data.time)

    # Lineout plots
    # -------------------------------------------------------------------------
    if plot_lineout_rho or plot_lineout_all:
        field = 'rho'
        fignames.append('lineout.{0:s}.{1:04d}'.format(field,step))
        plt.figure(ifig)
        ifig += 1
        figh = lineout(data,field,logy=False,ycut=2.0,zcut=0.0)
        plt.xlabel(r'$z$')
        plt.ylabel(r'$\rho$')
        plt.title('Step {0:d}, Time = {1:0.3f}'.format(step,time))

    if plot_lineout_f or plot_lineout_all:
        field = 'f'
        fignames.append('lineout.{0:s}.{1:04d}'.format(field,step))
        plt.figure(ifig)
        ifig += 1
        figh = lineout(data,field,logy=False,ycut=2.0,zcut=0.0)
        plt.xlabel(r'$z$')
        plt.ylabel(r'$f$')
        plt.title('Step {0:d}, Time = {1:0.3f}'.format(step,time))

    # Pseudocolor plots
    # -------------------------------------------------------------------------
    if plot_pcolor_rho or plot_pcolor_all:
        field = 'rho'
        fignames.append('pcolor.{0:s}.{1:04d}'.format(field,step))
        plt.figure(ifig)
        plt.clf()
        ifig += 1
        figh = pcolor(data,field,cmap,kcut=0,transpose=True,logz=True,clim=(1.0e-3,3.0))
        plt.xlabel(r'$x$')
        plt.ylabel(r'$z$')
        cb = plt.colorbar()
        cb.set_label(r'$\rho/\rho_*$')
        plt.title('Step {0:d}, Time = {1:0.3f}'.format(step,time))

    if plot_pcolor_Trad or plot_pcolor_all:
        field = 'Trad'
        fignames.append('pcolor.{0:s}.{1:04d}'.format(field,step))
        plt.figure(ifig)
        plt.clf()
        ifig += 1
        figh = pcolor(data,field,cmap,kcut=0,transpose=True)
        plt.xlabel(r'$x$')
        plt.ylabel(r'$z$')
        cb = plt.colorbar()
        cb.set_label(r'$T_\mathrm{rad}/T_*$')
        plt.title('Step {0:d}, Time = {1:0.3f}'.format(step,time))

    if saveflag:
        for i in range(nfig):
            filename = '{0:s}/{1:s}.{2:s}.{3:s}'.format(savepath,basename,fignames[i],savetype)
            if   savetype=='png':
                plt.savefig(filename,dpi=dpi,bbox_inches='tight')
            elif savetype=='eps':
                plt.savefig(filename)
            print 'Figure {0:s} saved!'.format(filename)
    else:
        print 'done'
        plt.show()