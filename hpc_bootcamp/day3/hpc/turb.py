from ath_read import *
from ath_plots import *
from ath_styles import *
import matplotlib.pyplot as plt
import matplotlib as mpl

datapath = '..'
savepath = datapath
basename = 'Turb_joined'
stepmin = 50
stepmax = 50
skip = 1

saveflag = 0
savetype = 'png'
dpi = 400
cmap = 'afmhot'

# Set colormap
nfig = 2
for i in range(nfig):
    plt.figure(i)
cmap  = mpl.cm.get_cmap('jet')
cnorm = mpl.colors.Normalize(vmin=stepmin,vmax=stepmax)
cm    = mpl.cm.ScalarMappable(norm=cnorm,cmap=cmap)

for step in range(stepmin,stepmax+1,skip):
    ifig = 0
    print 'step',step
    
    data = AthenaData(datapath,basename,step,units)

    plt.figure(ifig)
    ifig += 1
    pcolor(data,field='density',cmap=cmap,logv=True,vmin=1.0e-6,vmax=1.3,zcut=0.0)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    if saveflag:
        figname = '{0:s}.density.{1:04d}.{2:s}'.format(basename,step,savetype)
        if   savetype=='png':
            plt.savefig(savepath+'/'+figname,dpi=dpi,bbox_inches='tight')
        elif savetype=='eps':
            plt.savefig(savepath+'/'+figname)
        print 'Figure {0:s} saved!'.format(figname)
        plt.clf()
    else:
        print 'Close figure to continue...'
        plt.show()

#    plt.figure(ifig)
#    ifig += 1
#    line = lineout(data,field='density',logy=True,ycut=0.0,zcut=0.0)
#    plt.setp(line,color=cm.to_rgba(step),antialiased=True)
#    Tgas = data.get_field('Tgas')
#    rho = data.get_field('rho')
#    print 'Tgas[0]=',Tgas[0][0]
#    print 'rho [0]=',rho[0][0]
#    z,Tgas = data.get_lineout_xyz('Tgas',xcut=5.0,zcut=0.0)
#    z,rho  = data.get_lineout_xyz('rho' ,xcut=5.0,zcut=0.0)
#    print 'Tgas[-1]=',Tgas[0]
#    print 'rho [-1]=',rho[0]
##    print Tgas
##    print rho
#    plt.xlim(0,40)
#    plt.ylim(1.0e-8,2.0)
#    plt.xlabel(r'$z$')
#    plt.ylabel(r'$\rho$')
#
#    plt.figure(ifig)
#    ifig += 1
#    line = lineout(data,field='pressure',logy=True,ycut=0.0,zcut=0.0)
#    plt.setp(line,color=cm.to_rgba(step),antialiased=True)
#    plt.xlim(0,40)
#    plt.ylim(1.0e-8,2.0)
#    plt.xlabel(r'$z$')
#    plt.ylabel(r'$P$')
#
#    plt.figure(ifig)
#    ifig += 1
#    line = lineout(data,field='temperature',logy=True,ycut=0.0,zcut=0.0)
#    plt.setp(line,color=cm.to_rgba(step),antialiased=True)
##    plt.xlim(0,40)
##    plt.ylim(0,8.0)
#    plt.xlabel(r'$x$')
#    plt.ylabel(r'$T_\mathrm{gas}$')
#
#    plt.figure(ifig)
#    ifig += 1
#    line = lineout(data,field='f',logy=False,ycut=0.0,zcut=0.0)
#    plt.setp(line,color=cm.to_rgba(step),antialiased=True)
#    plt.xlim(0,40)
#    plt.ylim(0,1)
#    plt.xlabel(r'$z$')
#    plt.ylabel(r'$f$')


plt.show()
print 'done'
