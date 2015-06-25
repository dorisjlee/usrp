import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,Normalize

def pcolor(data,field='density',cmap='cubehelix',tight=True,transpose=False,logz=False,**kwargs):
    if   any(x in kwargs.keys() for x in ['icut','jcut','kcut']):
        X,Y,Z = data.get_slice_ijk(field,transpose,**kwargs)
    elif any(x in kwargs.keys() for x in ['xcut','ycut','zcut']):
        X,Y,Z = data.get_slice_xyz(field,transpose,**kwargs)
    else:
        raise ValueError('Exactly one of icut,jcut,kcut or xcut,ycut,zcut must be defined!')
    if 'clim' not in kwargs.keys():
        kwargs['clim'] = [Z.min(),Z.max()]
    if logz:
        norm = LogNorm(vmin=kwargs['clim'][0],vmax=kwargs['clim'][1])
    else:
        norm = Normalize(vmin=kwargs['clim'][0],vmax=kwargs['clim'][1])
    figh = plt.pcolor(X,Y,Z,cmap=cmap,norm=norm)
    if tight:
        ax = plt.gca()
        ax.set_aspect('equal')
        ax.autoscale(tight=True)
    return figh
        
def lineout(data,field='density',transpose=False,logx=False,logy=False,**kwargs):
    if   any(x in kwargs.keys() for x in ['icut','jcut','kcut']):
        x,y = data.get_lineout_ijk(field,**kwargs)
    elif any(x in kwargs.keys() for x in ['xcut','ycut','zcut']):
        x,y = data.get_lineout_xyz(field,**kwargs)
    else:
        raise ValueError('Exactly two of icut,jcut,kcut or xcut,ycut,zcut must be defined!')
    figh = plt.plot(x,y)
    if logx:
        plt.xscale('log')
    if logy:
        plt.yscale('log')
    return figh