import h5py
import matplotlib.pyplot as plt
import numpy as np
import hdfdict
import sys
from pathlib import Path
import glob


def mask_outside_polygon(poly_verts, ax=None):
    """
    Plots a mask on the specified axis ("ax", defaults to plt.gca()) such that
    all areas outside of the polygon specified by "poly_verts" are masked.  

    "poly_verts" must be a list of tuples of the verticies in the polygon in
    counter-clockwise order.

    Returns the matplotlib.patches.PathPatch instance plotted on the figure.
    """
    import matplotlib.patches as mpatches
    import matplotlib.path as mpath

    if ax is None:
        ax = plt.gca()

    # Get current plot limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    # Verticies of the plot boundaries in clockwise order
    bound_verts = [(xlim[0], ylim[0]), (xlim[0], ylim[1]), 
                   (xlim[1], ylim[1]), (xlim[1], ylim[0]), 
                   (xlim[0], ylim[0])]

    # A series of codes (1 and 2) to tell matplotlib whether to draw a line or 
    # move the "pen" (So that there's no connecting line)
    bound_codes = [mpath.Path.MOVETO] + (len(bound_verts) - 1) * [mpath.Path.LINETO]
    poly_codes = [mpath.Path.MOVETO] + (len(poly_verts) - 1) * [mpath.Path.LINETO]

    # Plot the masking patch
    path = mpath.Path(bound_verts + poly_verts, bound_codes + poly_codes)
    patch = mpatches.PathPatch(path, facecolor='white', edgecolor='none')
    patch = ax.add_patch(patch)

    # Reset the plot limits to their original extents
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    return patch




def get_sol_d(filename = 'plasma_1MW.h5', spec='/spec1/',name = 'n'):
    field = spec+name
    iphi = 1
    is2D = True

    # colorstyle
    colormap = 'jet'
    #plt.style.use('dark_background')

    fluctuations = 0 #0 : pp plane, 1: fluctuations, 2: phi average
    logplot = False
    logplotdecades = 3



    # tune colorbar range (auto if cmin=cmax)
    cmin=0
    cmax=0

    # reference parameters
    ref=h5py.File('refParam_raptorX.h5','r')
    n0=ref['/n0'][...]
    T0=ref['/T0'][...]
    c0=ref['/c0'][...]
    W0=ref['/W0'][...]
    rho0=ref['/rho0'][...]
    R0=ref['/R0'][...]
    a0=ref['/a0'][...]
    ref.close()
    if(field[-2:]=='/n'):
        X0 = n0
    elif(field[-2:]=='/T'):
        X0 = T0
    elif(field[-2:]=='/G'):
        X0 = n0*c0
    elif(field[-2:]=='/k'):
        X0 = c0*c0
    elif(field[-2:]=='/D'):
        X0 = c0*rho0
    elif(field[-2:]=='/v'):
        X0 = c0
    elif(field[-2:]=='HI'):
        X0 = T0
    elif(field[-2:]=='/W'):
        X0 = W0
    elif(field[-2:]=='Ee'):
        X0 = n0*T0*1.6e-19
    elif(field[-2:]=='Ei'):
        X0 = n0*T0*1.6e-19
    else:
        X0 = 1

    if is2D:
        NG=0
    else:
        NG=2

    # wall and psi
    mesh=h5py.File('mesh.h5','r')
    Rwall=mesh['/wall/R'][...]
    Zwall=mesh['/wall/Z'][...]
    r2D=mesh['/config/r'][...]
    z2D=mesh['/config/z'][...]
    psi2D=mesh['/config/psi'][...]
    psisep=mesh['/config/psisep1'][...]
    #rescale
    R0init=mesh['/R0'][...]
    a0init=mesh['/a0'][...]
    Rwall = (Rwall-R0init)/a0init*a0+R0
    Zwall = (Zwall)/a0init*a0
    r2D = (r2D-R0init)/a0init*a0+R0
    z2D = (z2D)/a0init*a0
    # restrict to points within wall limits
    S=np.shape(r2D)
    for i in range(0,S[0]):
        if(r2D[i,0]>np.min(Rwall)):
            break
    iRmin=i
    for i in range(0,S[0]):
        if(r2D[i,0]>np.max(Rwall)):
            break
    iRmax=i
    for i in range(0,S[1]):
        if(z2D[0,i]>np.min(Zwall)):
            break
    iZmin=i
    for i in range(0,S[0]):
        if(z2D[0,i]>np.max(Zwall)):
            break
    iZmax=i
    r2D=r2D[iRmin:iRmax,iZmin:iZmax]
    z2D=z2D[iRmin:iRmax,iZmin:iZmax]
    psi2D=psi2D[iRmin:iRmax,iZmin:iZmax]
    mesh.close()


    mesh=h5py.File('mesh_raptorX.h5','r')
    metric=h5py.File('metric_raptorX.h5','r')
    data=h5py.File(filename,'r')

    Nzones = int(mesh['/NZones'][...])

    #fig = plt.figure()

    # get min max 
    Xmin=1e100
    Xmax=-1e100
    for izone in range(0,Nzones):
        if(field[-2:]=='/v'):
            G = data['/zone'+str(izone+1)+'/spec1/G'][...]
            n = data['/zone'+str(izone+1)+'/spec1/n'][...]
            X=G/n
        elif(field[-2:]=='/M'):
            G = data['/zone'+str(izone+1)+'/spec1/G'][...]
            n = data['/zone'+str(izone+1)+'/spec1/n'][...]
            Te = data['/zone'+str(izone+1)+'/spec0/T'][...]
            Ti = data['/zone'+str(izone+1)+'/spec0/T'][...]
            mi = data['/zone'+str(izone+1)+'/spec1/mass'][...]
            cs=np.sqrt((Te+Ti)/mi)
            X=G/(n*cs)
        elif(field[-2:]=='/Mtheta'):
            G = data['/zone'+str(izone+1)+'/spec1/G'][...]
            n = data['/zone'+str(izone+1)+'/spec1/n'][...]
            Te = data['/zone'+str(izone+1)+'/spec0/T'][...]
            Ti = data['/zone'+str(izone+1)+'/spec0/T'][...]
            uEtheta = data['/zone'+str(izone+1)+'/ExB/theta'][...]
            uBtheta = data['/zone'+str(izone+1)+'/spec1/uGradB/theta'][...]
            b2 = metric['/zone'+str(izone+1)+'/b2'][...]
            mi = data['/zone'+str(izone+1)+'/spec1/mass'][...]
            cs=np.sqrt((Te+Ti)/mi)
            X=(G/n*b2+uEtheta+uBtheta)/(cs*b2)
        elif(field[-2:]=='Ee'):
            n = data['/zone'+str(izone+1)+'/spec1/n'][...]
            Te = data['/zone'+str(izone+1)+'/spec0/T'][...]
            X=3/2*n*Te
        elif(field[-2:]=='Ei'):
            n = data['/zone'+str(izone+1)+'/spec1/n'][...]
            G = data['/zone'+str(izone+1)+'/spec1/G'][...]
            Ti = data['/zone'+str(izone+1)+'/spec1/T'][...]
            mi = data['/zone'+str(izone+1)+'/spec1/mass'][...]
            X=3/2*n*Ti+1/2*mi*n*(G/n)**2
        elif(field[-2:]=='/D'):
            k = data['/zone'+str(izone+1)+'/k'][...]
            Te = data['/zone'+str(izone+1)+'/spec0/T'][...]
            Ti = data['/zone'+str(izone+1)+'/spec1/T'][...]
            mi = data['/zone'+str(izone+1)+'/spec1/mass'][...]
            cs=np.sqrt((Te+Ti)/mi)
            R = mesh['/zone'+str(izone+1)+'/Rcells'][...]
            X = k*R/cs
        else:
            X = data['/zone'+str(izone+1)+field][...]
        try:
            if(field[-6:]=='/theta'):
                g22=metric['/zone'+str(izone+1)+'/g22'][...]
                X=X/np.sqrt(g22)*c0
            if(field[-4:]=='/psi'):
                g11=metric['/zone'+str(izone+1)+'/g11'][...]
                X=X/np.sqrt(g11)*c0
        except:
            print("catched")

        if(fluctuations==1):
            X = X[2:-2,2:-2,iphi-1+NG] - np.mean(X[2:-2,2:-2,NG:-NG],2)
        elif(fluctuations==2):
            X = np.mean(X[2:-2,2:-2,NG:-NG],2)*X0
        else:
            X = X[2:-2,2:-2,iphi-1+NG]
        Xmin=min(Xmin,np.min(X))
        Xmax=max(Xmax,np.max(X))

    if(cmin==cmax):
        if(fluctuations==1):
            cmin=min(Xmin,-Xmax)
            cmax=max(Xmax,-Xmin)
        else:
            cmin=Xmin
            cmax=Xmax

    cmin = cmin*X0
    cmax = cmax*X0
    if(field[-2:]=='/M'):
        cmin=-1
        cmax=1

    if(logplot):
        cmax = int(np.log10(cmax))+1
        cmin = cmax - logplotdecades

    #sm = plt.cm.ScalarMappable(cmap=colormap, norm=plt.Normalize(vmin=cmin, vmax=cmax))
    # plot
    r_tot = []
    z_tot = []
    n_tot = []

    for izone in range(0,Nzones):
        Rc = mesh['/zone'+str(izone+1)+'/Rcorners'][...]
        Zc = mesh['/zone'+str(izone+1)+'/Zcorners'][...]
        Rc = Rc[:,:,0]*rho0
        Zc = Zc[:,:,0]*rho0
        if(field[-2:]=='/v'):
            G = data['/zone'+str(izone+1)+'/spec1/G'][...]
            n = data['/zone'+str(izone+1)+'/spec1/n'][...]
            X=G/n
        elif(field[-2:]=='/M'):
            G = data['/zone'+str(izone+1)+'/spec1/G'][...]
            n = data['/zone'+str(izone+1)+'/spec1/n'][...]
            Te = data['/zone'+str(izone+1)+'/spec0/T'][...]
            Ti = data['/zone'+str(izone+1)+'/spec0/T'][...]
            mi = data['/zone'+str(izone+1)+'/spec1/mass'][...]
            cs=np.sqrt((Te+Ti)/mi)
            X=G/(n*cs)
        elif(field[-2:]=='/Mtheta'):
            G = data['/zone'+str(izone+1)+'/spec1/G'][...]
            n = data['/zone'+str(izone+1)+'/spec1/n'][...]
            Te = data['/zone'+str(izone+1)+'/spec0/T'][...]
            Ti = data['/zone'+str(izone+1)+'/spec0/T'][...]
            uEtheta = data['/zone'+str(izone+1)+'/ExB/theta'][...]
            uBtheta = data['/zone'+str(izone+1)+'/spec1/uGradB/theta'][...]
            b2 = metric['/zone'+str(izone+1)+'/b2'][...]
            mi = data['/zone'+str(izone+1)+'/spec1/mass'][...]
            cs=np.sqrt((Te+Ti)/mi)
            X=(G/n*b2+uEtheta+uBtheta)/(cs*b2)
        elif(field[-2:]=='Ee'):
            n = data['/zone'+str(izone+1)+'/spec1/n'][...]
            Te = data['/zone'+str(izone+1)+'/spec0/T'][...]
            X=3/2*n*Te
        elif(field[-2:]=='Ei'):
            n = data['/zone'+str(izone+1)+'/spec1/n'][...]
            G = data['/zone'+str(izone+1)+'/spec1/G'][...]
            Ti = data['/zone'+str(izone+1)+'/spec0/T'][...]
            mi = data['/zone'+str(izone+1)+'/spec1/mass'][...]
            X=3/2*n*Ti+1/2*mi*n*(G/n)**2
        elif(field[-2:]=='/D'):
            k = data['/zone'+str(izone+1)+'/k'][...]
            Te = data['/zone'+str(izone+1)+'/spec0/T'][...]
            Ti = data['/zone'+str(izone+1)+'/spec1/T'][...]
            mi = data['/zone'+str(izone+1)+'/spec1/mass'][...]
            cs=np.sqrt((Te+Ti)/mi)
            R = mesh['/zone'+str(izone+1)+'/Rcells'][...]
            X = k*R/cs
        else:
            X = data['/zone'+str(izone+1)+field][...]

        try:
            if(field[-6:]=='/theta'):
                g22=metric['/zone'+str(izone+1)+'/g22'][...]
                X=X/np.sqrt(g22)*c0
            if(field[-4:]=='/psi'):
                g11=metric['/zone'+str(izone+1)+'/g11'][...]
                X=X/np.sqrt(g11)*c0
        except:
            print("catched")

        if(fluctuations==1):
            X = X[2:-2,2:-2,iphi-1+NG] - np.mean(X[2:-2,2:-2,NG:-NG],2)
            X = X*X0
        elif(fluctuations==2):
            X = np.mean(X[2:-2,2:-2,NG:-NG],2)*X0
        else:
            X = X[2:-2,2:-2,iphi-1+NG]*X0        
        if(logplot):
            X = np.log10(X)
        #plt.pcolor(Rc,Zc,X,vmin=cmin,vmax=cmax,cmap=colormap)
        r_tot.append(Rc)
        z_tot.append(Zc)
        #n_tot.append(X)












    #field = '/spec0/T'
    iphi = 1
    is2D = True

    # colorstyle
    colormap = 'jet'
    #plt.style.use('dark_background')

    fluctuations = 0 #0 : pp plane, 1: fluctuations, 2: phi average
    logplot = False
    logplotdecades = 3



    # tune colorbar range (auto if cmin=cmax)
    cmin=0
    cmax=0

    # reference parameters
    ref=h5py.File('refParam_raptorX.h5','r')
    n0=ref['/n0'][...]
    T0=ref['/T0'][...]
    c0=ref['/c0'][...]
    W0=ref['/W0'][...]
    rho0=ref['/rho0'][...]
    R0=ref['/R0'][...]
    a0=ref['/a0'][...]

    if(field[-2:]=='/n'):
        X0 = n0
    elif(field[-2:]=='/T'):
        X0 = T0
    elif(field[-2:]=='/G'):
        X0 = n0*c0
    elif(field[-2:]=='/k'):
        X0 = c0*c0
    elif(field[-2:]=='/D'):
        X0 = c0*rho0
    elif(field[-2:]=='/v'):
        X0 = c0
    elif(field[-2:]=='HI'):
        X0 = T0
    elif(field[-2:]=='/W'):
        X0 = W0
    elif(field[-2:]=='Ee'):
        X0 = n0*T0*1.6e-19
    elif(field[-2:]=='Ei'):
        X0 = n0*T0*1.6e-19
    else:
        X0 = 1

    if is2D:
        NG=0
    else:
        NG=2

    # wall and psi
    mesh=h5py.File('mesh.h5','r')
    Rwall=mesh['/wall/R'][...]
    Zwall=mesh['/wall/Z'][...]
    r2D=mesh['/config/r'][...]
    z2D=mesh['/config/z'][...]
    psi2D=mesh['/config/psi'][...]
    psisep=mesh['/config/psisep1'][...]
    #rescale
    R0init=mesh['/R0'][...]
    a0init=mesh['/a0'][...]
    Rwall = (Rwall-R0init)/a0init*a0+R0
    Zwall = (Zwall)/a0init*a0
    r2D = (r2D-R0init)/a0init*a0+R0
    z2D = (z2D)/a0init*a0
    # restrict to points within wall limits
    S=np.shape(r2D)
    for  i in range(0,S[0]):
        if(r2D[i,0]>np.min(Rwall)):
            break
    iRmin=i
    for i in range(0,S[0]):
        if(r2D[i,0]>np.max(Rwall)):
            break
    iRmax=i
    for i in range(0,S[1]):
        if(z2D[0,i]>np.min(Zwall)):
            break
    iZmin=i
    for i in range(0,S[0]):
        if(z2D[0,i]>np.max(Zwall)):
            break
    iZmax=i
    r2D=r2D[iRmin:iRmax,iZmin:iZmax]
    z2D=z2D[iRmin:iRmax,iZmin:iZmax]
    psi2D=psi2D[iRmin:iRmax,iZmin:iZmax]
    mesh.close()


    mesh=h5py.File('mesh_raptorX.h5','r')
    metric=h5py.File('metric_raptorX.h5','r')
    data=h5py.File(filename,'r')

    Nzones = int(mesh['/NZones'][...])

    #fig = plt.figure()

    # get min max 
    Xmin=1e100
    Xmax=-1e100
    for izone in range(0,Nzones):
        if(field[-2:]=='/v'):
            G = data['/zone'+str(izone+1)+'/spec1/G'][...]
            n = data['/zone'+str(izone+1)+'/spec1/n'][...]
            X=G/n
        elif(field[-2:]=='/M'):
            G = data['/zone'+str(izone+1)+'/spec1/G'][...]
            n = data['/zone'+str(izone+1)+'/spec1/n'][...]
            Te = data['/zone'+str(izone+1)+'/spec0/T'][...]
            Ti = data['/zone'+str(izone+1)+'/spec0/T'][...]
            mi = data['/zone'+str(izone+1)+'/spec1/mass'][...]
            cs=np.sqrt((Te+Ti)/mi)
            X=G/(n*cs)
        elif(field[-2:]=='/Mtheta'):
            G = data['/zone'+str(izone+1)+'/spec1/G'][...]
            n = data['/zone'+str(izone+1)+'/spec1/n'][...]
            Te = data['/zone'+str(izone+1)+'/spec0/T'][...]
            Ti = data['/zone'+str(izone+1)+'/spec0/T'][...]
            uEtheta = data['/zone'+str(izone+1)+'/ExB/theta'][...]
            uBtheta = data['/zone'+str(izone+1)+'/spec1/uGradB/theta'][...]
            b2 = metric['/zone'+str(izone+1)+'/b2'][...]
            mi = data['/zone'+str(izone+1)+'/spec1/mass'][...]
            cs=np.sqrt((Te+Ti)/mi)
            X=(G/n*b2+uEtheta+uBtheta)/(cs*b2)
        elif(field[-2:]=='Ee'):
            n = data['/zone'+str(izone+1)+'/spec1/n'][...]
            Te = data['/zone'+str(izone+1)+'/spec0/T'][...]
            X=3/2*n*Te
        elif(field[-2:]=='Ei'):
            n = data['/zone'+str(izone+1)+'/spec1/n'][...]
            G = data['/zone'+str(izone+1)+'/spec1/G'][...]
            Ti = data['/zone'+str(izone+1)+'/spec1/T'][...]
            mi = data['/zone'+str(izone+1)+'/spec1/mass'][...]
            X=3/2*n*Ti+1/2*mi*n*(G/n)**2
        elif(field[-2:]=='/D'):
            k = data['/zone'+str(izone+1)+'/k'][...]
            Te = data['/zone'+str(izone+1)+'/spec0/T'][...]
            Ti = data['/zone'+str(izone+1)+'/spec1/T'][...]
            mi = data['/zone'+str(izone+1)+'/spec1/mass'][...]
            cs=np.sqrt((Te+Ti)/mi)
            R = mesh['/zone'+str(izone+1)+'/Rcells'][...]
            X = k*R/cs
        else:
            X = data['/zone'+str(izone+1)+field][...]
        try:
            if(field[-6:]=='/theta'):
                g22=metric['/zone'+str(izone+1)+'/g22'][...]
                X=X/np.sqrt(g22)*c0
            if(field[-4:]=='/psi'):
                g11=metric['/zone'+str(izone+1)+'/g11'][...]
                X=X/np.sqrt(g11)*c0
        except:
            print("catched")

        if(fluctuations==1):
            X = X[2:-2,2:-2,iphi-1+NG] - np.mean(X[2:-2,2:-2,NG:-NG],2)
        elif(fluctuations==2):
            X = np.mean(X[2:-2,2:-2,NG:-NG],2)*X0
        else:
            X = X[2:-2,2:-2,iphi-1+NG]
        Xmin=min(Xmin,np.min(X))
        Xmax=max(Xmax,np.max(X))

    if(cmin==cmax):
        if(fluctuations==1):
            cmin=min(Xmin,-Xmax)
            cmax=max(Xmax,-Xmin)
        else:
            cmin=Xmin
            cmax=Xmax

    cmin = cmin*X0
    cmax = cmax*X0
    if(field[-2:]=='/M'):
        cmin=-1
        cmax=1

    if(logplot):
        cmax = int(np.log10(cmax))+1
        cmin = cmax - logplotdecades

    #sm = plt.cm.ScalarMappable(cmap=colormap, norm=plt.Normalize(vmin=cmin, vmax=cmax))
    # plot
    r_tot = []
    z_tot = []
    t_tot = []
    bphi_tot = []
    br_tot = []
    bz_tot = []
    b_tot = []
    parr_flow_tot = []
    temp_tot = []
    for izone in range(0,Nzones):
        Rc = mesh['/zone'+str(izone+1)+'/Rcorners'][...]
        Zc = mesh['/zone'+str(izone+1)+'/Zcorners'][...]
        
        br   = mesh['/zone'+str(izone+1)+'/Br'][...]
        bphi = mesh['/zone'+str(izone+1)+'/Bphi'][...]
        bz   = mesh['/zone'+str(izone+1)+'/Bz'][...]
        b    = mesh['/zone'+str(izone+1)+'/B'][...]


        
        Rc = Rc[:,:,0]*rho0
        Zc = Zc[:,:,0]*rho0
        if(field[-2:]=='/v'):
            G = data['/zone'+str(izone+1)+'/spec1/G'][...]
            n = data['/zone'+str(izone+1)+'/spec1/n'][...]
            X=G/n
        elif(field[-2:]=='/M'):
            G = data['/zone'+str(izone+1)+'/spec1/G'][...]
            n = data['/zone'+str(izone+1)+'/spec1/n'][...]
            Te = data['/zone'+str(izone+1)+'/spec0/T'][...]
            Ti = data['/zone'+str(izone+1)+'/spec0/T'][...]
            mi = data['/zone'+str(izone+1)+'/spec1/mass'][...]
            cs=np.sqrt((Te+Ti)/mi)
            X=G/(n*cs)
        elif(field[-2:]=='/Mtheta'):
            G = data['/zone'+str(izone+1)+'/spec1/G'][...]
            n = data['/zone'+str(izone+1)+'/spec1/n'][...]
            Te = data['/zone'+str(izone+1)+'/spec0/T'][...]
            Ti = data['/zone'+str(izone+1)+'/spec0/T'][...]
            uEtheta = data['/zone'+str(izone+1)+'/ExB/theta'][...]
            uBtheta = data['/zone'+str(izone+1)+'/spec1/uGradB/theta'][...]
            b2 = metric['/zone'+str(izone+1)+'/b2'][...]
            mi = data['/zone'+str(izone+1)+'/spec1/mass'][...]
            cs=np.sqrt((Te+Ti)/mi)
            X=(G/n*b2+uEtheta+uBtheta)/(cs*b2)
        elif(field[-2:]=='Ee'):
            n = data['/zone'+str(izone+1)+'/spec1/n'][...]
            Te = data['/zone'+str(izone+1)+'/spec0/T'][...]
            X=3/2*n*Te
        elif(field[-2:]=='Ei'):
            n = data['/zone'+str(izone+1)+'/spec1/n'][...]
            G = data['/zone'+str(izone+1)+'/spec1/G'][...]
            Ti = data['/zone'+str(izone+1)+'/spec0/T'][...]
            mi = data['/zone'+str(izone+1)+'/spec1/mass'][...]
            X=3/2*n*Ti+1/2*mi*n*(G/n)**2
        elif(field[-2:]=='/D'):
            k = data['/zone'+str(izone+1)+'/k'][...]
            Te = data['/zone'+str(izone+1)+'/spec0/T'][...]
            Ti = data['/zone'+str(izone+1)+'/spec1/T'][...]
            mi = data['/zone'+str(izone+1)+'/spec1/mass'][...]
            cs=np.sqrt((Te+Ti)/mi)
            R = mesh['/zone'+str(izone+1)+'/Rcells'][...]
            X = k*R/cs
        else:
            X = data['/zone'+str(izone+1)+field][...]
            parr_flow = data['/zone'+str(izone+1)+spec+'G'][...]
            temp = data['/zone'+str(izone+1)+spec+'T'][...]
            
        try:
            if(field[-6:]=='/theta'):
                g22=metric['/zone'+str(izone+1)+'/g22'][...]
                X=X/np.sqrt(g22)*c0
            if(field[-4:]=='/psi'):
                g11=metric['/zone'+str(izone+1)+'/g11'][...]
                X=X/np.sqrt(g11)*c0
        except:
            print("catched")

        if(fluctuations==1):
            X = X[2:-2,2:-2,iphi-1+NG] - np.mean(X[2:-2,2:-2,NG:-NG],2)
            X = X*X0
        elif(fluctuations==2):
            X = np.mean(X[2:-2,2:-2,NG:-NG],2)*X0
        else:
            X = X[2:-2,2:-2,iphi-1+NG]*X0
            bphi = bphi[2:-2,2:-2,iphi-1+NG]*ref['B0'][...]
            br = br[2:-2,2:-2,iphi-1+NG]*ref['B0'][...]
            bz = bz[2:-2,2:-2,iphi-1+NG]*ref['B0'][...]
            b = b[2:-2,2:-2,iphi-1+NG]*ref['B0'][...]
            parr_flow = parr_flow[2:-2,2:-2,iphi-1+NG]*ref['n0'][...]*ref['c0'][...]
            temp = temp[2:-2,2:-2,iphi-1+NG]*ref['T0'][...]
            
        if(logplot):
            X = np.log10(X)
        #plt.pcolor(Rc,Zc,X)#,vmin=cmin,vmax=cmax,cmap=colormap)
        r_tot.append(Rc)
        z_tot.append(Zc)
        t_tot.append(X)
        bphi_tot.append(bphi)
        br_tot.append(br)
        bz_tot.append(bz)
        b_tot.append(b)
        parr_flow_tot.append(parr_flow)
        temp_tot.append(temp)
    return r_tot, z_tot, t_tot, Rwall, Zwall,r2D,z2D,psi2D,psisep,bphi_tot,br_tot,bz_tot,b_tot,parr_flow_tot,temp_tot



def remesh_soledge_cherab_data(t_tot, name='out',tot=np.zeros((146,165)),save=False):
    kk = 0
    jj = 0
    for i in range(0,16):#16 zones for now should not hard code this
        for k in range(0,np.shape(t_tot[i])[0]):
            for j in range(0,np.shape(t_tot[i])[1]):
                if(jj==np.shape(tot)[1]):
                    jj = 0
                    kk = kk +1
                tot[kk,jj] = t_tot[i][k,j]
                jj = jj +1

    if(save):
        np.save(name+'.npy',tot.T)
    else:
        return tot.T

    
            
def remesh_soledge_cherab_grid(r_tot,z_tot, name='', cr_z = np.zeros((146,165,4)), cr_r = np.zeros((146,165,4)),save_npy=False):

    kk = 0
    jj = 0
    for i in range(0,16):
        for k in range(0,np.shape(r_tot[i])[0]-1):
            for j in range(0,np.shape(r_tot[i])[1]-1):

                if(jj==np.shape(cr_r)[1]):
                    jj = 0
                    kk = kk +1
                cr_z[kk,jj,2] = z_tot[i][k,j]
                cr_z[kk,jj:,3] = z_tot[i][k+1,j]
                cr_z[kk,jj,1] = z_tot[i][k,j+1]
                cr_z[kk,jj,0] = z_tot[i][k+1,j+1]

                cr_r[kk,jj,2] = r_tot[i][k,j]
                cr_r[kk,jj:,3] = r_tot[i][k+1,j]
                cr_r[kk,jj,1] = r_tot[i][k,j+1]
                cr_r[kk,jj,0] = r_tot[i][k+1,j+1]

                jj = jj +1

    neigh_i = np.ones_like(cr_r)*-1
    neigh_j = np.copy(neigh_i)
    for i in range(0,146):
        for j in range(0,165):
            for k in range(0,4):
                if(k ==0):
                    f = np.where( (np.abs(cr_r -cr_r[i,j,0]) ==0) & (np.abs(cr_z -cr_z[i,j,0]) == 0))
                    s = np.where( (np.abs(cr_r -cr_r[i,j,1]) ==0) & (np.abs(cr_z -cr_z[i,j,1]) == 0))
                    ss = np.vstack(s)
                    ff = np.vstack(f)

                    for ii in range(0,np.shape(ss)[1]):
                        if(np.all(ss[:,ii] == np.array([i,j,1]))):
                            ss = np.delete(ss,ii,axis=1)

                            break

                    for ii in range(0,np.shape(ff)[1]):
                        if(np.all(ff[:,ii] ==np.array([i,j,0]))):
                            ff = np.delete(ff,ii,axis=1)

                            break
                    ss = ss[0:2,:]
                    ff = ff[0:2,:]

                    for ij in range(0,np.shape(ff)[1]):
                        for ji in range(0,np.shape(ss)[1]):
                            if(np.all(ff[:,ij] == ss[:,ji])):
                                neigh_i[i,j,k] = ff[0,ij]
                                neigh_j[i,j,k] = ff[1,ij]



                if(k ==1):
                    f = np.where( (np.abs(cr_r -cr_r[i,j,1]) ==0) & (np.abs(cr_z -cr_z[i,j,1]) == 0))
                    s = np.where( (np.abs(cr_r -cr_r[i,j,2]) ==0) & (np.abs(cr_z -cr_z[i,j,2]) == 0))
                    ss = np.vstack(s)
                    ff = np.vstack(f)

                    for ii in range(0,np.shape(ss)[1]):
                        if(np.all(ss[:,ii] == np.array([i,j,2]))):
                            ss = np.delete(ss,ii,axis=1)

                            break

                    for ii in range(0,np.shape(ff)[1]):
                        if(np.all(ff[:,ii] == np.array([i,j,1]))):
                            ff = np.delete(ff,ii,axis=1)

                            break
                    ss = ss[0:2,:]
                    ff = ff[0:2,:]

                    for ij in range(0,np.shape(ff)[1]):
                        for ji in range(0,np.shape(ss)[1]):
                            if(np.all(ff[:,ij] == ss[:,ji])):
                                neigh_i[i,j,k] = ff[0,ij]
                                neigh_j[i,j,k] = ff[1,ij]

                if(k ==2):
                    f = np.where( (np.abs(cr_r -cr_r[i,j,2]) ==0) & (np.abs(cr_z -cr_z[i,j,2]) == 0))
                    s = np.where( (np.abs(cr_r -cr_r[i,j,3]) ==0) & (np.abs(cr_z -cr_z[i,j,3]) == 0))
                    ss = np.vstack(s)
                    ff = np.vstack(f)

                    for ii in range(0,np.shape(ss)[1]):
                        if(np.all(ss[:,ii] == np.array([i,j,3]))):
                            ss = np.delete(ss,ii,axis=1)

                            break

                    for ii in range(0,np.shape(ff)[1]):
                        if(np.all(ff[:,ii] == np.array([i,j,2]))):
                            ff = np.delete(ff,ii,axis=1)

                            break
                    ss = ss[0:2,:]
                    ff = ff[0:2,:]

                    for ij in range(0,np.shape(ff)[1]):
                        for ji in range(0,np.shape(ss)[1]):
                            if(np.all(ff[:,ij] == ss[:,ji])):
                                neigh_i[i,j,k] = ff[0,ij]
                                neigh_j[i,j,k] = ff[1,ij]

                if(k ==3):
                    f = np.where( (np.abs(cr_r -cr_r[i,j,3]) ==0) & (np.abs(cr_z -cr_z[i,j,3]) == 0))
                    s = np.where( (np.abs(cr_r -cr_r[i,j,0]) ==0) & (np.abs(cr_z -cr_z[i,j,0]) == 0))
                    ss = np.vstack(s)
                    ff = np.vstack(f)

                    for ii in range(0,np.shape(ss)[1]):
                        if(np.all(ss[:,ii] == np.array([i,j,0]))):
                            ss = np.delete(ss,ii,axis=1)

                            break

                    for ii in range(0,np.shape(ff)[1]):
                        if(np.all(ff[:,ii] == np.array([i,j,3]))):
                            ff = np.delete(ff,ii,axis=1)

                            break
                    ss = ss[0:2,:]
                    ff = ff[0:2,:]

                    for ij in range(0,np.shape(ff)[1]):
                        for ji in range(0,np.shape(ss)[1]):
                            if(np.all(ff[:,ij] == ss[:,ji])):
                                neigh_i[i,j,k] = ff[0,ij]
                                neigh_j[i,j,k] = ff[1,ij]


    if(save_npy):
        np.save(name +'_r.npy', np.transpose(cr_r, (2,1,0)))
        np.save(name +'_z.npy', np.transpose(cr_z, (2,1,0)))
        np.save(name +'_neigh_i.npy', np.transpose(neigh_i, (2,1,0)))
        np.save(name +'_neigh_j.npy', np.transpose(neigh_j, (2,1,0)))
    else:
        return(cr_r,cr_z,neigh_i,neigh_j)



    
fields = np.array(['/spec0/','/spec1/','/spec2/','/spec3/','/spec4/','/spec5/','/spec6/','/spec7/','/spec8/','/spec9/'])
names = np.array(['n_e','n_i','o_1','o_2','o_3','o_4','o_5','o_6','o_7','o_8'])





name = 'plasma_1MW.h5'

sav = {}
for i in range(0,len(fields)):
    print(i)
    r_tot, z_tot, o_tot,Rwall, Zwall,r2D,z2D,psi2D,psisep,bphi,br,bz,b,parr_flow,temp = get_sol_d(filename=name, spec=fields[i],name='n')
    if(i ==0):
        cr_r,cr_z,neigh_i,neigh_j = remesh_soledge_cherab_grid(r_tot,z_tot, name='', cr_z = np.zeros((146,165,4)), cr_r = np.zeros((146,165,4)),save_npy=False)
    
        sav['r_wall_points'] = Rwall
        sav['z_wall_points'] = Zwall
        sav['bfield'] ={}
        sav['bfield']['b_tot'] = np.copy(remesh_soledge_cherab_data(b))
        sav['bfield']['b_phi'] = np.copy(remesh_soledge_cherab_data(bphi))
        sav['bfield']['b_r'] = np.copy(remesh_soledge_cherab_data(br))
        sav['bfield']['b_z'] = np.copy(remesh_soledge_cherab_data(bz))
        sav['solps_like'] = {}
        sav['solps_like']['r'] = cr_r
        sav['solps_like']['z'] = cr_z
        sav['solps_like']['neighbor_i'] = neigh_i
        sav['solps_like']['neighbor_j'] = neigh_j
        
    sav[names[i]] = {}
    sav[names[i]]['dens'] = np.copy(remesh_soledge_cherab_data(o_tot))
    sav[names[i]]['temp'] = np.copy(remesh_soledge_cherab_data(temp))
    sav[names[i]]['parr_flow'] = np.copy(remesh_soledge_cherab_data(parr_flow))



    
hdfdict.dump(sav,h5py.File('psep1p0mw.h5','w'))

# xv, yv = np.meshgrid(np.arange(sav['n_e']['dens'].shape[0]), np.arange(sav['n_e']['dens'].shape[1]))
# indices = np.stack(( yv.flatten(), xv.flatten()), axis=1)
# rr = (sav['solps_like']['r'].reshape(-1,sav['solps_like']['r'].shape[-1]))
# zz = (sav['solps_like']['z'].reshape(-1,sav['solps_like']['z'].shape[-1]))
# verts = np.dstack((rr,zz))


# fig, ax1 = plt.subplots(1,1,figsize = (3,2.9),dpi=100)
# cm = PolyCollection(verts,linewidth=1,facecolor='none',edgecolor='b')
# ax1.add_collection(cm)
# ax1.set_xlim(1.5,4)
# ax1.set_ylim(-1.3,1.3)









# plt.rc('font',size=6)
# plt.rcParams['font.weight'] = 'semibold'

# cm = PolyCollection(verts,antialiaseds=16,edgecolors='face',linewidth=0.1,norm=matplotlib.colors.LogNorm())

# cm.set_array(sav['o_4']['dens'][indices[:,1],indices[:,0]])
# cm.set_clim(vmin=10**15,vmax=10**17)
# cm.set_edgecolors('face')
# cm.set_cmap('plasma')




# fig, ax1 = plt.subplots(1,1,sharex=True,figsize = (3,2.9),dpi=1200)
# fig.subplots_adjust(bottom=0.113,top=0.97,left=0.25,right=0.97,hspace=0.1,wspace=.35)

# ax1.add_collection(cm)
# plt.gca().set_aspect('equal')

# ax1.set_xlabel('R (m)',weight='semibold',size=6)
# ax1.set_ylabel('Z (m)',weight='semibold',size=6,labelpad=-.5)

# ax1.set_xlim(np.min(sav['r_wall_points'])-.1,np.max(sav['r_wall_points'])+.1)
# ax1.set_ylim(np.min(sav['z_wall_points'])-.1,np.max(sav['z_wall_points'])+.1)
# ax1.plot(sav['r_wall_points'][:],sav['z_wall_points'][:],color='g',lw=.8)



# plt.tight_layout()


# tups = []
# for i in range(0,len(sav['r_wall_points'])):
#     tups.append((sav['r_wall_points'][i],sav['z_wall_points'][i]))

# mask_outside_polygon(tups)

# #cm.autoscale_None()
# cbar = fig.colorbar(cm,ax=ax1)

# #cbar.ax1.set_ytick_labels(ax1,['$10^{15}$', '$10^{16}$', '$10^{17}$'])

# #cbar.formatter = LogFormatterExponent(base=10)

# ax1.text(0.05,.05,r'n$_{O^{1+}}$ (m$^{-3}$)',transform=ax1.transAxes,fontsize=9)

# plt.savefig('/home/69j/a.png',format='png')














# cm = PolyCollection(verts,antialiaseds=16,edgecolors='face',linewidth=0.1)

# cm.set_array(sav['bfield']['b_tot'][indices[:,1],indices[:,0]])
# cm.set_clim(vmin=2.7,vmax=6)
# cm.set_edgecolors('face')
# cm.set_cmap('plasma')




# fig, ax1 = plt.subplots(1,1,sharex=True,figsize = (3,2.9),dpi=1200)
# fig.subplots_adjust(bottom=0.113,top=0.97,left=0.25,right=0.97,hspace=0.1,wspace=.35)

# ax1.add_collection(cm)
# plt.gca().set_aspect('equal')

# ax1.set_xlabel('R (m)',weight='semibold',size=6)
# ax1.set_ylabel('Z (m)',weight='semibold',size=6,labelpad=-.5)

# ax1.set_xlim(np.min(sav['r_wall_points'])-.1,np.max(sav['r_wall_points'])+.1)
# ax1.set_ylim(np.min(sav['z_wall_points'])-.1,np.max(sav['z_wall_points'])+.1)
# ax1.plot(sav['r_wall_points'][:],sav['z_wall_points'][:],color='g',lw=.8)



# plt.tight_layout()


# tups = []
# for i in range(0,len(sav['r_wall_points'])):
#     tups.append((sav['r_wall_points'][i],sav['z_wall_points'][i]))

# mask_outside_polygon(tups)

# #cm.autoscale_None()
# cbar = fig.colorbar(cm,ax=ax1)

# #cbar.ax1.set_ytick_labels(ax1,['$10^{15}$', '$10^{16}$', '$10^{17}$'])

# #cbar.formatter = LogFormatterExponent(base=10)

# ax1.text(0.05,.05,r'B$_{tot}$ (T)',transform=ax1.transAxes,fontsize=9)

# plt.savefig('/home/69j/aa.png',format='png')


