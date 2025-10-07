import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp

class Geometry():
    pass

def get_next_data(f,ndat,ndat_per_line=5,ncol_per_dat=16):
    data = []
    for idat in range(0,ndat):
        data.append(float(f.read(ncol_per_dat)))
        if ((idat+1)%ndat_per_line == 0) or (idat == ndat-1):
            dummy=f.read(1)

    return np.array(data)

# Needed from eqdsk file:
# - rgrid
# - zgrid
# - psirz
# - ssimag
# - rmaxis
# - zmaxis
# - ssibry
# - rmid
# - B0, R0
def read_geqdsk(gfilename,plot=False):
    g = Geometry()

    f = open(gfilename,"r")
    info=f.read(48)
    dummy=f.read(4)
    mw=int(f.read(4))
    mh=int(f.read(4))
    dummy=f.read(1)

    rdim=float(f.read(16))
    zdim=float(f.read(16))
    rzero=float(f.read(16))
    rmin=float(f.read(16))
    zmid=float(f.read(16))
    dummy=f.read(1)

    rmaxis=float(f.read(16))
    zmaxis=float(f.read(16))
    ssimag=float(f.read(16))
    ssibry=float(f.read(16))
    bcenter=float(f.read(16))
    dummy=f.read(1)

    current=float(f.read(16))
    simag=float(f.read(16))
    dummy=float(f.read(16))
    rmaxis=float(f.read(16))
    dummy=float(f.read(16))
    dummy=f.read(1)

    zmaxis=float(f.read(16))
    dummy=float(f.read(16))
    sibry=float(f.read(16))
    while dummy != '\n':
        dummy=f.read(1)

    fpol=get_next_data(f,mw)
    pres=get_next_data(f,mw)
    workk=get_next_data(f,mw)
    workk=get_next_data(f,mw)
    pprime=-workk

    psirz=get_next_data(f,mw*mh)
#    psirz=psirz.reshape([mw,mh]).transpose()
    psirz=psirz.reshape([mh,mw])

    q=get_next_data(f,mw)
    
    dummy="0"
    nsep=int(f.read(5))
    nlim=int(f.read(5))
    while dummy != '\n':
        dummy=f.read(1)

    sep = get_next_data(f,2*nsep)
    sep = sep.reshape([nsep,2])
    lim = get_next_data(f,2*nlim)
    lim = lim.reshape([nlim,2])

    f.close()

    if ssibry < ssimag:
        psirz = -psirz
        temp = ssibry
        ssibry = ssimag
        ssimag = temp


    psi1 = np.array(range(0,mw))/(mw-1)

    dr=rdim/(mw-1)
    dz=zdim/(mh-1)
    zmin=zmid-0.5*zdim

    rgrid=rmin+np.array(range(0,mw))*dr
    zgrid=zmin+np.array(range(0,mh))*dz

    if plot:
        plt.contour(rgrid,zgrid,psirz.transpose(),levels=150,linewidths=0.4) 
        plt.plot(lim[:,0],lim[:,1])
        plt.plot(sep[:,0],sep[:,1])
        plt.gca().set_aspect("equal")
        plt.xlabel("R")
        plt.ylabel("Z")
        plt.tight_layout()
        plt.savefig("psicontour.pdf")
        plt.clf()

    nrmid = 1000
    rmid_max = max(lim[:,0])
    drmid = (rmid_max - rmaxis)/(nrmid-1)
    rmid = np.array(rmaxis + drmid*range(0,nrmid))
    zmid = np.array([zmaxis]*nrmid)

    psi_interp = interp.RectBivariateSpline(rgrid,zgrid,np.transpose(psirz))
    psimid = psi_interp(rmid,zmaxis)

    R0=rmaxis
    B0=fpol[0]/R0


    g.rgrid = np.array(rgrid)
    g.zgrid = np.array(zgrid)
    g.psirz = np.array(psirz)
    g.ssimag = ssimag
    g.ssibry = ssibry
    g.rmid = np.array(rmid)
    g.psimid = np.array(psimid)
    g.rmaxis = rmaxis
    g.zmaxis = zmaxis
    g.B0 = B0
    g.R0 = R0
    g.lim = lim
    g.sep = sep
    g.qpsi = q
    g.current = float(current)

    return g
