from __future__ import print_function 
import numpy as np
import glob,sys
from math import pi

# phis is the average order parameter for an atom inside bulk solid
# phil is the average order parameter for an atom inside bulk liquid
# vs is the molar volume (per atom) for bulk solid
# These quantities were computed from the MD simulations of the bulk phases
bulkproperties = {
"0.6178": {"phis": 0.9851, "phil": 0.0095, "vs": 1.05865},
"0.60": {"phis": 0.9878, "phil": 0.0094, "vs": 1.05238},
"0.58": {"phis": 0.98996, "phil": 0.0094, "vs": 1.046},
"0.56": {"phis": 0.99159, "phil": 0.00926, "vs": 1.04}
}

def read_frame(filedesc):
    natoms = int(filedesc.readline())
    comment = filedesc.readline()
    cell = np.zeros(3,float)
    cell[:] = comment.split()[0:3]
    
    names = np.zeros(natoms,np.dtype('|S6'))
    q = np.zeros((natoms,3),float)
    orderp = np.zeros(natoms,float)

    for i in range(natoms):
        line = filedesc.readline().split();
        names[i] = line[0]
        q[i] = line[1:4]
        orderp[i] = line[5] # The 6th column is the order parameter of each atom

    return [cell, names, q, orderp]

def FT_gibbs(phi, qxy, kgrid):
    # This is the un-normalized FT
    ng = len(kgrid)
    ak = np.zeros(ng,dtype=complex)
    
    n=0
    for k in kgrid:
        ak[n] = np.sum(phi[:]*np.exp(-1j*(qxy[:,0]*k[0]+qxy[:,1]*k[1])))
        n += 1
    return ak

def output_FT(outfile, kgrid, Ak, msg):
    # output kx, ky, the real and imaginary part of the Fourier coefficient
    outfile.write("#! Output the values of %s\n#! FIELDS k_x  k_y  A_real A_imag" % (msg))
    for i in range(len(Ak)):
        outfile.write("%10.6f  %10.6f  %12.10e  %12.10e\n" % (kgrid[i,0],kgrid[i,1],Ak[i].real,Ak[i].imag))
    return 0

def main(temperature):
    # collect all the input filenames
    filelist = []
    filelist.append(glob.glob('smapfcc*.xyz'))
    print(filelist)
    # Outputs
    ofile=open("allak.dat","w")
    avgfile=open("avgaksqr.dat","w")

    # set the reference values
    [phis, phil, vs] = [bulkproperties[temperature]["phis"], bulkproperties[temperature]["phil"], bulkproperties[temperature]["vs"]]
    #print(phis, phil, vs)

    # get total number of bins and initialize the grid
    binv = [20,20]; ndim=int(len(binv))
    if (ndim == 2):
        bins=np.prod(binv)
    else:
        print("Has to be 2-dimensional")
        sys.exit()
    kgrid = np.zeros((bins,2),float)

    # the main thing
    nframe = 0
    for ifile in filelist[0]:
        print("Reading file:", ifile)
        traj = open(ifile,"r")
        while True:
            try:
                [ cell, names, q, phi ] = read_frame(traj)
            except:
                break
            nframe += 1
            print("Frame No:", nframe)

            # normalize the order parameters so that <phis>=1 amd <phil>=0
            phi -= phil
            phi /= (phis-phil)

            if (nframe == 1):
                surfacearea = cell[0]*cell[1]
                # 1/(Area*phi_s*rho_s)
                normfactor = vs/surfacearea
                # initialize the k grid
                dkx = (2.*pi)/cell[0]
                dky = (2.*pi)/cell[1]
                n=0
                for i in range(binv[0]):
                    for j in range(binv[1]):
                        kgrid[n,0] = dkx * i
                        kgrid[n,1] = dky * j
                        n+=1
                ksqr = np.zeros(bins,float)
                for i in range(bins):
                    ksqr[i]=kgrid[i,0]**2.+kgrid[i,1]**2.
    
                # FT analysis
                ak = normfactor*FT_gibbs(phi, q[:,0:2], kgrid)
                allak = np.matrix(ak)
            else:
                ak = normfactor*FT_gibbs(phi, q[:,0:2], kgrid)
                allak = np.append(allak,np.matrix(ak),axis=0)
            output_FT(ofile, kgrid, ak, "A(kx,ky)")

    print("A total of data points ", len(allak))

    # compute also <A(k)^2>
    meanreal=np.mean(allak[:].real, axis=0, dtype=np.float64)
    meanimag=np.mean(allak[:].imag, axis=0, dtype=np.float64)
    varreal=np.var(allak[:].real, axis=0, dtype=np.float64, ddof=1)
    varimag=np.var(allak[:].imag, axis=0, dtype=np.float64, ddof=1)
    avgaksqr = np.zeros(bins,dtype=complex)

    for i in range(bins):
        # Bear in mind that for each capillary mode, each sin and cos component on average stores energy of kbt/2
        # so by analyzing the two components seperately one can accumulate more statistics
        avgaksqr[i] = (varreal[0,i]+meanreal[0,i]**2.)+1j*(varimag[0,i]+meanimag[0,i]**2.)
    output_FT(avgfile, kgrid, avgaksqr,"<A(kx,ky)^2>")

    ofile.close()
    avgfile.close()
    sys.exit()


if __name__ == '__main__':
    main(sys.argv[1])

# to use: python quick-direct-ft.py [temperature]
