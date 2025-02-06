import numpy as np
import os
import sys
import gsd
import gsd.hoomd
import math
import scipy
import scipy.special
#from numba import jit
import matplotlib.pyplot as plt
import matplotlib

font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 22}

matplotlib.rc('font', **font)

traj_file = sys.argv[1]
traj = gsd.hoomd.open(traj_file, "rb")

#print(len(traj))
#numsnap = len(traj)

traj_file2 = sys.argv[2]
traj2 = gsd.hoomd.open(traj_file2, "rb")

traj_file3 = sys.argv[3]
traj3 = gsd.hoomd.open(traj_file3, "rb")


#frame = 0
#snap = traj[frame]
#box = snap.configuration.box
#pos = snap.particles.position
#print("number of prticles",len(pos))
#boxL = box[0]

rcut = 300.0
l = 6
pi = 3.14
ncut = 11
m = np.linspace(-l,l,num=2*l+1)

#print(m)


#istart = 0
#iend = 3025

def get_avg_q(positions,l,m,rcut, ncut,pos_avgq):
    npart = len(positions)
    avgqt = 0
    ncount = 0
    sph_real = []
    sph_imag = []
    part_index = []
    nbor_index = []

    for j in range(0,npart):
        ylmall = []
        #nborind_init = []
        posj = positions[j]
        posjbcst = np.broadcast_to(posj,(npart,3))
        dr = posjbcst - positions
        #dr = np.delete(dr,j)
        r = np.sqrt(np.sum(dr**2,axis=1))
        indn = np.where(r < rcut)[0]
        nborind_init = np.delete(indn,np.where(indn==j))
        #print(j,nborind_init)
        numnbor = len(nborind_init)
        dr_nbor = dr[nborind_init]
        r_nbor = r[nborind_init]
        x = dr_nbor[:,2] / r_nbor
        x[x>1.0] = 1.0
        x[x<-1.0] = -1.0
        theta = np.arccos(x)
        phi = np.arctan2(dr_nbor[:,1], dr_nbor[:,0])
        for k in range(numnbor):
            phik = phi[k]
            thetak = theta[k]
            ylm = scipy.special.sph_harm(m, l, phik, thetak)
            ylm.real[np.isnan(ylm.real)==True]=0
            ylm.imag[np.isnan(ylm.imag)==True]=0
            ylmall.append(ylm)
        ylmall = np.array(ylmall)
    
        if len(ylmall) != 0:
            ncount = ncount + 1
            num_neighbor, num_m = ylmall.shape
            part_index.append(j)
            nbor_index.append(nborind_init)
            qlm_real = np.zeros(num_m)
            qlm_imag = np.zeros(num_m)
            qlm = np.zeros(num_m)
            for n in range(num_m):
                qlm_real[n] = np.sum(ylmall[:,n].real)
                qlm_imag[n] = np.sum(ylmall[:,n].imag)
                qlm_real[n] = qlm_real[n] / num_neighbor
                qlm_imag[n] = qlm_imag[n] / num_neighbor

            sph_real.append(qlm_real)
            sph_imag.append(qlm_imag)
            #print(sph_real)
            #print(sph_imag)


    nbor_index = np.array(nbor_index,dtype=object)
    sph_real = np.array(sph_real)
    sph_imag = np.array(sph_imag)
    #print("shape of nbor_index", nbor_index.shape)
    #print("shape of sph_real", sph_real.shape)
    #print("shape of sph_imag", sph_imag.shape)
    #print(nbor_index)
    for i in range(ncount):
        num_nbor = len(nbor_index[i])
        indices = nbor_index[i]
        if num_nbor > ncut:
            ic = part_index[i]
            real_ = np.zeros(num_m)
            imag_ = np.zeros(num_m)
            indices2 = []
            for j in range(len(indices)):
                indices2.append(np.where(part_index==indices[j])[0])
            indices2 = np.array(indices2)
            for j in range(num_m):
                real_[j] = sph_real[i,j]+np.sum(sph_real[indices2,j])
                imag_[j] = sph_imag[i,j]+np.sum(sph_imag[indices2,j])
                real_[j] = real_[j] / (num_nbor + 1)
                imag_[j] = imag_[j] / (num_nbor + 1)
            avgq = np.sum(real_**2+imag_**2)
            avgq = math.sqrt((4 * 3.14 * avgq) / num_m) 
            pos_avgq.append(avgq)

    return(pos_avgq)

#meanqall = []
#fout = open("meanqall_0_100.out", "w")
pos_avgq = []
for i in range(-5,-1):
    print(i)
    snap = traj[i]
    positions = snap.particles.position
    get_avg_q(positions,l,m,rcut,ncut,pos_avgq)
    #pos_avgq = np.array(pos_avgq)
    #np.savetxt("avgq_snap{0}.txt".format(int(i)),pos_avgq)
    print(i, "done")

pos_avgq = np.array(pos_avgq)
hist, bin_edges = np.histogram(pos_avgq, range=(0,0.8), bins=50)
hist = hist / len(pos_avgq)
bins = (bin_edges[:-1] + bin_edges[1:]) / 2

pos_avgq2 = []
for i in range(-5,-1):
    print(i)
    snap = traj2[i]
    positions = snap.particles.position
    get_avg_q(positions,l,m,rcut,ncut,pos_avgq2)
#    #pos_avgq = np.array(pos_avgq)
#    #np.savetxt("avgq_snap{0}.txt".format(int(i)),pos_avgq)
#    print(i, "done")

pos_avgq2 = np.array(pos_avgq2)
hist2, bin_edges2 = np.histogram(pos_avgq2, range=(0,0.8), bins=50)
bins2 = (bin_edges2[:-1] + bin_edges2[1:]) / 2
hist2 = hist2 / len(pos_avgq2)

pos_avgq3 = []
for i in range(1570,1571):
    print(i)
    snap = traj3[i]
    positions = snap.particles.position
    get_avg_q(positions,l,m,rcut,ncut,pos_avgq3)
#    #pos_avgq = np.array(pos_avgq)
#    #np.savetxt("avgq_snap{0}.txt".format(int(i)),pos_avgq)
#    print(i, "done")

pos_avgq3 = np.array(pos_avgq3)
hist3, bin_edges3 = np.histogram(pos_avgq3, range=(0,0.8), bins=50)
bins3 = (bin_edges3[:-1] + bin_edges3[1:]) / 2
hist3 = hist3 / len(pos_avgq3)

fig = plt.figure(figsize=(8,6))
plt.plot(bins, hist, label = "CsCl-like", linewidth=5)
plt.plot(bins2, hist2, label = "Liquid", linewidth=5)
plt.plot(bins3, hist3, linestyle='--',linewidth=5)
plt.xlabel(r'$\bar{q_{6}}$', fontweight="bold")
plt.ylabel("Distribution", fontweight="bold")
plt.legend(fontsize=12)
plt.savefig("q6_bar.eps", bbox_inches='tight', format='eps', dpi=2000)








