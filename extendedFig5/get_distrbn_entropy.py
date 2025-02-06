import numpy as np
import freud
import gsd
import gsd.hoomd
import sys
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import entropy

font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 22}

matplotlib.rc('font', **font)
#plt.rcParams.update({'font.size': 22})

trajin1 = sys.argv[1]
trajin2 = sys.argv[2]
trajin3 = sys.argv[3]
trajectory1 = gsd.hoomd.open(trajin1, mode = "rb")
trajectory2 = gsd.hoomd.open(trajin2, mode = "rb")
trajectory3 = gsd.hoomd.open(trajin3, mode = "rb")

#snap = trajectory[-1]
#box = snap.configuration.box
#npart = snap.particles.N
"""
rdf = freud.density.RDF(r_max = 1000.0,bins=200,normalize=True)
ref_points = np.array([positions[0]])
#print(ref_points)
rdf.compute(system=(box,positions))#,query_points=ref_points)
#print(rdf.bin_counts)
#print(rdf.bin_centers)
a = np.concatenate((rdf.bin_centers.reshape(-1,1),rdf.bin_counts.reshape(-1,1)),axis=1)
#print(a)
plt.plot(rdf.bin_centers,rdf.bin_counts)
plt.show()
"""

def get_distributions(istart, iend,traj):
    r_state = []
    for snap in range(istart, iend):
        positions = traj[snap].particles.position
        npart = len(positions)
        for i in range(npart):
            pos = positions[i]
            pos_brd = np.broadcast_to(pos,(npart-1,3))
            position_ = np.delete(positions, i, axis=0)
            dr = pos_brd - position_
            r = np.sqrt(np.sum(dr**2,axis=1))
            r_nbr = r[r<245]
            if len(r_nbr) > 6:
                r_state.extend(r[r<1000])
    r_state = np.array(r_state)
    hist, bin_edges = np.histogram(r_state,range=(0,1000),bins=100,density=True)
    bins = (bin_edges[1:] + bin_edges[:-1])/2
    return hist, bins 




#def check_distributions(i,positions, npart):
#    r_state = []
#    pos = positions[i]
#    pos_brd = np.broadcast_to(pos,(npart-1,3))
#    position_ = np.delete(positions, i, axis=0)
#    dr = pos_brd - position_
#    r = np.sqrt(np.sum(dr**2,axis=1))
#    r_nbr = r[r<245]
#    if len(r_nbr) > 6:
#        r_state.extend(r[r<1000])
#        r_state = np.array(r_state)
#        hist, bin_edges = np.histogram(r_state,range=(0,1000),bins=100,density=True)
#        bins = (bin_edges[1:] + bin_edges[:-1])/2
#        return hist, bins
#"""

l_start = -10
l_end = -1
c_start = -10
c_end = -1

hist_liq, bins = get_distributions(l_start, l_end, trajectory1)
hist_crs, bins = get_distributions(c_start, c_end, trajectory2)
np.savetxt("reference_liquid_1000nm.txt", np.concatenate((bins.reshape(-1,1),hist_liq.reshape(-1,1)),axis=1))
np.savetxt("reference_crystal_1000nm.txt", np.concatenate((bins.reshape(-1,1),hist_crs.reshape(-1,1)),axis=1))

fig = plt.figure(figsize =(8,6))
plt.plot(bins,hist_liq,label="Liquid",linewidth=5)
plt.plot(bins,hist_crs,label="Crystal",linewidth=5)
plt.xlabel("Distance between particles (nm)", fontweight = "bold")
plt.ylabel("Distance distributions", fontweight = "bold")
plt.legend()
plt.savefig("reference_distributions_new.eps",format='eps', bbox_inches='tight',dpi=2000)
#plt.show()

#fout = open("snap_4820_ver2.xyz", mode = "w")

hist_liq[hist_liq==0]=0.000001
hist_crs[hist_crs==0]=0.000001
snaps = [-1]
for snap in snaps:
    fout = open("entropy_distrbn{0}_liquid.out".format(int(snap)),mode = "w")
    positions = trajectory3[snap].particles.position
    fig = plt.figure(figsize=(7,6))
    entropy_liq = []
    entropy_crys = []
    npart = len(positions)
    j = 0
    for i in range(0,npart,100):
        j = j + 1
        r_state = []
        pos = positions[i]
        pos_brd = np.broadcast_to(pos,(npart-1,3))
        position_ = np.delete(positions, i, axis=0)
        dr = pos_brd - position_
        r = np.sqrt(np.sum(dr**2,axis=1))
        r_nbr = r[r<245]
        if len(r_nbr) > 6:
            r_state.extend(r[r<1000])
            r_state = np.array(r_state)
            hist_i, bin_edges = np.histogram(r_state,range=(0,1000),bins=100,density=True)
            bins = (bin_edges[1:] + bin_edges[:-1])/2
            x = hist_i / hist_liq
            y = hist_i / hist_crs
            x[x==0] = 0.0000001
            y[y==0] = 0.0000001

            ent_l = np.sum(hist_i*np.log(x))
            ent_c = np.sum(hist_i*np.log(y))
            #if ent_l < ent_c:
            #    if trajectory[snap].particles.typeid[i]==0:
            #        print("P", pos[0], pos[1], pos[2], file = fout)
            #    if trajectory[snap].particles.typeid[i]==1:
            #        print("N", pos[0], pos[1], pos[2], file = fout)
            #if ent_l > ent_c:
            #    if trajectory[snap].particles.typeid[i]==0:
            #        print("Cs", pos[0], pos[1], pos[2], file = fout)
            #    if trajectory[snap].particles.typeid[i]==1:
            #        print("Cl", pos[0], pos[1], pos[2], file = fout)


            entropy_liq.append(np.array([j,ent_l]))
            entropy_crys.append(np.array([j,ent_c]))
            #entropy_liq.append(np.array([i, entropy(hist_i, hist_liq)]))
            #entropy_crys.append(np.array([i, entropy(hist_i, hist_crs)]))
            #print("entropy between part", i,"and liq state", entropy(hist_i, hist_liq))
            #print("entropy between part", i,"and crs state", entropy(hist_i, hist_crs))
            #plt.plot(bins,hist_i)
        #ent_liq = entropy(hist_liq, hist_i)
        #ent_crs = entropy(hist_crs, hist_i)
        #print(i, ent_liq, ent_crs, file = fout)


entropy_liq = np.array(entropy_liq)
entropy_crys = np.array(entropy_crys)
print(entropy_liq)
print(entropy_crys)
fig = plt.figure(figsize=(8,6))
plt.scatter(entropy_liq[:,0], entropy_liq[:,1],label = "Entropy_liquid")
plt.scatter(entropy_crys[:,0], entropy_crys[:,1],label="Entropy_crystal")
plt.xlabel("Particle index", fontweight="bold")
plt.ylabel("Entropy between particle's and \nreference distribution (distance)", fontweight="bold")
plt.ylim(0,0.04)
plt.legend(fontsize=12)
plt.savefig("entrpy_perpart_crystal_new.eps", bbox_inches='tight',format='eps', dpi=2000)
#plt.show()

#plt.plot(bins,hist_liq,label="Liquid")
#plt.plot(bins,hist_crs,label="Crystalline")
#plt.legend()
#plt.show()


