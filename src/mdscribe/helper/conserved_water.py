import glob
import sys
import os
import numpy as np
from scipy.spatial import distance

grid_space = 1.0
threshold = 80

"""
Preprocessing:

Align Desmond MD trajectory to a reference
$ bash trj_align.sh

Extract water molecules within 2A from protein and write to pdb files
$ vmdcli ~/bucket/namd/vmd_water.tcl
"""

""" preparing water coordinates for density analysis """
wat = []
dirname = sys.argv[1]
for filename in glob.glob(os.path.join(dirname,"water_*.pdb")):
    with open(filename,"r") as pdb:
        for line in pdb:
            if line.startswith('ATOM'):
                line = line.strip()
                # PDB format version 2.3
                #serial = line[6:12]
                name = line[12:16].strip()
                #altLoc = line[16:17]
                #if altLoc=="": altLoc=" "
                resName = line[17:21].strip()
                #chainId = line[21:22]
                #if chainId=="": chainId=" "
                # resSeq = line[22:26].strip()
                #iCode = line[26:27]
                #if iCode=="": iCode=" "
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                #occupancy = line[54:60]
                #tempFactor = float(line[60:66])
                #segId = line[72:76]
                #element = line[76:78]
                #charge = line[78:80]
                if resName == "T3P" and name == "O":
                    wat.append([x,y,z])

""" generating 3D histogram (binning) """
centers = []
frequency = []
xyz = np.array(wat)
bins = np.ceil((np.amax(xyz, axis=0)-np.amin(xyz, axis=0))/grid_space)
H, edges = np.histogramdd(xyz, bins=bins)
# print(H.shape)
# print(edges[0].shape, edges[1].shape, edges[2].shape)
# print(edges[0])
for idx in np.argsort(H, axis=None)[-1::-1]:
    ix,iy,iz = np.unravel_index(idx, shape=H.shape)
    if H[(ix,iy,iz)] > threshold:
        frequency.append((
            edges[0][ix], edges[0][ix+1], 
            edges[1][iy], edges[1][iy+1], 
            edges[2][iz], edges[2][iz+1],
            H[(ix,iy,iz)]
        ))
        centers.append((
            (edges[0][ix]+edges[0][ix+1])/2.0,
            (edges[1][iy]+edges[1][iy+1])/2.0,
            (edges[2][iz]+edges[2][iz+1])/2.0,
        ))
        print("{:2} X {:7.2f} {:7.2f} Y {:7.2f} {:7.2f} Z {:7.2f} {:7.2f} N {:6}".format(
            len(frequency),
            edges[0][ix], edges[0][ix+1], 
            edges[1][iy], edges[1][iy+1], 
            edges[2][iz], edges[2][iz+1],
            H[(ix,iy,iz)]
        ))

""" distance analysis """
# points = np.array(center)
# D = distance.cdist(points, points, 'euclidean')
# for idx in np.argsort(D, axis=None):
#     ir, ic = np.unravel_index(idx, shape=D.shape)
#     if ir != ic and D[(ir,ic)] < 2:
#         print("{:6} {:6} {:5.2f} {} {}".format(ir, ic, D[(ir,ic)],frequency[ir][-1], frequency[ic][-1]))
        # print("X {:7.2f} {:7.2f} Y {:7.2f} {:7.2f} Z {:7.2f} {:7.2f} N {:6}".format(frequency[ir]))
        # print("X {:7.2f} {:7.2f} Y {:7.2f} {:7.2f} Z {:7.2f} {:7.2f} N {:6}".format(frequency[ic]))
        # print()

with open("critical_water.pdb","w") as pdbout:
    altLoc = " "
    chainId = "W"
    iCode = " "
    charge = " "
    segId = "WAT"
    resName = "HOH"
    name = "O"
    element = "O"
    for idx, (x,y,z) in enumerate(centers):
        serial = idx
        resSeq = 9000 + idx
        occupancy = frequency[idx][-1]
        tempFactor = frequency[idx][-1]
        pdbout.write("%-6s%5d %-4s%c%3s " % ("ATOM",serial,name,altLoc,resName))
        pdbout.write("%c%4d%c   " % (chainId,resSeq,iCode))
        pdbout.write("%8.3f%8.3f%8.3f" % (x,y,z))
        pdbout.write("%6.2f%6.2f      " % (occupancy,tempFactor))
        pdbout.write("%4s%2s%2s" % (segId,element,charge))
        pdbout.write("\n")