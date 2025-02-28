import mdtraj as md
import numpy as np
import copy
import os
from sklearn.decomposition import PCA


atoms_for_PCA = "protein and element != 'H'"
output_csv = "pca_data.csv"


if os.path.exists(output_csv):
    os.remove(output_csv)
  
with open("pca_data.csv","a") as f:
    f.write("PC1,PC2,time,ligand\n") # header
    ref = None
    pca = PCA(n_components=2) # for 2D representation
    # trajectory frames will be superimposed on to the first frame of the first trajectory
    for i,(trj,top,name) in enumerate([
        ("md_1/r00/desmond_md_job_00_md2_m_naloxone_trj",
            "md_1/r00/desmond_md_job_00_md2_m_naloxone.pdb","m_Naloxone"),
        ("md_1/r00/desmond_md_job_00_md2_az617_trj",
            "md_1/r00/desmond_md_job_00_md2_az617.pdb","AZ617"),
        ]):
        
        print("processing", trj, "...")
        
        t = md.load_dtr(trj,top=top)
        # print(f"number of frames= {t.n_frames} residues= {t.n_residues} atoms= {t.n_atoms}")
        # print(atoms_for_PCA)
        
        aidx =  t.topology.select(atoms_for_PCA) # list of atom indices
        s = t.atom_slice(aidx)
        
        if i == 0: # set reference trajectory
            ref = copy.deepcopy(s)
        s.superpose(ref,frame=0)
        # print(f"number of frames= {s.n_frames} residues= {s.n_residues} atoms= {s.n_atoms}")

        if i == 0:
            # fit a PCA model
            # rd: reduced dimension in PC space
            rd = pca.fit_transform(s.xyz.reshape(s.n_frames, s.n_atoms*3))
        else:
            # transform coordinates according to the previously fitted PCA model
            rd = pca.transform(s.xyz.reshape(s.n_frames, s.n_atoms*3))
        
        ns = s.time.reshape(s.time.shape[0],1)
        ns = ns * 0.001 # (ns)
        lg = np.repeat(name, ns.shape[0]).reshape(ns.shape[0],1)
        
        dataframe  = np.concatenate([rd, ns, lg], axis=1)
        
        for pc1,pc2,ns,ligand in dataframe:
            f.write("{:.3f},{:.3f},{:.1f},{}\n".format(
                float(pc1),
                float(pc2),
                float(ns), 
                ligand ))