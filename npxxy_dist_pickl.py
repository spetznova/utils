import mdtraj
import mdtraj.io
import numpy as np
import pickl

ref = mdtraj.load('/home/willf/analysis/ago2rh1/2rh1_common_residues_crystal.pdb')

def resid_to_index(resid,atom_wanted):
    '''
    function to return atom indices associated with a atom(CA,N etc)
    '''

    for a in ref.topology.residue(resid).atoms:
        if a.name == atom_wanted:
            return a.index
            break

    return Exception

def rmsd(traj, ref, idx):
    traj = traj.superpose(ref,atom_indices=idx)
    return np.sqrt(np.sum(np.square(traj.xyz[:,idx,:] - ref.xyz[:,idx,:]),axis=(1,2))/len(idx))

npxxy_list=[]
for i in range(251,257):
    npxxy_list += resid_to_index(i,'CA'),resid_to_index(i,'N'),resid_to_index(i,'C'),resid_to_index(i,'O')

arg_131_ca = resid_to_index(100,'CA')
leu_272_ca = resid_to_index(201,'CA')

assignments = mdtraj.io.loadh('/home/willf/analysis/ago2rh1/tica/dist20E/data_s1n5k100_d10_dist20EAssignments.h5')['arr_0']

npxxy_rmsd = -1*np.ones(shape(assignments))
dist = -1*np.ones(shape(assignments))
for n in range(39999):
    file = "/home/willf/Trajectories/trj" + str(n) + ".lh5"
    trj = mdtraj.load(file)
    
    npxxy_rmsd[n,:len(trj)] = rmsd(trj,ref,npxxy_list)
    
    r131 = trj.xyz[:,arg_131_ca,:]
    l272 = trj.xyz[:,leu_272_ca,:]
    for i in range(len(trj)):
        dist[n,i] = np.linalg.norm(r131[i,:]-l272[i,:])

pickle.dump(npxxy_rmsd, open( "npxxy_trajs.pickl", "wb"))
pickle.dump(dist, open( "dist_trajs.pickl", "wb"))
