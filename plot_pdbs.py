import mdtraj as mdt
import os
from glob import glob
import sys

import matplotlib
matplotlib.use('Agg')
#import matplotlib.pyplot as plt
from pylab import *

if len(sys.argv) > 1:
	os.chdir(sys.argv[1])

trj = mdt.load(glob('*.pdb'))
ref = mdt.load('/home/willf/analysis/ago2rh1/2rh1_common_residues_crystal_apo.pdb')

pwd = os.getcwd()
tit = pwd[pwd.find('data'):pwd.find('/PDB')]


def resid_to_index(resid,atom_wanted):
    '''
    function to return atom indices associated with a atom(CA,N etc)
    '''
    for a in trj.topology.residue(resid).atoms:
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

npxxy_rms = rmsd(trj,ref,npxxy_list)

arg_131_ca = resid_to_index(100,'CA')
leu_272_ca = resid_to_index(201,'CA')
r131 = trj.xyz[:,arg_131_ca,:]
l272 = trj.xyz[:,leu_272_ca,:]
dis = np.zeros(len(trj))
for i in np.arange(len(trj)):
    dis[i] = np.linalg.norm(r131[i,:]-l272[i,:])

color_list = np.zeros(len(trj))
for i in np.arange(len(trj)/100):
    color_list[i*100:(i+1)*100] = i 

#fig = figure()
figure(figsize=(10,8))
scatter(dis*10,npxxy_rms*10,c=color_list,s=100)
#fig.plot(dis*10,npxxy_rms*10,'.',c=color_list,s=100)
xlabel("Helix 6- Helix3")
ylabel("NPxxY to inactive")
title(tit)
#savefig(pwd,dpi=300)
savefig('npxxy_vs_h3h6.png', bbox_inches='tight', dpi=300)

