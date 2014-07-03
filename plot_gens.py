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

pwd = os.getcwd()
tit = pwd[pwd.find('data'):]

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

ago = 1
if pwd.find('ago') < 0:
	ago = 0

if not ago:

	print 'APO TYPE DETECTED'

	ref = mdt.load("/home/willf/xtals_diwakar/2rh1_common_residues_crystal_apo.pdb")

	npxxy_list=[]
	for i in range(250,256):
    		npxxy_list += resid_to_index(i,'CA'),resid_to_index(i,'N'),resid_to_index(i,'C'),resid_to_index(i,'O')

	#gens = mdt.io.loadh('Gens.h5')
	#gens = mdt.load('Gens.h5')
	#mytraj = mdt.Trajectory(gens.xyz,ref.topology)
	gens = mdt.load('Gens.h5')

	arg_131_ca = resid_to_index(100,'CA')
	leu_272_ca = resid_to_index(200,'CA')

	npxxy_rms = rmsd(gens,ref,npxxy_list)
	r131 = gens.xyz[:,arg_131_ca,:]
	l272 = gens.xyz[:,leu_272_ca,:]
	dis = np.zeros(len(npxxy_rms))
	for i in np.arange(len(npxxy_rms)):
    		dis[i] = np.linalg.norm(r131[i,:]-l272[i,:])

	figure(figsize=(10,8))
	scatter(dis*10,npxxy_rms*10,s=1)#s=np.arange(0,5000)/100)
	xlabel("Helix 6- Helix3")
	ylabel("NPxxY to Inactive")
	title(tit)
	savefig('gens_npxxy_vs_h3h6.png', bbox_inches='tight', dpi=300)

if ago:
	
	print 'AGONIST-BOUND DETECTED'

	ref = mdt.load("/home/willf/xtals_diwakar/2rh1_common_residues_crystal.pdb")

	npxxy_list=[]
	for i in range(251,257):
    		npxxy_list += resid_to_index(i,'CA'),resid_to_index(i,'N'),resid_to_index(i,'C'),resid_to_index(i,'O')

	#gens = mdt.io.loadh('Gens.h5')
	#gens = mdt.load('Gens.h5')
	#mytraj = mdt.Trajectory(gens.xyz,ref.topology)
	gens = mdt.load('Gens.h5')

	arg_131_ca = resid_to_index(100,'CA')
	leu_272_ca = resid_to_index(201,'CA')

	npxxy_rms = rmsd(gens,ref,npxxy_list)
	r131 = gens.xyz[:,arg_131_ca,:]
	l272 = gens.xyz[:,leu_272_ca,:]
	dis = np.zeros(len(npxxy_rms))
	for i in np.arange(len(npxxy_rms)):
    		dis[i] = np.linalg.norm(r131[i,:]-l272[i,:])

	figure(figsize=(10,8))
	scatter(dis*10,npxxy_rms*10,s=1)#s=np.arange(0,5000)/100)
	xlabel("Helix 6- Helix3")
	ylabel("NPxxY to Inactive")
	title(tit)
	savefig('gens_npxxy_vs_h3h6.png', bbox_inches='tight', dpi=300)



