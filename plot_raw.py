import mdtraj as mdt
import os
import matplotlib
matplotlib.use('Agg')
from pylab import *
import random
from glob import glob
import sys



if len(sys.argv) > 1:
	os.chdir(sys.argv[1])

#p = xrange(len(glob('*.lh5'))
#s = random.sample(p,100)
s = glob('*.lh5')
random.shuffle(s)
s = s[:100]



count = 0
trj = mdt.load("trj" + str(s[0]) + ".lh5")
for n in s[1:]:
    trj = trj + mdt.load("trj" + str(n) + ".lh5")
    count += 1
    if count%100 == 0:
        print count


def resid_to_index(resid,atom_wanted):
from pylab import *    '''
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

ago = 0

if not ago:

	print 'APO TYPE DETECTED'

	ref = mdt.load("/home/willf/xtals_diwakar/2rh1_common_residues_crystal_apo.pdb")

	npxxy_list=[]
	for i in range(250,256):
    		npxxy_list += resid_to_index(i,'CA'),resid_to_index(i,'N'),resid_to_index(i,'C'),resid_to_index(i,'O')

	arg_131_ca = resid_to_index(100,'CA')
	leu_272_ca = resid_to_index(200,'CA')



if ago:
	
	print 'AGONIST-BOUND DETECTED'

	ref = mdt.load("/home/willf/xtals_diwakar/2rh1_common_residues_crystal.pdb")

	npxxy_list=[]
	for i in range(251,257):
    		npxxy_list += resid_to_index(i,'CA'),resid_to_index(i,'N'),resid_to_index(i,'C'),resid_to_index(i,'O')

	arg_131_ca = resid_to_index(100,'CA')
	leu_272_ca = resid_to_index(201,'CA')

npxxy_rms = rmsd(trj,ref,npxxy_list)
r131 = trj.xyz[:,arg_131_ca,:]
l272 = trj.xyz[:,leu_272_ca,:]
dis = np.zeros(len(trj))
for i in np.arange(len(trj)):
	dis[i] = np.linalg.norm(r131[i,:]-l272[i,:])

figure(figsize=(10,8))
scatter(dis*10,npxxy_rms*10,s=1)#s=np.arange(0,5000)/100)
xlabel("Helix 6- Helix3")
ylabel("NPxxY to Inactive")
#title(tit)
savefig('raw_npxxy_vs_h3h6.png', bbox_inches='tight', dpi=300)
