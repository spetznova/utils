from sklearn.pipeline import Pipeline
#from mixtape.featurizer import DihedralFeaturizer
from mixtape.featurizer import ContactFeaturizer
from mixtape.tica import tICA
from mixtape.utils import verbosedump, verboseload
from glob import glob
import mdtraj as md
from msmbuilder import Project
import numpy as np
from mixtape.ghmm import GaussianFusionHMM

prj  =  Project.load_from("ProjectInfo.yaml")
#feat = DihedralFeaturizer(['phi', 'psi'], sincos=True)
feat = ContactFeaturizer(contacts='all', scheme='closest-heavy')
tica = tICA(n_components=5, gamma=0,lag_time=20)
paths = glob('*.lh5')

output = {}

for path in np.arange(prj.n_trajs):
        featurized_path =  feat.partial_transform(prj.load_traj(path))
        try:
                tica.partial_fit(featurized_path)
        except:
                print "skipping",path

for path in np.arange(prj.n_trajs):
        featurized_path =  feat.partial_transform(prj.load_traj(path))
        output[path] = tica.partial_transform(featurized_path)

# save output
verbosedump(output, 'my-tics.pkl')
verbosedump(tica, 'tica-obj.pkl')

#Ignore below; the model fitting is very quick and can be done in iPython

#get the data
#X =  [ output[i] for i in output.iterkeys()]

#model = GaussianFusionHMM(n_states=5, n_features=5)
#model.fit(X)

#save model
#verbosedump(model,"ghmm_s5_n5.pkl")


