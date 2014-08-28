from sklearn.pipeline import Pipeline
from mixtape.tica import tICA
from mixtape.featurizer import ContactFeaturizer
import pickle


pipeline = Pipeline([
    ('featurizer', ContactFeaturizer('closest-heavy')),
    ('tica', tICA(n_components=5)),
])

with open('contact-tica-n5.pkl', 'w') as f:
    pickle.dump(pipeline, f)
