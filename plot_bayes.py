import matplotlib
matplotlib.use('Agg')
#import matplotlib.pyplot as plt
from pylab import *


bace = np.loadtxt('bayesFactors.dat')
figure(figsize=(10,8))
plot(bace[:,0],log(bace[:,1]),'.-')
xlabel("# Macrostates")
ylabel("log(Bayes factor)")
savefig('bayesfactor.png', bbox_inches='tight')
