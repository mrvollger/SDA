#!/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
import argparse 


parser = argparse.ArgumentParser()
parser.add_argument("--nucfreq", help="assembly.consensus.nucfreq", default="assembly.consensus.nucfreq")
parser.add_argument("--automin", help="autoMinCoverage", default="autoMinCoverage")
parser.add_argument("--automax", help="autoMaxCoverage", default="autoMaxCoverage")
args = parser.parse_args()

nucfreq=args.nucfreq
autoMin=args.automin
autoMax=args.automax


colnames = ["contig", "pos", "A", "C", "G", "T", "deletion", "insertion"]

f = open(nucfreq)
first  = []
second = []
truepos= []
for line in f:
    line = line.split()
    truepos.append(int(line[1]))
    bases = []
    for basepair in line[2:6]:
        bases.append(int(basepair))
    bases = sorted(bases, reverse=True)
    first.append(bases[0])
    second.append(bases[1])
pos = np.array( range(0,len(second)) ) 
second = np.array(second)
first = np.array(first)
truepos = np.array(truepos)




plt.plot(truepos, first, 'bo')
plt.plot(truepos, second, 'ro')
plt.savefig('threshold.png')

exit(0)
# after this stuff is gmm stuff that did not work

#second[second > 300 ] = 15
gmm = GaussianMixture(n_components=2, covariance_type="diag", tol=0.001, weights_init= [.5,.5])
gmm = gmm.fit(X=np.expand_dims(second, 1))

varrs = gmm.covariances_
stds = np.sqrt(varrs)
mean = max(gmm.means_)[0]
sd = max(stds)[0]

print("mean: {} std: {} weights {}".format(mean, sd, gmm.weights_))

# but if there is a duplicaiton or something the number of errors doubles because errors from
# both regiosn map there. So I have to make the threshold much higher
minCoverage = (sd  + mean)  
numAbove = sum(second > minCoverage)
print("Threshhold: {} \t Number of points above threshhold: {}".format(minCoverage, numAbove))


# reset stuff 
psvs = second[second > minCoverage]
psvPos = pos[second > minCoverage]
noise = second[second <= minCoverage]
noisePos = pos[second <= minCoverage]


psvGmm = GaussianMixture(n_components=1, covariance_type="diag", tol=0.001)
psvGmm = psvGmm.fit(X=np.expand_dims(psvs, 1))
varrs = psvGmm.covariances_
stds = np.sqrt(varrs)
mean = max(psvGmm.means_)[0]
sd = max(stds)[0]
# if I look at things 5 or more std distributions away I should only have a .0001  error rate 
maxCoverage = mean + sd * 5 
print("mean: {} std: {}".format(mean, sd))


plt.plot(noisePos, noise, 'ro')
plt.plot(psvPos, psvs, 'go')
plt.plot(pos, first, 'bo')
plt.axhline(y=minCoverage, color='g', linestyle='--')
plt.axhline(y=maxCoverage, color='g', linestyle='--')

#plt.show()
plt.savefig('threshold.png')


minCoverage = int(minCoverage)
maxCoverage = int(maxCoverage)
maxFile = open(autoMax, "w+")
minFile = open(autoMin, "w+")
maxFile.write(str(maxCoverage) + "\n" )
minFile.write(str(minCoverage) + "\n" )

print("min/maxCoverage: {}/{}".format((minCoverage), (maxCoverage)))
print( "PLEASE LOOK AT threshold.png TO CONFIRM AUTO THRESHOLDS" )


