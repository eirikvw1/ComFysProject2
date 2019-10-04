import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from math import *

# Reading from file
f=open("Egenverdier.txt", "r")
lines =f.readlines()[1:]
n = int(lines[0].split()[0])
iter = int(lines[0].split()[1])


f.close()
eigval = np.zeros(n)
eigval_rot =np.zeros(n)
analytisk=np.zeros(n)
x = np.zeros(n)


# Reading out values
i = 0
an=0
for line in lines[1:]:
    words = line.split()
    eigval[i] = float(words[0])
    eigval_rot[i] = float(words[1])
    analytisk[i]=3+an
    an = an+4
    i = i + 1

# Plotting
# plt.hold(True)
plt.plot(eigval,'.r')
plt.plot(eigval_rot,'.b')
plt.plot(analytisk,'.g')
plt.xlabel('points')
plt.ylabel('Eigenvalues')
plt.legend(['eigval', 'eigval_rot','analytisk'])
title = 'Numerical solutions to eigenvalus, n ='+str(n)+', # of iterations: '+str(iter)
plt.title(title)
plt.savefig('project2_100__100000_5png')
