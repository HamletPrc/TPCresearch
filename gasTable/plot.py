import matplotlib.pyplot as plt
import numpy as numpy

filename = 'temp.dat'
X,Y = [],[]
with open(filename, 'r') as f:#1
    lines = f.readlines()#2
    for line in lines:#3
        value = [float(s) for s in line.split()]#4
        X.append(value[0])#5
        Y.append(value[1])

print(X)
print(Y)

plt.loglog(X, Y, marker='o')
plt.xlabel('Field Strength/(V/cm)')
plt.ylabel('Electron Drift Velocity/(cm/us)')
plt.title('Electron Drift Velocity for T2K')
plt.grid()
plt.legend()
plt.show()

