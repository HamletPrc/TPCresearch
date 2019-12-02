import matplotlib.pyplot as plt
import numpy as numpy

filename1 = 'Ar-CO2aa'
filename2 = 'Ar-CO2ab'
filename3 = 'Ar-CO2ac'
filename4 = 'Ar-CO2ad'
filename5 = 'Ar-CO2ae'
filename6 = 'Ar-CO2af'
filename7 = 'Ar-CO2ag'
filename8 = 'Ar-CO2ah'
filename9 = 'Ar-CO2ai'
filename10 = 'Ar-CO2aj'
filename11 = 'Ar-CO2ak'
X, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9, Y10, Y11 = [],[],[],[],[],[],[],[],[],[],[],[]
with open(filename1, 'r') as f:#1
    lines = f.readlines()#2
    for line in lines:#3
        value = [float(s) for s in line.split()]#4
        X.append(value[0])#5
        Y1.append(value[1])
with open(filename2, 'r') as f:#1
    lines = f.readlines()#2
    for line in lines:#3
        value = [float(s) for s in line.split()]#4
        Y2.append(value[1])
with open(filename3, 'r') as f:#1
    lines = f.readlines()#2
    for line in lines:#3
        value = [float(s) for s in line.split()]#4
        Y3.append(value[1])
with open(filename4, 'r') as f:#1
    lines = f.readlines()#2
    for line in lines:#3
        value = [float(s) for s in line.split()]#4
        Y4.append(value[1])
with open(filename5, 'r') as f:#1
    lines = f.readlines()#2
    for line in lines:#3
        value = [float(s) for s in line.split()]#4
        Y5.append(value[1])
with open(filename6, 'r') as f:#1
    lines = f.readlines()#2
    for line in lines:#3
        value = [float(s) for s in line.split()]#4
        Y6.append(value[1])
with open(filename7, 'r') as f:#1
    lines = f.readlines()#2
    for line in lines:#3
        value = [float(s) for s in line.split()]#4
        Y7.append(value[1])
with open(filename8, 'r') as f:#1
    lines = f.readlines()#2
    for line in lines:#3
        value = [float(s) for s in line.split()]#4
        Y8.append(value[1])
with open(filename9, 'r') as f:#1
    lines = f.readlines()#2
    for line in lines:#3
        value = [float(s) for s in line.split()]#4
        Y9.append(value[1])
with open(filename10, 'r') as f:#1
    lines = f.readlines()#2
    for line in lines:#3
        value = [float(s) for s in line.split()]#4
        Y10.append(value[1])
with open(filename11, 'r') as f:#1
    lines = f.readlines()#2
    for line in lines:#3
        value = [float(s) for s in line.split()]#4
        Y11.append(value[1])


print(X)
print(Y1)
print(Y2)
print(Y3)
print(Y4)
print(Y5)
print(Y6)
print(Y7)
print(Y8)
print(Y9)
print(Y10)
print(Y11)


plt.plot(X, Y1, marker='o', label='89/11')
plt.plot(X, Y2, marker='o', label='90/10')
plt.plot(X, Y3, marker='o', label='91/9')
plt.plot(X, Y4, marker='o', label='92/8')
plt.plot(X, Y5, marker='o', label='93/7')
plt.plot(X, Y6, marker='o', label='94/6')
plt.plot(X, Y7, marker='o', label='95/5')
plt.plot(X, Y8, marker='o', label='96/4')
plt.plot(X, Y9, marker='o', label='97/3')
plt.plot(X, Y10, marker='o', label='98/2')
plt.plot(X, Y11, marker='o', label='99/1')
plt.gca().set_xlim([100,1000])

plt.xlabel('Field Strength/(V/cm)')
plt.ylabel('Electron Drift Velocity/(cm/us)')
plt.title('Electron Drift Velocity for Ar-CO2')
plt.grid()
plt.legend(loc=3)
plt.show()
