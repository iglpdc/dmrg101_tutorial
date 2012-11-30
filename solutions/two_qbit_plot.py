import matplotlib.pyplot as plt

fin = open('two_qbit_entropies.dat','r')
x = []
y = []
for lines in fin:
    line = lines.split()
    x.append(float(line[0]))
    y.append(float(line[1]))
plt.plot(x,y)
plt.show()
