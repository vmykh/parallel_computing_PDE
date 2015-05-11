import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def drange(start, n, step):
	result = float(start)
	for i in xrange(0, int(n)):
		yield result
		result += float(step)

def readDataFromFile():
	return [line.strip().split() for line in open('data.dat')]

def getDataFromFile():
	lines = readDataFromFile()
	xs = list(drange(lines[0][0], lines[0][1], lines[0][2]))
	ts = list(drange(lines[1][0], lines[1][1], lines[1][2]))

	Ui = []

	for ind in xrange(2, len(lines)):
		line = [float(i) for i in lines[ind]]
		Ui.append(line)

	realXs, realTs, realUis = [], [], []

	for x in xrange(0, len(xs)):
		for t in xrange(0, len(ts)):
			realXs.append(xs[x])
			realTs.append(ts[t])
			realUis.append(Ui[t][x])

	return realXs, realTs, realUis

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

xs, ts, Ui = getDataFromFile()

ax.scatter(xs, ts, Ui, c='r', marker='o')

ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('U`')

plt.show()

