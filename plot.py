import matplotlib.pyplot as plt
import numpy as np
import sys, csv

# get file arguments
files = sys.argv
dataFile = files[1]

# get data from file
with open(dataFile, newline='') as f:
    reader = csv.reader(f)
    dataMat = np.array(list(reader))

dataMat = dataMat[1:]
# convert dataMat to format [gene, PSM, #Tyr/PSM]
genes = [i[2] for i in dataMat]
PSMs = [float(i[6]) for i in dataMat]
tyrPSMs = [float(i[3])/float(i[6]) for i in dataMat]

r = np.linspace(0,1,len(genes))
theta = 2 * np.pi * r
dotsize = PSMs
distance = tyrPSMs

fig = plt.figure()
ax = fig.add_subplot(projection="polar",facecolor="whitesmoke")
ax.tick_params(grid_color="lightgray", grid_linestyle="--")
ax.scatter(theta, r, s=dotsize, label="data")

angle = np.deg2rad(67.5)
ax.legend(loc="lower left",
          bbox_to_anchor=(.5 + np.cos(angle)/2, .5 + np.sin(angle)/2))

fig.savefig("data.svg", bbox_inches='tight', format='svg')
plt.show()
