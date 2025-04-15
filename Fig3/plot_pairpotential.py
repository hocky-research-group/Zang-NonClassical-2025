import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import argparse

parser=argparse.ArgumentParser()
filename=parser.add_argument("--filename",help="Filename (required)",type=str,required=True)
args = parser.parse_args()
locals().update(vars(args))

plt.figure(figsize=(8,6))
X, Y = [], []
for line in open(filename, 'r'):
    row=line.split()
    if(row[0]!='#'):
        values = [float(s) for s in row]
        X.append(values[0])
        Y.append(values[1])

plt.plot(X, Y)
plt.ylim(-30,30)
plt.xlabel(r"$r$"+" ("+r"$nm$"+")",fontsize=13)
plt.ylabel(r"$V$"+" ("+r"$k_{B}T$"+")",fontsize=13)
plt.show()
