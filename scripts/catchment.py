# CRC catchment code

    # this script takes cdc crc incidence data from nation and plots density
    # by state so that we can see where burden is

# libraries
import os
import pandas as pd
import matplotlib as mtl
import matplotlib.pyplot as plt

# variables
aneuploidy_home = "/scratch/chd5n/aneuploidy/"

# read data
os.chdir(aneuploidy_home)
x = pd.read_csv("crc-cdc-incidence.txt", sep="\t")

# prepare data
y = x.reindex(columns=['States', 'Age-Adjusted Rate'])
y = y.dropna()

# plot incidence
y.plot.kde()
plt.show()

# explore incidence
y.describe()
y.iloc[46,:]
y.iloc[48,:]
