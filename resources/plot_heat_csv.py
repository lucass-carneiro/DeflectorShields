import numpy as np
import pandas as pd

csv = pd.read_csv("C:/Users/steve/repos/DeflectorShields/plot.csv")
df_pivoted = csv.pivot(columns='x', index='y', values='val')

print(max(csv["x"]),min(csv["x"]),max(csv["y"]),min(csv["y"]))


import matplotlib.pyplot as plt
import seaborn as sns
ax = sns.heatmap(data=df_pivoted, yticklabels=40, xticklabels=40, annot=False, fmt='f', cmap='RdYlGn', cbar=True, cbar_kws={'label': 'val'}, square=True)
ax.tick_params(labelrotation=0)

plt.gca().invert_yaxis()
plt.savefig("C:/Users/steve/Downloads/deflector.png")
plt.show()

#plt.pcolor(csv["x"], csv["y"], csv["v"])
#plt.colorbar()
#plt.colorbar()
#plt.show()