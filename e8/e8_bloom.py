import numpy as np, colorsys, itertools
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
P=np.load('/tmp/e8_P.npy'); R=np.load('/tmp/e8_R.npy'); n=len(R); G=np.round(R@R.T).astype(int)
radii=np.round(np.hypot(P[:,0],P[:,1]),4); rings=sorted(set(radii)); ring_of={r:i for i,r in enumerate(rings)}
cols=[colorsys.hsv_to_rgb(0.02+0.80*i/7,0.75,1.0) for i in range(8)]
E=[(i,j) for i,j in itertools.combinations(range(n),2) if G[i,j]==0]
segs=[[P[i],P[j]] for i,j in E]
ecol=[cols[max(ring_of[radii[i]],ring_of[radii[j]])] for i,j in E]
fig,ax=plt.subplots(figsize=(8,8)); fig.patch.set_facecolor('#0a0d17'); ax.set_facecolor('#0a0d17')
ax.add_collection(LineCollection(segs,colors=ecol,linewidths=0.18,alpha=0.10))
ax.scatter(P[:,0],P[:,1],s=6,c=[cols[ring_of[r]] for r in radii],zorder=3)
ax.set_aspect('equal'); ax.axis('off')
ax.set_title("E_8 orbital ip=0: the bloom (126-regular, 15120 edges)",color='#cccccc',style='italic',fontsize=11,y=-0.06)
m=max(rings)*1.05; ax.set_xlim(-m,m); ax.set_ylim(-m,m)
fig.savefig('/Users/robin/git/graphviz/e8/e8_ip0.pdf',facecolor=fig.get_facecolor(),bbox_inches='tight')
print("bloom rendered:", len(E), "edges")
