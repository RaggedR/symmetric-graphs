import numpy as np, colorsys, itertools
P=np.load('/tmp/e7_P.npy'); W=np.load('/tmp/e7_W.npy'); n=len(W)
H=np.round(W@W.T).astype(int)
radii=np.round(np.hypot(P[:,0],P[:,1]),4); rings=sorted(set(radii))
ring_of={r:i for i,r in enumerate(rings)}
s=7.3/max(rings)
cols=[colorsys.hsv_to_rgb(0.03+0.72*i/(len(rings)-1),0.78,1.0) for i in range(len(rings))]
E=[(i,j) for i,j in itertools.combinations(range(n),2) if H[i,j]==1]   # Gosset graph, 27-regular
o=[r"\documentclass[tikz,border=6pt]{standalone}",r"\usepackage{tikz}",r"\begin{document}",r"\begin{tikzpicture}"]
for i,(r,g,b) in enumerate(cols): o.append(fr"\definecolor{{ring{i}}}{{rgb}}{{{r:.3f},{g:.3f},{b:.3f}}}")
o.append(r"\definecolor{bg}{rgb}{0.04,0.05,0.09}")
o.append(r"\fill[bg] (-8.2,-8.9) rectangle (8.2,8.2);")
for i,j in E:
    ro=max(ring_of[radii[i]],ring_of[radii[j]])
    o.append(fr"\draw[ring{ro},line width=0.12pt,opacity=0.30] ({P[i,0]*s:.3f},{P[i,1]*s:.3f})--({P[j,0]*s:.3f},{P[j,1]*s:.3f});")
for k in range(n):
    # nudge the 2 central points apart so both are visible
    x,y=P[k,0]*s,P[k,1]*s
    if radii[k]<1e-6: x+= (0.14 if k%2==0 else -0.14)
    o.append(fr"\fill[ring{ring_of[radii[k]]}] ({x:.3f},{y:.3f}) circle (1.7pt);")
o.append(r"\node[white!80!black,font=\small\itshape] at (0,-8.4) {$E_7$: 56 minuscule weights on the Coxeter plane (18-fold); the Gosset graph};")
o+=[r"\end{tikzpicture}",r"\end{document}"]
open('/Users/robin/git/graphviz/e8/e7_mandala.tex','w').write("\n".join(o))
print("E7 Gosset-graph edges:",len(E),"(expect 56*27/2 = 756)")
