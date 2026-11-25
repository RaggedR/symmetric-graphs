import numpy as np, colorsys, itertools, sys
ip = int(sys.argv[1]); fname = sys.argv[2]; op=float(sys.argv[3]); lw=float(sys.argv[4]); title=sys.argv[5]
P=np.load('/tmp/e8_P.npy'); R=np.load('/tmp/e8_R.npy'); n=len(R); G=np.round(R@R.T).astype(int)
radii=np.round(np.hypot(P[:,0],P[:,1]),4); rings=sorted(set(radii)); ring_of={r:i for i,r in enumerate(rings)}
s=7.3/max(rings)
cols=[colorsys.hsv_to_rgb(0.02+0.80*i/7,0.75,1.0) for i in range(8)]
E=[(i,j) for i,j in itertools.combinations(range(n),2) if G[i,j]==ip]
o=[r"\documentclass[tikz,border=6pt]{standalone}",r"\usepackage{tikz}",r"\begin{document}",r"\begin{tikzpicture}"]
for i,(r,g,b) in enumerate(cols): o.append(fr"\definecolor{{ring{i}}}{{rgb}}{{{r:.3f},{g:.3f},{b:.3f}}}")
o.append(r"\definecolor{bg}{rgb}{0.04,0.05,0.09}")
o.append(r"\fill[bg] (-8.2,-8.9) rectangle (8.2,8.2);")
for i,j in E:
    ro=max(ring_of[radii[i]],ring_of[radii[j]])
    o.append(fr"\draw[ring{ro},line width={lw}pt,opacity={op}] ({P[i,0]*s:.3f},{P[i,1]*s:.3f})--({P[j,0]*s:.3f},{P[j,1]*s:.3f});")
for k in range(n):
    o.append(fr"\fill[ring{ring_of[radii[k]]}] ({P[k,0]*s:.3f},{P[k,1]*s:.3f}) circle (1.4pt);")
o.append(fr"\node[white!80!black,font=\small\itshape] at (0,-8.4) {{{title}}};")
o+=[r"\end{tikzpicture}",r"\end{document}"]
open(f'/Users/robin/git/graphviz/e8/{fname}.tex','w').write("\n".join(o))
print(fname, "edges:", len(E))
