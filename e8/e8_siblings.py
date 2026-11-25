import numpy as np, colorsys, itertools
P = np.load('/tmp/e8_P.npy'); R = np.load('/tmp/e8_R.npy'); n=len(R)
G = np.round(R@R.T).astype(int)
radii = np.round(np.hypot(P[:,0],P[:,1]),4); rings=sorted(set(radii))
ring_of={r:i for i,r in enumerate(rings)}
s = 3.3/max(rings)                       # panel radius ~3.3cm
cols=[colorsys.hsv_to_rgb(0.02+0.80*i/7,0.75,1.0) for i in range(8)]
panels=[(1,"ip = +1  (56-reg): the mandala",(-4.4, 4.4),0.22,0.08),
        (-1,"ip = -1  (56-reg): the weave",( 4.4, 4.4),0.22,0.08),
        (0,"ip = 0  (126-reg): the bloom",(-4.4,-4.4),0.09,0.07),
        (-2,"ip = -2  (120 K2): antipodal diameters",( 4.4,-4.4),0.60,0.22)]
def edges(ip): return [(i,j) for i,j in itertools.combinations(range(n),2) if G[i,j]==ip]
out=[r"\documentclass[tikz,border=6pt]{standalone}",r"\usepackage{tikz}",r"\begin{document}",r"\begin{tikzpicture}"]
for i,(r,g,b) in enumerate(cols): out.append(fr"\definecolor{{ring{i}}}{{rgb}}{{{r:.3f},{g:.3f},{b:.3f}}}")
out.append(r"\definecolor{bg}{rgb}{0.04,0.05,0.09}")
out.append(r"\fill[bg] (-8.6,-8.8) rectangle (8.6,8.6);")
for ip,label,(cx,cy),op,lw in panels:
    for i,j in edges(ip):
        ro=max(ring_of[radii[i]],ring_of[radii[j]])
        x1,y1=cx+P[i,0]*s,cy+P[i,1]*s; x2,y2=cx+P[j,0]*s,cy+P[j,1]*s
        out.append(fr"\draw[ring{ro},line width={lw}pt,opacity={op}] ({x1:.3f},{y1:.3f})--({x2:.3f},{y2:.3f});")
    for k in range(n):
        x,y=cx+P[k,0]*s,cy+P[k,1]*s
        out.append(fr"\fill[ring{ring_of[radii[k]]}] ({x:.3f},{y:.3f}) circle (0.7pt);")
    out.append(fr"\node[white!75!black,font=\scriptsize] at ({cx},{cy-3.9}) {{{label}}};")
out.append(r"\node[white!85!black,font=\small\itshape] at (0,8.1) {The four orbital graphs of $H(W(E_8),W(E_7))$ on the 240 roots};")
out+= [r"\end{tikzpicture}",r"\end{document}"]
open('/Users/robin/git/graphviz/e8/e8_siblings.tex','w').write("\n".join(out))
print("edge counts:", {ip:len(edges(ip)) for ip,_,_,_,_ in panels})
