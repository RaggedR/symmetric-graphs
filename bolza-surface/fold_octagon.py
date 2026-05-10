#!/usr/bin/env python3
"""
Fold the Bolza octagon into a double torus.

Each pentagon folds into a torus via bilinear map. The 48 triangles
get placed on two tori, with the gash visible between them.
"""

import numpy as np
import json
from itertools import product
from collections import defaultdict, deque

# ══════════════════════════════════════════════════════════════════════
#  Build tiling (same compact setup as before)
# ══════════════════════════════════════════════════════════════════════

def mm3(A,B):return(A@B)%3
def mk3(A):return tuple((A%3).flatten())
def md3(A):return(A[0,0]*A[1,1]-A[0,1]*A[1,0])%3
def mi3(A):
    d=md3(A);di=1 if d==1 else 2
    return(di*np.array([[A[1,1],(-A[0,1])%3],[(-A[1,0])%3,A[0,0]]]))%3
el={}
for a,b,c,d in product(range(3),repeat=4):
    A=np.array([[a,b],[c,d]])
    if md3(A)==0:continue
    k=mk3(A)
    if k not in el:el[k]=A%3
pm=list(el.values());pk=list(el.keys());ki={k:i for i,k in enumerate(pk)}
def gm(i,j):return ki[mk3(mm3(pm[i],pm[j]))]
def ginv(i):return ki[mk3(mi3(pm[i]))]
I=ki[mk3(np.eye(2,dtype=int))]
def order(i):
    c=i
    for n in range(1,49):
        if c==I:return n
        c=gm(c,i)
ords=[order(i) for i in range(48)]
ag=bg=cg=None
for a in range(48):
    if ords[a]!=2:continue
    for b in range(48):
        if ords[b]!=3:continue
        if ords[gm(a,b)]==8:ag,bg,cg=a,b,gm(a,b);break
    if ag:break

a_real=np.array([[0,1],[-1,0]],dtype=float)
cp8=np.cos(np.pi/8);sp3=np.sin(np.pi/3);cp3=np.cos(np.pi/3)
rhs=2*cp8/sp3;fp=(rhs+np.sqrt(rhs**2-4))/2
b_real=np.array([[cp3,-fp*sp3],[sp3/fp,cp3]])
b_inv_real=np.linalg.inv(b_real);c_real=a_real@b_real
def fxp(M):
    a,b,c,d=M[0,0],M[0,1],M[1,0],M[1,1];disc=4-(a+d)**2
    z1=((a-d)+1j*np.sqrt(disc))/(2*c);z2=((a-d)-1j*np.sqrt(disc))/(2*c)
    return z1 if z1.imag>0 else z2
z_a=fxp(a_real);z_b=fxp(b_real);z_c=fxp(c_real)
def mob(M,z):a,b,c,d=M[0,0],M[0,1],M[1,0],M[1,1];return(a*z+b)/(c*z+d)
def to_d(z):return(z-1j)/(z+1j)
wc=to_d(z_c)
def rc(w):return(w-wc)/(1-np.conj(wc)*w)
def mkk(M,tol=1e-6):return tuple(round(x/tol)*tol for x in M.flatten())
def hd(w):r=abs(w);return 2*np.arctanh(r) if r<1 else 1e9

print("Generating tiling...")
at=[];vis=set();q=deque([(np.eye(2),)])
while q and len(at)<5000:
    (M,)=q.popleft();m=mkk(M)
    if m in vis:continue
    vis.add(m)
    try:wa=rc(to_d(mob(M,z_a)));wb=rc(to_d(mob(M,z_b)));wcc=rc(to_d(mob(M,z_c)))
    except:continue
    if max(abs(wa),abs(wb),abs(wcc))>0.997:continue
    ar=abs((wb-wa).real*(wcc-wa).imag-(wb-wa).imag*(wcc-wa).real)
    if ar<1e-10:continue
    at.append((wa,wb,wcc))
    for gm2 in[a_real,b_real,b_inv_real]:q.append((M@gm2,))

td=[(i,hd((at[i][0]+at[i][1]+at[i][2])/3))for i in range(len(at))]
td.sort(key=lambda x:x[1])
fund=[at[idx] for idx,_ in td[:48]]

def ck(z):return(round(z.real,7),round(z.imag,7))
vm={};vl=[]
def vi(z):
    k=ck(z)
    if k not in vm:vm[k]=len(vl);vl.append(z)
    return vm[k]
it=[]
for wa,wb,wcc in fund:it.append((vi(wa),vi(wb),vi(wcc)))
N=len(vl)
print(f"  {N} vertices, 48 triangles")

# ══════════════════════════════════════════════════════════════════════
#  Octagon vertices, handle classification, bilinear maps
# ══════════════════════════════════════════════════════════════════════

# Octagon vertices
R_oct = max(abs(vl[i]) for i in range(N))
ov = sorted([i for i in range(N) if abs(abs(vl[i]) - R_oct) < 0.01],
            key=lambda i: np.angle(vl[i]))

# Handle classification: cut is V₀→V₄ (horizontal)
V0 = vl[ov[0]]; V4 = vl[ov[4]]
cut_n = 1j * (V4 - V0)
def wh(w): return 0 if ((w - V0) * np.conj(cut_n)).real >= 0 else 1

th = []
for ia, ib, ic in it:
    cen = (vl[ia] + vl[ib] + vl[ic]) / 3
    th.append(wh(cen))
print(f"  Handles: {sum(1 for h in th if h==0)} + {sum(1 for h in th if h==1)}")

# Bilinear maps: quad V₀V₁V₂V₃ for handle 0, V₄V₅V₆V₇ for handle 1
def bilinear_inv(w, Q):
    q = [np.array([c.real, c.imag]) for c in Q]
    p = np.array([w.real, w.imag])
    s, t = 0.5, 0.5
    for _ in range(30):
        fwd = (1-s)*(1-t)*q[0] + s*(1-t)*q[1] + s*t*q[2] + (1-s)*t*q[3]
        F = fwd - p
        if np.linalg.norm(F) < 1e-14: break
        dFs = -(1-t)*q[0] + (1-t)*q[1] + t*q[2] - t*q[3]
        dFt = -(1-s)*q[0] - s*q[1] + s*q[2] + (1-s)*q[3]
        J = np.column_stack([dFs, dFt])
        try: delta = np.linalg.solve(J, -F)
        except: break
        s += delta[0]; t += delta[1]
    return s, t

quads = {
    0: [vl[ov[0]], vl[ov[1]], vl[ov[2]], vl[ov[3]]],
    1: [vl[ov[4]], vl[ov[5]], vl[ov[6]], vl[ov[7]]]
}

# Verify
for h in [0, 1]:
    base = 0 if h == 0 else 4
    for k in range(4):
        s, t = bilinear_inv(vl[ov[base+k]], quads[h])
        print(f"  Handle {h}: V{base+k} → ({s:.3f}, {t:.3f})")

# ══════════════════════════════════════════════════════════════════════
#  Map to 3D
# ══════════════════════════════════════════════════════════════════════

print("\nMapping to 3D...")

R_tor = 1.3   # torus major radius
r_tor = 0.5   # torus minor radius
sep = 2.2     # separation

def torus_3d(u_frac, v_frac, handle):
    """Map (u,v) fractions to 3D. Wraps to [0,1] automatically."""
    u = 2 * np.pi * (u_frac % 1.0)
    v = 2 * np.pi * (v_frac % 1.0)
    offset = -sep if handle == 0 else sep
    x = offset + (R_tor + r_tor * np.cos(v)) * np.cos(u)
    y = r_tor * np.sin(v)
    z = (R_tor + r_tor * np.cos(v)) * np.sin(u)
    return np.array([x, y, z])

# Map each vertex using its handle's bilinear map
pos = np.zeros((N, 3))
assigned = np.zeros(N, dtype=bool)

for ti, (ia, ib, ic) in enumerate(it):
    h = th[ti]
    for vi2 in [ia, ib, ic]:
        if not assigned[vi2]:
            s, t = bilinear_inv(vl[vi2], quads[h])
            pos[vi2] = torus_3d(s, t, h)
            assigned[vi2] = True

print(f"  Positions computed for {assigned.sum()} / {N} vertices")
print(f"  x: [{pos[:,0].min():.2f}, {pos[:,0].max():.2f}]")
print(f"  y: [{pos[:,1].min():.2f}, {pos[:,1].max():.2f}]")
print(f"  z: [{pos[:,2].min():.2f}, {pos[:,2].max():.2f}]")

# Normals
normals = np.zeros_like(pos); counts = np.zeros(N)
for ia, ib, ic in it:
    fn = np.cross(pos[ib]-pos[ia], pos[ic]-pos[ia])
    n = np.linalg.norm(fn)
    if n > 1e-10: fn /= n
    for v2 in [ia, ib, ic]: normals[v2] += fn; counts[v2] += 1
for i in range(N):
    if counts[i] > 0:
        normals[i] /= counts[i]
        n = np.linalg.norm(normals[i])
        if n > 1e-10: normals[i] /= n

# Colors by face
def hsl(h2, s2, l):
    if s2 == 0: return (l, l, l)
    def f(p, q2, t):
        if t<0:t+=1
        if t>1:t-=1
        if t<1/6:return p+(q2-p)*6*t
        if t<1/2:return q2
        if t<2/3:return p+(q2-p)*(2/3-t)*6
        return p
    q2=l*(1+s2) if l<0.5 else l+s2-l*s2; p=2*l-q2
    return (f(p,q2,h2+1/3), f(p,q2,h2), f(p,q2,h2-1/3))

fg = defaultdict(list)
for ti, (ia, ib, ic) in enumerate(it): fg[ib].append(ti)
tc = [(0.5,0.5,0.5)] * len(it)
for fi, (fv, tl) in enumerate(fg.items()):
    h2 = (fi * 0.618034) % 1.0
    r2, g, b = hsl(h2, 0.65, 0.55)
    for ti in tl: tc[ti] = (r2, g, b)

# Gash curve data (precomputed)
gash_uv = [
    [0.0000, 0.0000], [-0.0707, 0.0679], [-0.1414, 0.1306], [-0.2121, 0.1887],
    [-0.2828, 0.2426], [-0.3536, 0.2929], [-0.4243, 0.3398], [-0.4950, 0.3837],
    [-0.5657, 0.4249], [-0.6364, 0.4636], [-0.7071, 0.5000], [-0.7778, 0.5343],
    [-0.8485, 0.5668], [-0.9192, 0.5975], [-0.9899, 0.6266], [-1.0607, 0.6542],
    [-1.1314, 0.6804], [-1.2021, 0.7054], [-1.2728, 0.7291], [-1.3435, 0.7518],
    [-1.4142, 0.7735]
]


# ══════════════════════════════════════════════════════════════════════
#  HTML
# ══════════════════════════════════════════════════════════════════════

print("\nGenerating HTML...")

html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>Bolza Surface — 48 Triangles on a Double Torus</title>
<style>
*{{margin:0;padding:0;box-sizing:border-box}}
body{{background:#08080c;overflow:hidden;font-family:'Courier New',monospace}}
#info{{position:absolute;top:12px;left:16px;color:#8899aa;font-size:13px;line-height:1.6;pointer-events:none;max-width:400px}}
#info h2{{color:#ccd0d8;font-size:16px;margin-bottom:4px;font-weight:normal}}
#info .n{{color:#e8b84d}}
</style></head><body>
<div id="info">
<h2>Bolza Surface — Folded Octagon</h2>
<span class="n">48</span> triangles on two wounded tori<br>
Each pentagon → torus via bilinear map<br>
<span style="color:#ff4444">Red</span> = gash (to be sewn together)
</div>
<script type="importmap">
{{"imports":{{"three":"https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js","three/addons/":"https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/"}}}}
</script>
<script type="module">
import*as THREE from'three';
import{{OrbitControls}}from'three/addons/controls/OrbitControls.js';

const V={json.dumps(pos.tolist())};
const NR={json.dumps(normals.tolist())};
const T={json.dumps([(int(a),int(b),int(c)) for a,b,c in it])};
const C={json.dumps(tc)};
const gashUV={json.dumps(gash_uv)};

const scene=new THREE.Scene();scene.background=new THREE.Color(0x08080c);
const cam=new THREE.PerspectiveCamera(40,innerWidth/innerHeight,0.1,100);
cam.position.set(0, 4, 9);
const ren=new THREE.WebGLRenderer({{antialias:true}});
ren.setSize(innerWidth,innerHeight);ren.setPixelRatio(Math.min(devicePixelRatio,2));
ren.toneMapping=THREE.ACESFilmicToneMapping;ren.toneMappingExposure=1.2;
document.body.appendChild(ren.domElement);
const ctl=new OrbitControls(cam,ren.domElement);
ctl.enableDamping=true;ctl.dampingFactor=0.05;ctl.autoRotate=true;ctl.autoRotateSpeed=0.4;

scene.add(new THREE.AmbientLight(0x556677,0.6));
let l;
l=new THREE.DirectionalLight(0xffeedd,1.2);l.position.set(5,8,5);scene.add(l);
l=new THREE.DirectionalLight(0x99aacc,0.5);l.position.set(-4,2,-3);scene.add(l);
l=new THREE.DirectionalLight(0x8899bb,0.3);l.position.set(0,-4,5);scene.add(l);

const R={R_tor}, r={r_tor}, sep={sep};

function torusPoint(u, v, h) {{
  const offset = h === 0 ? -sep : sep;
  return new THREE.Vector3(
    offset + (R + r*Math.cos(v)) * Math.cos(u),
    r * Math.sin(v),
    (R + r*Math.cos(v)) * Math.sin(u)
  );
}}

// ── Ghost tori (semi-transparent background) ──
for (let h = 0; h < 2; h++) {{
  const tg = new THREE.TorusGeometry(R, r, 48, 96);
  const tm = new THREE.Mesh(tg, new THREE.MeshPhongMaterial({{
    color: 0x3a4555, side: THREE.DoubleSide, shininess: 15,
    transparent: true, opacity: 0.25
  }}));
  tm.position.x = h === 0 ? -sep : sep;
  scene.add(tm);
}}

// ── Gash curves ──
for (let h = 0; h < 2; h++) {{
  const pts = [];
  for (let i = 0; i < gashUV.length - 1; i++) {{
    for (let k = 0; k < 8; k++) {{
      const t = k/8;
      const u = ((1-t)*gashUV[i][0] + t*gashUV[i+1][0]) * 2*Math.PI;
      const v = ((1-t)*gashUV[i][1] + t*gashUV[i+1][1]) * 2*Math.PI;
      pts.push(torusPoint(u, v, h));
    }}
  }}
  const last = gashUV[gashUV.length-1];
  pts.push(torusPoint(last[0]*2*Math.PI, last[1]*2*Math.PI, h));
  const curve = new THREE.CatmullRomCurve3(pts, false);
  scene.add(new THREE.Mesh(
    new THREE.TubeGeometry(curve, 160, 0.025, 8, false),
    new THREE.MeshPhongMaterial({{ color: 0xff3333, emissive: 0x441111 }})
  ));
}}

// ── The 48 triangles ──
const ps=[],ns=[],cs=[];
for(let i=0;i<T.length;i++){{
  const[a,b,c]=T[i];const[cr,cg,cb]=C[i];
  for(const idx of[a,b,c]){{ps.push(...V[idx]);ns.push(...NR[idx]);cs.push(cr,cg,cb)}}
}}
const geo=new THREE.BufferGeometry();
geo.setAttribute('position',new THREE.Float32BufferAttribute(ps,3));
geo.setAttribute('normal',new THREE.Float32BufferAttribute(ns,3));
geo.setAttribute('color',new THREE.Float32BufferAttribute(cs,3));
scene.add(new THREE.Mesh(geo,new THREE.MeshPhongMaterial({{
  vertexColors:true,side:THREE.DoubleSide,shininess:30,specular:new THREE.Color(0x445566)
}})));

// ── Wireframe ──
const wp=[];
for(const[a,b,c]of T)wp.push(...V[a],...V[b],...V[b],...V[c],...V[c],...V[a]);
const wg=new THREE.BufferGeometry();
wg.setAttribute('position',new THREE.Float32BufferAttribute(wp,3));
scene.add(new THREE.LineSegments(wg,new THREE.LineBasicMaterial({{
  color:0xeedd88,transparent:true,opacity:0.5
}})));

window.onresize=()=>{{cam.aspect=innerWidth/innerHeight;cam.updateProjectionMatrix();ren.setSize(innerWidth,innerHeight)}};
(function a(){{requestAnimationFrame(a);ctl.update();ren.render(scene,cam)}})();
</script></body></html>"""

with open('/Users/robin/git/graphviz/bolza-surface/bolza_folded.html', 'w') as f:
    f.write(html)
print("Done! → bolza_folded.html")
