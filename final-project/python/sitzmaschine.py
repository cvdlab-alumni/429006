from pyplasm import *
import scipy
from scipy import *

#---------------------------------------------------------
def VERTEXTRUDE((V,coords)):
    """
        Utility function to generate the output model vertices in a 
        multiple extrusion of a LAR model.
        V is a list of d-vertices (each given as a list of d coordinates).
        coords is a list of absolute translation parameters to be applied to 
        V in order to generate the output vertices.
        
        Return a new list of (d+1)-vertices.
    """
    return CAT(AA(COMP([AA(AR),DISTR]))(DISTL([V,coords])))

def cumsum(iterable):
    # cumulative addition: list(cumsum(range(4))) => [0, 1, 3, 6]
    iterable = iter(iterable)
    s = iterable.next()
    yield s
    for c in iterable:
        s = s + c
        yield s

def larExtrude(model,pattern):
    V,FV = model
    d = len(FV[0])
    offset = len(V)
    m = len(pattern)
    outcells = []
    for cell in FV:
        # create the indices of vertices in the cell "tube"
        tube = [v + k*offset for k in range(m+1) for v in cell]
        # take groups of d+1 elements, via shifting by one
        rangelimit = len(tube)-d
        cellTube = [tube[k:k+d+1] for k in range(rangelimit)]
        outcells += [scipy.reshape(cellTube,newshape=(m,d,d+1)).tolist()]
    outcells = AA(CAT)(TRANS(outcells))
    outcells = [group for k,group in enumerate(outcells) if pattern[k]>0 ]
    coords = list(cumsum([0]+(AA(ABS)(pattern))))
    outVerts = VERTEXTRUDE((V,coords))
    newModel = outVerts, CAT(outcells)
    return newModel

def GRID(args):
    model = ([[]],[[0]])
    for k,steps in enumerate(args):
        model = larExtrude(model,steps*[1])
    V,cells = model
    verts = AA(list)(scipy.array(V) / AA(float)(args))
    return MKPOL([verts, AA(AA(lambda h:h+1))(cells), None])

#---------------------------------------------------------

dom1D = INTERVALS(1)(32)
dom2D = GRID([20,20])

c1 = CUBICHERMITE(S1)([[5,5.5,5],[4.85,8,2.5],[0,4,0],[0,0,-4]])
c2 = CUBICHERMITE(S1)([[4.85,8,2.5],[4.7,5.5,0],[0,0,-4],[0,-4,0]])
c3 = CUBICHERMITE(S1)([[5,5.5,4.8],[4.85,7.8,2.5],[0,3.4,0],[0,0,-3.4]])
c4 = CUBICHERMITE(S1)([[4.85,7.8,2.5],[4.7,5.5,0.2],[0,0,-3.4],[0,-3.4,0]])
c11 = CUBICHERMITE(S1)([[4.5,5.5,5],[4.5,8,2.5],[0,4,0],[0,0,-4]])
c22 = CUBICHERMITE(S1)([[4.5,8,2.5],[4.5,5.5,0],[0,0,-4],[0,-4,0]])
c33 = CUBICHERMITE(S1)([[4.5,5.5,4.8],[4.5,7.8,2.5],[0,3.4,0],[0,0,-3.4]])
c44 = CUBICHERMITE(S1)([[4.5,7.8,2.5],[4.5,5.5,0.2],[0,0,-3.4],[0,-3.4,0]])
out1 = MAP(BEZIER(S2)([c1,c11]))(dom2D);
out2 = MAP(BEZIER(S2)([c2,c22]))(dom2D);
out3 = MAP(BEZIER(S2)([c3,c33]))(dom2D);
out4 = MAP(BEZIER(S2)([c4,c44]))(dom2D);
out5 = MAP(BEZIER(S2)([c1,c3]))(dom2D);
out6 = MAP(BEZIER(S2)([c2,c4]))(dom2D);
out7 = MAP(BEZIER(S2)([c11,c33]))(dom2D);
out8 = MAP(BEZIER(S2)([c44,c22]))(dom2D);
s = STRUCT ([out1,out2,out3,out4,out5,out6,out7,out8])

c1b = CUBICHERMITE(S1)([[-0.5,5.5,5],[-0.35,8,2.5],[0,4,0],[0,0,-4]]);
c2b = CUBICHERMITE(S1)([[-0.35,8,2.5],[-0.2,5.5,0],[0,0,-4],[0,-4,0]]);
c3b = CUBICHERMITE(S1)([[-0.5,5.5,4.8],[-0.35,7.8,2.5],[0,3.4,0],[0,0,-3.4]]);
c4b = CUBICHERMITE(S1)([[-0.35,7.8,2.5],[-0.2,5.5,0.2],[0,0,-3.4],[0,-3.4,0]]);
c11b = CUBICHERMITE(S1)([[0,5.5,5],[0,8,2.5],[0,4,0],[0,0,-4]]);
c22b = CUBICHERMITE(S1)([[0,8,2.5],[0,5.5,0],[0,0,-4],[0,-4,0]]);
c33b = CUBICHERMITE(S1)([[0,5.5,4.8],[0,7.8,2.5],[0,3.4,0],[0,0,-3.4]]);
c44b = CUBICHERMITE(S1)([[0,7.8,2.5],[0,5.5,0.2],[0,0,-3.4],[0,-3.4,0]]);
out1b = MAP(BEZIER(S2)([c11b,c1b]))(dom2D);
out2b = MAP(BEZIER(S2)([c22b,c2b]))(dom2D);
out3b = MAP(BEZIER(S2)([c33b,c3b]))(dom2D);
out4b = MAP(BEZIER(S2)([c44b,c4b]))(dom2D);
out5b = MAP(BEZIER(S2)([c3b,c1b]))(dom2D);
out6b = MAP(BEZIER(S2)([c4b,c2b]))(dom2D);
out7b = MAP(BEZIER(S2)([c11b,c33b]))(dom2D);
out8b = MAP(BEZIER(S2)([c44b,c22b]))(dom2D);

b = STRUCT ([out1b,out2b,out3b,out4b,out5b,out6b,out7b,out8b])

cuboid = CUBOID([0.5,5.5,0.2])
ct = T([1,3])([4.5,4.8])(cuboid)
ctb = T([1,3])([-0.5,4.8])(cuboid)

cuboid1 = CUBOID([0.2,4.95,0.2])
ct1 = T([1,2])([4.5,0.55])(cuboid1)
ct1b = T([1,2])([-0.2,0.55])(cuboid1)

c5 = CUBICHERMITE(S1)([[4.7,0.15,0.4],[4.7,0.55,0],[0,0,-0.6],[0,0.6,0]]);
c6 = CUBICHERMITE(S1)([[4.7,0.35,0.4],[4.7,0.55,0.2],[0,0,-0.2],[0,0.2,0]]);

out9 = MAP(CYLINDRICALSURFACE([c5,[-0.2,0,0]]))(dom2D);
out10 = MAP(CYLINDRICALSURFACE([c6,[-0.2,0,0]]))(dom2D);
out11 = MAP(BEZIER(S2)([c5,c6]))(dom2D);
out12 = T([1])([-0.2])(out11);
s1 = STRUCT([out9,out10,out11,out12])

cuboid2 = CUBOID([0.2,0.2,4.4])
ct2 = T([1,2,3])([4.5,0.15,0.4])(cuboid2)

cs = STRUCT([ct2,s1,ct1])
csb = T([1])([-4.7])(cs)

#sfera
#domain_sfera = DOMAIN([[0,2*PI],[0,2*PI]])([24,36]);

model = SPHERE(0.2)([24,36]);
mt = T([1,2,3])([4.85,0.23,4.6])(model)
mt1 = T([1,2,3])([-0.35,0.23,4.6])(model)
mt2 = T([1,2,3])([4.25,0.23,2.1])(model)
mt3 = T([1,2,3])([0.2,0.23,2.1])(model)

struct = COLOR([184.0/255, 115.0/255, 51.0/255])(STRUCT([cs,csb,ctb,ct,s,b,mt,mt1,mt2,mt3]))

SIMPLEX_GRID = COMP([INSR(PROD),AA(QUOTE)])

#rifinitura cornicione

def wall_y_with_hole(dx,dy,dz,hole_y,hole_z,hole_dy,hole_dz):
	g1 = SIMPLEX_GRID (([[dx],[hole_y,-hole_dy, dy-hole_y-hole_dy],[dz]]))
	g2 = SIMPLEX_GRID ([[dx],[-hole_y, hole_dy],[hole_z, -hole_dz,dz-hole_z-hole_dz]])
	g = STRUCT([g1,g2])
	return g;



def wall_y_with_N_hole(dx,dy,dz,hole_y,hole_z,hole_dy,hole_dz,n):
	f = wall_y_with_hole(dx,dy,dz,hole_y,hole_z,hole_dy,hole_dz)
	t = f
	for i in range(n):
		t = STRUCT([t,T([2])([i*dy])(f)])
	return t;

f1 = wall_y_with_N_hole(0.1,0.4,3,0.1,1,0.2,1.5,5)
f2 = wall_y_with_N_hole(0.1,0.4,1.4,0.1,0.5,0.2,0.2,5)
f1t = T([1,2,3])([-0.15,2.2,1.7])(f1)
f2t = T([1,2,3])([-0.15,2.2,0.3])(f2)



def semicerchio(p1,p2,vettore):
	var1 = CUBICHERMITE(S1)(p1)
	var2 = CUBICHERMITE(S1)(p2)
	v12 = MAP(BEZIER(S2)([var1,var2]))(dom2D)
	v21 = T([1,2,3])(vettore)(v12)
	var1l = MAP(CYLINDRICALSURFACE([var1,vettore]))(dom2D)
	var2l = MAP(CYLINDRICALSURFACE([var2,vettore]))(dom2D)
	s = STRUCT([v12,v21,var1l,var2l])
	return s;



p0 = [[0,1.5,0.4],[0,1.7,0.2],[0,0,-0.3],[0,0.3,0]]
p00 = [[0,1.6,0.4],[0,1.7,0.3],[0,0,-0.2],[0,0.2,0]]
sm0 = semicerchio(p0,p00,[-0.2,0,0])
p1 = [[0,1.5,4.6],[0,1.7,4.8],[0,0,0.3],[0,0.3,0]]
p11 = [[0,1.6,4.6],[0,1.7,4.7],[0,0,0.2],[0,0.2,0]]
sm1 = semicerchio(p1,p11,[-0.2,0,0])
t1 = CUBICHERMITE(S1)([[-0.05,1.6,0.4],[-0.05,1.7,0.3],[0,0,-0.2],[0,0.2,0]]);
t0 = CUBICHERMITE(S1)([[-0.05,1.6,4.6],[-0.05,1.7,4.7],[0,0,0.2],[0,0.2,0]]);
t10 = MAP(BEZIER(S2)([t1,t0]))(dom2D)
t01 = T([1])([-0.1])(t10)
cx =  CUBOID([0.2,0.1,4.2])
cxt = T([1,2,3])([-0.2,1.5,0.4])(cx)
cx1 = CUBOID([0.1,0.5,4.4])
cx1t = T([1,2,3])([-0.15,1.7,0.3])(cx1)
rif = STRUCT([sm0,sm1,cxt,cx1t, t10,t01])
sr = R([1,2])(PI)(rif)
srt = T([1,2])([-0.2,6.4])(sr)
cx2 = CUBOID([0.2,3,0.1])
cx2t = T([1,2,3])([-0.2,1.7,0.2])(cx2)
cx2tb = T([1,2,3])([-0.2,1.7,4.7])(cx2)

struct_rif = STRUCT([rif,srt,cx2t,cx2tb,f1t,f2t])
rif2 = T([1])([4.7])(struct_rif)
rifinitura = STRUCT([struct_rif,rif2])

sedia1 = COLOR([184.0/255, 115.0/255, 51.0/255])(STRUCT([rifinitura,struct]))

#schienale

p2 = [[0,5.5,3],[0.4,5.5,2.6],[0,0,-0.8],[0.8,0,0]]
p22 = [[0.2,5.5,3],[0.4,5.5,2.8],[0,0,-0.4],[0.4,0,0]]
sm2 = semicerchio(p2,p22,[0,0.2,0])
p3 = [[0,5.5,9.7],[0.4,5.5,10.1],[0,0,0.8],[0.8,0,0]]
p33 = [[0.2,5.5,9.7],[0.4,5.5,9.9],[0,0,0.4],[0.4,0,0]]
sm3 = semicerchio(p3,p33,[0,0.2,0])
cil = T([2,3])([5.5,3])(CUBOID([0.2,0.2,6.7]))
sc = STRUCT([sm2,sm3,cil])
scr = T([1,3])([4.3,12.7])(R([1,3])(PI)(sc))
cil1 = T([1,2,3])([0.3,5.5,2.6])(CUBOID([3.7,0.2,0.2]))
cil1t = T([3])([7.3])(cil1)

schienale_ext = STRUCT([sc,scr,cil1,cil1t])

#parte interna

p3 = [[0.75,5.5,3.1],[1.05,5.5,2.8],[0,0,-0.6],[0.6,0,0]]
p33 = [[0.9,5.5,3.1],[1.05,5.5,2.95],[0,0,-0.3],[0.3,0,0]]
sm3 = semicerchio(p3,p33,[0,0.2,0])
p4 = [[0.75,5.5,9.6],[1.05,5.5,9.9],[0,0,0.6],[0.6,0,0]]
p44 = [[0.9,5.5,9.6],[1.05,5.5,9.75],[0,0,0.3],[0.3,0,0]]
sm4 = semicerchio(p4,p44,[0,0.2,0])
cil2 = T([1,2,3])([0.75,5.5,3])(CUBOID([0.15,0.2,6.7]))
ti1 = CUBICHERMITE(S1)([[0.9,5.575,3.1],[1.05,5.575,2.95],[0,0,-0.3],[0.3,0,0]]);
ti0 = CUBICHERMITE(S1)([[0.9,5.575,9.6],[1.05,5.575,9.75],[0,0,0.3],[0.3,0,0]]);
ti10 = MAP(BEZIER(S2)([ti1,ti0]))(dom2D)
ti01 = T([2])([0.05])(ti10)
sc1 = STRUCT([sm3,sm4,cil2,ti10,ti01])
scr1 = T([1,3])([4.3,12.7])(R([1,3])(PI)(sc1))
cil3 = T([1,2,3])([1.05,5.5,2.8])(CUBOID([2.4,0.2,0.15]))
cil3t = T([3])([6.95])(cil3)
interno1 = SIMPLEX_GRID([[-1.05,2.4],[-5.575,0.05],
	[-2.8,0.375,-0.5,0.23,-0.5,0.23,-0.5,0.23,-0.5,0.23,-0.5,0.23,-0.5,0.23,-0.5,0.23,-0.5,0.23,-0.5,0.385]]) 
interno2 = SIMPLEX_GRID([[-1.05,0.5,-0.5,0.4,-0.5,0.5],[-5.575,0.05],[-2.8,7.1]]) 
schienale_int = STRUCT([sc1,scr1,cil3,cil3t,interno1,interno2])
schienale = COLOR([184.0/255, 115.0/255, 51.0/255])(T([2,3])([-0.4,2.1])(R([2,3])(-PI/9)(STRUCT([schienale_int,schienale_ext]))))

#sedile

sedile = CUBOID([4.5,5.5,0.3])
sed_t = T([3])([2.3])(sedile)

l1 = CUBICHERMITE(S1)([[0,0.8,2.3],[0.8,0,2.3],[0,-1.6,0],[1.6,0,0]])
l2 = CUBICHERMITE(S1)([[0,4.7,2.3],[0.8,5.5,2.3],[0,1.6,0],[1.6,0,0]])
l3 = BEZIER(S1)([[0,0.8,2.3],[0,4.7,2.3]])
v12 = MAP(BEZIER(S2)([l1,l2]))(dom2D)
v21 = T([1,2,3])([0,0,0.3])(v12)
var3l = MAP(CYLINDRICALSURFACE([l2,[0,0,0.3]]))(dom2D)
var2l = MAP(CYLINDRICALSURFACE([l2,[0,0,0.3]]))(dom2D)
var1l = MAP(CYLINDRICALSURFACE([l1,[0,0,0.3]]))(dom2D)
sx = STRUCT([v12,v21,var1l,var2l,var3l])
srt = T([1,2])([4.5,5.5])(R([1,2])(PI)(sx))
centro = T([1,3])([0.8,2.3])(CUBOID([2.9,5.5,0.3]))
sedile = COLOR([184.0/255, 115.0/255, 51.0/255])(STRUCT([centro, srt,sx]))

#funghetti

domain_semisphere = PROD([INTERVALS(PI)(24),INTERVALS(PI)(36)])
def partial_sphere (radius,domain):
	fx  = lambda p: radius * math.cos(p[0])  * math.sin  (p[1])
	fy  = lambda p: radius * math.cos(p[0]) * math.cos (p[1])
	fz  = lambda p: radius * math.sin(p[0]) 
	ret=  MAP([fx, fy, fz])(domain)
	return ret

model1 = T([3])([0.1])(partial_sphere(0.15,domain_semisphere))
d = CYLINDER([0.1,0.1])(24);
d1 = T([3])([0.09])(CYLINDER([0.15,0.01])(24));
fung = STRUCT([model1,d,d1]);

#da vedere ancora

dt1 = T([1,2,3])([-0.25,5.5,5])(fung)
dt2 = T([2,3])([0.1,0.75])(R([2,3])(-PI/20)(dt1))
dt3 = T([1,2,3])([0.02,0.2,1.4])(R([2,3])(-2*PI/20)(dt1))
dt4 = T([1,2,3])([0.03,0.285,2])(R([2,3])(-3*PI/20)(dt1))
dt5 = T([1,2,3])([0.04,0.45,2.45])(R([2,3])(-4*PI/20)(dt1))
funghetti_left = STRUCT([dt1,dt2,dt3,dt4,dt5])
funghetti_right = T([1])([4.95])(funghetti_left)

funghetti = COLOR([184.0/255, 115.0/255, 51.0/255])(STRUCT([funghetti_left,funghetti_right]))

cilindro = CYLINDER([0.07,6])(24)
cilindro_r = R([1,3])(PI/2)(cilindro)
cilindro_t = COLOR([0,0,0])(T([1,2,3])([5.3,6.8,4.6])(cilindro_r))
sfera1 = COLOR([184.0/255, 115.0/255, 51.0/255])(T([1,2,3])([-0.7,6.8,4.6])(model))
sfera2 = COLOR([184.0/255, 115.0/255, 51.0/255])(T([1,2,3])([5.3,6.8,4.6])(model))
cil_sf = STRUCT([cilindro_t,sfera2,sfera1])

#cerniera

cil4 = CYLINDER([0.1,0.6])(24)
cil4r = R([1,3])(PI/2)(cil4)
cil4t1 = T([1,2,3])([1.4,5.6,2.65])(cil4r)
cil4t2 = T([1,2,3])([3.5,5.6,2.65])(cil4r)
cerniera = COLOR([0.3,0.3,0.3])(STRUCT([cil4t2,cil4t1]))
VIEW(cerniera)
#ruote

domain_wheel = PROD([INTERVALS(PI/2)(24),INTERVALS(5*PI/8-3*PI/8)(36)])
model1 = COLOR([184/255, 115/255, 51/255])(partial_sphere(0.3,domain_wheel))
disk = CYLINDER([0.15,0.2])(24)
disk_r = R([1,3])(PI/2)(disk)
disk_t = T([1,2,3])([0.1,0.08,0.08])(disk_r)
ruota = STRUCT([COLOR([0.3,0.3,0.3])(disk_t),T([1])([-0.1])(model1)])
ruota1 = T([1,2])([-0.1,6.5])(ruota)
ruota2 = T([1,2])([4.6,6.5])(ruota)
ruote = STRUCT([ruota2,ruota1])
sedia = STRUCT([sedia1,sedile,schienale,funghetti,cil_sf,cerniera,ruote])
VIEW(sedia)