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

light_brown = COLOR([194/255, 115/255, 51/255])
brown = COLOR([150/255, 75/255, 0/255])


def table_top(c0,c1,c00,c11,w,h,s):
	cmap0 = MAP(CYLINDRICALSURFACE([c0,[0,0,h]]))(dom2D);
	cmap1 = MAP(CYLINDRICALSURFACE([c1,[0,0,h]]))(dom2D);
	cmap00 = MAP(CYLINDRICALSURFACE([c00,[0,0,h]]))(dom2D);
	cmap11 = MAP(CYLINDRICALSURFACE([c11,[0,0,h]]))(dom2D);
	cuboid = brown(T([1,2])([w-s,2*s])(CUBOID([s,6.1-1.2,h])))
	cuboid1 = brown(T([1])([0.6])(CUBOID([w-1.2,s,h])))
	cuboid2 = T([2])([6.1-s])(cuboid1)
	cuboid3 = light_brown(T([1,2,3])([0.6,0.3,h-0.11])(CUBOID([w-1.2,6.1-0.6,0.01])))
	bez1 = MAP(BEZIER(S2)([c0,c1]))(dom2D);
	bez2 = MAP(BEZIER(S2)([c00,c11]))(dom2D);
	bez3 = T([3])([h])(bez1);
	bez4 = T([3])([h])(bez2);
	bez5 = light_brown(T([3])([h-0.1])(MAP(BEZIER(S1)([c1,c11]))(dom2D)));
	border_r0 = brown(STRUCT([bez1,bez2,bez3,bez4,cmap0,cmap1,cmap00,cmap11,cuboid]))
	border_r = STRUCT([border_r0,bez5])
	border_l = T([1,2])([4,6.1])(R([1,2])(PI)(border_r))
	border = STRUCT([border_r,border_l,cuboid1,cuboid2,cuboid3])
	return border;

def trasla (p,v):
	q = []
	length=len(p)
	for i in range(length):
		q += [ADD([p[i],v])]
	return q;



def table(punti0,punti1,punti2,n,w,d,s):
	tavoli = []
	for i in range(n):
		cyl = CYLINDER([0.1,6.8-i*s])(24)
		p0 = trasla(punti0,[+i*s,+i*s,0])
		p1 = trasla(punti1,[-i*s,+i*s,0])
		p2 = trasla(punti2,[-i*s,+i*s,0])
		c0 = CUBICHERMITE(S1)([p0[0],p0[1],[-s*2,0,0],[s*2,0,0]]);
		c1 = CUBICHERMITE(S1)([p1[0],p1[1],[s*4,0,0],[0,s*4,0]]);
		c2 = CUBICHERMITE(S1)([p2[0],p2[1],[s*2,0,0],[0,s*2,0]]);
		p00[0] = [punti0[0]+i*s,d-(punti0[1]+i*s),punti0[2]]
		p00[1] = [punti0[0]+i*s,d-(punti0[1]+i*s),punti0[2]]
		p11[0] = [punti1[0]-i*s,d-(punti1[1]+i*s),punti1[2]]
		p11[1] = [punti1[0]-i*s,d-(punti1[1]+i*s),punti1[2]]
		p22[0] = [punti2[0]-i*s,d-(punti2[1]+i*s),punti2[2]]
		p22[1] = [punti2[0]-i*s,d-(punti2[1]+i*s),punti2[2]]
		p3 = trasla(p0,[i*s,+i*s,+s])
		p4 = trasla(p00,[i*s,i*s,s])
		c3 = CUBICHERMITE(S1)([p00[0],p00[1],[-s*2,0,0],[s*2,0,0]]);
		c4 = CUBICHERMITE(S1)([p11[0],p11[1],[s*4,0,0],[0,-s*4,0]]);
		c5 = CUBICHERMITE(S1)([p22[0],p22[1],[s*2,0,0],[0,-s*2,0]]);
		c00 = BEZIER(S1)(p0);
		c33 = BEZIER(S1)(p00);
		cmap0 = MAP(CYLINDRICALSURFACE(c0,[0,0,s]))(dom2D);
		cmap1 = MAP(CYLINDRICALSURFACE(c1,[0,0,s]))(dom2D);
		cmap2 = MAP(CYLINDRICALSURFACE(c2,[0,0,s]))(dom2D);
		cmap3 = MAP(CYLINDRICALSURFACE(c3,[0,0,s]))(dom2D);
		cmap4 = MAP(CYLINDRICALSURFACE(c4,[0,0,s]))(dom2D);
		cmap5 = MAP(CYLINDRICALSURFACE(c5,[0,0,s]))(dom2D);
		bez1 = MAP(BEZIER(S2)([c1,c2]))(dom2D);
		bez2 = MAP(BEZIER(S2)([c4,c5]))(dom2D);
		bez3 = T([3])([s])(bez1);
		bez4 = T([3])([s])(bez2);
		bez5 = MAP(BEZIER(S2)([c0,c00]))(dom2D);
		bez6 = MAP(BEZIER(S2)([c3,c33]))(dom2D);
		bez7 = T([3])([s])(bez5);
		bez8 = T([3])([s])(bez6);
		cuboid0 = T([1,2])([i*s+s,i*s])(CUBOID([p1[0][0]-p0[0][0],s,s]))
		cuboid1 = T([2])([p00[0][1]-i*s-s])(cuboid0)
		cuboid2 = T([1,2])([w-i*s-s,i*s+2*s])(CUBOID([s,p22[1][1]-p2[1][1],s]))
		basis = brown(STRUCT([cmap0,cmap1,cmap2,cmap3,cmap4,cmap5,bez1,bez2,bez3,bez4,bez5,bez6,bez7,bez8,cuboid0,cuboid1,cuboid2]))
		#table leg
		leg1 = T([1,2])([i*s+s+0.1,i*s+0.15])(cyl)
		leg2 = T([1])([p1[0][0]-p0[0][0]])(leg1)
		leg12 = STRUCT([leg1,leg2])
		leg4 = T([2])([p00[0][1]-i*s-s])(leg12)
		table_leg = brown(STRUCT([leg12,leg4]))
		#
		if i!=0:
			cuboid3 = brown(T([1,2,3])([i*s+s,i*s,6.8-((i-1)*s)-s])(CUBOID([p1[0][0]-p0[0][0]+0.2,0.2,0.3])))
			cuboid4 = brown(T([2])([p00[0][1]-i*s-s])(cuboid3))
			cuboid45 = T([1,2,3])([p1[0][0]-p0[0][0]+i*s+s+0.1,i*s,6.8-((i-1)*s)-s])(CUBOID([0.2,p11[0][1]-p1[0][1]-0.1,0.3]))
			cuboid55 = T([1,2,3])([i*s+s-0.1,i*s,6.8-((i-1)*s)-0.15])(CUBOID([0.2,p11[0][1]-p1[0][1]-0.1,0.15]))
			cuboid5 = brown(STRUCT([cuboid45,cuboid55]))
			cuboid6 = light_brown(T([1,2,3])([i*s+s+0.1,i*s+0.2,6.8-((i-1)*s)-0.1])(CUBOID([p1[0][0]-p0[0][0],p11[0][1]-p1[1][1]+0.1,0.05])))
		else:
			cuboid3 = brown(T([1,2,3])([i*s+s,i*s+0.1,5.9])(CUBOID([p1[0][0]-p0[0][0]+0.05,0.2,0.9])))
			cuboid4 = brown(T([2])([p00[0][1]-i*s-s-0.05])(cuboid3))
			cuboid5 = brown(T([1,2,3])([p1[0][0]-p0[0][0]+i*s+s,i*s+0.15,5.9])(CUBOID([0.2,p11[0][1]-p1[0][1]-0.2,0.9])))
			cuboid6 = T([3])([6.8])(table_top(c1,c2,c4,c5,w,0.6,s))
		cub = STRUCT([cuboid3,cuboid4,cuboid5,cuboid6])
		table = STRUCT([basis,table_leg,cub])
		#tavoli += [table]
	return cuboid6;

p1 = [[0.3,0,0],[0.3,0.3,0]]
p2 = [[3.4,0,0],[4,0.6,0]]
p3 = [[3.4,0.3,0],[3.7,0.6,0]]
t = table(p1,p2,p3,4,4,6.1,0.3)

function sphere(r){
	return function(v){
	return [r*SIN(v[0])*COS(v[1]),r*SIN(v[0])*SIN(v[1]),r*COS(v[0])]}
}

function spheres(r,n){
	var domain_sfera = DOMAIN([[0,2*PI],[0,2*PI]])([24,36]);
	var mapping = sphere(r);
	var model = MAP(mapping)(domain_sfera);
	var model1 = model;
	for (i=1; i<n;i++){
		s = T([0])([i*2*r])(model)
		model1 = STRUCT([s,model1])
	}
	return brown(model1);
}
var spheres_r= T([0,1,2])([1.4,-0.2,7.1])(spheres(0.2,4))
var spheres_l = T([0,1,2])([1.4,6.3,7.1])(spheres(0.2,4))
var sp = STRUCT([spheres_l,spheres_r])
var table = STRUCT([t[0],t[1],t[2],t[3],sp])
DRAW(table)

