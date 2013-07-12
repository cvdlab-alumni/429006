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

#sugar bowl
Su0 = CUBICHERMITE(S1)([[COS(PI/6),SIN(PI/12),0],[COS(PI/6),-SIN(PI/12),0],[TAN(PI/6),0,0],[0,TAN(PI/6),0]]);
Su1 = CUBICHERMITE(S1)([[0.2+COS(PI/6),SIN(PI/12),1.9],[0.2+COS(PI/6),(-SIN(PI/12)),1.9],[TAN(PI/6),0,0],[0,TAN(PI/6),0]]);

sb = CUBICHERMITE(S1)([[2*COS(PI/6),2*SIN(PI/12),1],[2*COS(PI/6),-2*SIN(PI/12),1],[2*TAN(PI/6),0,0],[0,2*TAN(PI/6),0]]);
sb1 = MAP(BEZIER(S2)([Su1,sb,Su0]))(dom2D)

#sugar bowl bottom
Su2 = CUBICHERMITE(S1)([[1*COS(PI/6),1*SIN(PI/12),-0.6],[1*COS(PI/6),-1*SIN(PI/12),-0.6],[1*TAN(PI/6),0,0],[0,1*TAN(PI/6),0]]);
Su22 = CUBICHERMITE(S1)([[1.3*COS(PI/6),1.5*SIN(PI/12),-0.5],[1.5*COS(PI/6),-1.5*SIN(PI/12),-0.5],[1.5*TAN(PI/6),0,0],[0,1.5*TAN(PI/6),0]]);
sb2 = MAP(BEZIER(S2)([Su0,Su22,Su2]))(dom2D)

partial_sugarbowl_form = STRUCT([sb1,sb2])

partial_sb_black = COLOR([0,0,0])(partial_sugarbowl_form)
partial_sb_white = COLOR([1,1,1])(R([1,2])(PI/6)(partial_sugarbowl_form))

partial_sb = STRUCT([partial_sb_white,partial_sb_black])

sugarbowl_form = STRUCT([partial_sb,R([1,2])(PI/3)]*6)

#knob
k0 = CUBICHERMITE(S2)([[-0.15,0,0],[0.15,0,0.5],[0,0,0.5],[0.5,0,0]])
k1 = CUBICHERMITE(S2)([[0,0,0],[0.15,0,0.3],[0,0,0.3],[0.3,0,0]])
k2 = CUBICHERMITE(S1)([[-0.15,0,0],[0,0,0],[0,-0.15,0],[0,0.15,0]])
k3 = CUBICHERMITE(S1)([[0.15,0,0.5],[0.15,0,0.3],[0,-0.15,0],[0,0.15,0]])
knob_right = MAP(COONSPATCH([k2,k3,k0,k1]))(dom2D);

k4 = CUBICHERMITE(S1)([[-0.15,0,0],[0,0,0],[0,0.15,0],[0,-0.15,0]])
k5 = CUBICHERMITE(S1)([[0.15,0,0.5],[0.15,0,0.3],[0,0.15,0],[0,-0.15,0]])
knob_left = MAP(COONSPATCH([k4,k5,k0,k1]))(dom2D);

knob = COLOR([0,0,0])(T([1,3])([0.1,2.6])(S([1,2,3])([0.36,0.36,0.36])(STRUCT([knob_left,knob_right]))))

domain_sphere = PROD([INTERVALS(PI/2)(24),INTERVALS(PI/6)(36)])
def partial_sphere (radius,domain):
	fx  = lambda p: radius * math.cos(p[0])  * math.sin  (p[1])
	fy  = lambda p: radius * math.cos(p[0]) * math.cos (p[1])
	fz  = lambda p: radius * math.sin(p[0]) 
	ret=  MAP([fx, fy, fz])(domain)
	return ret;


partial_cover = R([1,2])(PI/12)(partial_sphere(1.15,domain_sphere));
partial_cover_black = COLOR([0,0,0])(partial_cover)
partial_cover_white = COLOR([1,1,1])(R([1,2])(PI/6)(partial_cover))
partial_cover_bw = STRUCT([partial_cover_white,partial_cover_black])
partial_cv = STRUCT([partial_cover_bw,R([1,2])(PI/3)] *6)
cover = T([3])([1.95])(R([1,2])(PI/6)(S([3])([0.56])(partial_cv)))
cov = STRUCT([cover,knob])

def arc(alpha,r,R):
	domain = PROD([INTERVALS(alpha)(36),INTERVALS(R-r)(1)])
	def mapping(v):
		a = v[0];
		r = v[1];
		return [r*COS(a), r*SIN(a)];
	model = MAP(mapping)(domain);
	return model;


arc1 = arc(PI/6,0,1.0);
arc_white = COLOR([1,1,1])(R([1,2])(PI/12)(arc1))
arc_black = COLOR([0,0,0])(R([1,2])(-PI/12)(arc1))
partial_sb_bottom = STRUCT([arc_black,arc_white])
sb_bottom = STRUCT([partial_sb_bottom,R([1,2])(PI/3)]*6)
sb = T([3])([-0.6])(sb_bottom)

sugarbowl = STRUCT([sugarbowl_form,cov,sb])
VIEW(sugarbowl)