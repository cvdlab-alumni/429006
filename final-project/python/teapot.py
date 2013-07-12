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

#teapot's form
Su0 = CUBICHERMITE(S1)([[COS(PI/6),SIN(PI/12),0],[COS(PI/6),-SIN(PI/12),0],[TAN(PI/6),0,0],[0,TAN(PI/6),0]]);
Su1 = CUBICHERMITE(S1)([[0.2+COS(PI/6),SIN(PI/12),1.9],[0.2+COS(PI/6),(-SIN(PI/12)),1.9],[TAN(PI/6),0,0],[0,TAN(PI/6),0]]);
sb = CUBICHERMITE(S1)([[2*COS(PI/6),2*SIN(PI/12),1],[2*COS(PI/6),-2*SIN(PI/12),1],[2*TAN(PI/6),0,0],[0,2*TAN(PI/6),0]]);
tp1 = MAP(BEZIER(S2)([Su1,sb,Su0]))(dom2D)

#top
Su2 = CUBICHERMITE(S1)([[COS(PI/6),SIN(PI/12),3.1],[COS(PI/6),-SIN(PI/12),3.1],[TAN(PI/6),0,0],[0,TAN(PI/6),0]]);
Su3 = CUBICHERMITE(S1)([[0.1+COS(PI/6),SIN(PI/12),2.1],[0.1+COS(PI/6),(-SIN(PI/12)),2.1],[TAN(PI/6),0,0],[0,TAN(PI/6),0]]);
Su23 = CUBICHERMITE(S1)([[1.5*COS(PI/6),1.5*SIN(PI/12),2.6],[1.5*COS(PI/6),(-1.5*SIN(PI/12)),2.6],[1.5*TAN(PI/6),0,0],[0,1.5*TAN(PI/6),0]]);
tp2 = MAP(BEZIER(S2)([Su2,Su23,Su3]))(dom2D)

#center part
Su13 = CUBICHERMITE(S1)([[0.2+COS(PI/6),SIN(PI/12),1.9],[0.2+COS(PI/6),(-SIN(PI/12)),1.9],[-TAN(PI/6)/2,0,0],[0,TAN(PI/6)/2,0]]);
tp3 = MAP(BEZIER(S2)([Su3,Su13,Su1]))(dom2D)

#bottom
Su4 = CUBICHERMITE(S1)([[1*COS(PI/6),1*SIN(PI/12),-0.6],[1*COS(PI/6),-1*SIN(PI/12),-0.6],[1*TAN(PI/6),0,0],[0,1*TAN(PI/6),0]]);
Su40 = CUBICHERMITE(S1)([[1.3*COS(PI/6),1.5*SIN(PI/12),-0.5],[1.5*COS(PI/6),-1.5*SIN(PI/12),-0.5],[1.5*TAN(PI/6),0,0],[0,1.5*TAN(PI/6),0]]);
tp4 = MAP(BEZIER(S2)([Su0,Su40,Su4]))(dom2D)



partial_teapot_form = STRUCT([tp1,tp2,tp3,tp4])

partial_tp_black = COLOR([0,0,0])(partial_teapot_form)
partial_tp_white = COLOR([1,1,1])(R([1,2])(PI/6)(partial_teapot_form))

partial_tp = STRUCT([partial_tp_white,partial_tp_black])

teapot_form = STRUCT([partial_tp,R([1,2])(PI/3)]*6)

#knob
k0 = CUBICHERMITE(S2)([[-0.15,0,0],[0.15,0,0.5],[0,0,0.5],[0.5,0,0]])
k1 = CUBICHERMITE(S2)([[0,0,0],[0.15,0,0.3],[0,0,0.3],[0.3,0,0]])
k2 = CUBICHERMITE(S1)([[-0.15,0,0],[0,0,0],[0,-0.15,0],[0,0.15,0]])
k3 = CUBICHERMITE(S1)([[0.15,0,0.5],[0.15,0,0.3],[0,-0.15,0],[0,0.15,0]])
knob_right = MAP(COONSPATCH([k2,k3,k0,k1]))(dom2D);

k4 = CUBICHERMITE(S1)([[-0.15,0,0],[0,0,0],[0,0.15,0],[0,-0.15,0]])
k5 = CUBICHERMITE(S1)([[0.15,0,0.5],[0.15,0,0.3],[0,0.15,0],[0,-0.15,0]])
knob_left = MAP(COONSPATCH([k4,k5,k0,k1]))(dom2D);

knob = COLOR([0,0,0])(T([3])([3.6])(S([1,2,3])([0.35,0.35,0.35])(STRUCT([knob_left,knob_right]))))

domain_sphere = PROD([INTERVALS(PI/2)(24),INTERVALS(PI/6)(36)])

def partial_sphere (radius,domain):
	fx  = lambda p: radius * math.cos(p[0])  * math.sin  (p[1])
	fy  = lambda p: radius * math.cos(p[0]) * math.cos (p[1])
	fz  = lambda p: radius * math.sin(p[0]) 
	ret=  MAP([fx, fy, fz])(domain)
	return ret;


partial_cover = R([1,2])(PI/12)(partial_sphere(0.9,domain_sphere));
partial_cover_black = COLOR([0,0,0])(partial_cover)
partial_cover_white = COLOR([1,1,1])(R([1,2])(PI/6)(partial_cover))
partial_cover_bw = STRUCT([partial_cover_white,partial_cover_black])
partial_cv = STRUCT([partial_cover_bw,R([1,2])(PI/3)] *6)
cov = T([3])([3.1])(R([1,2])(PI/6)(S([3])([0.56])(partial_cv)))
cover = STRUCT([cov,knob])


#handle
#external_handle

Su0 = CUBICHERMITE(S1)([[-1,0.15,0.5],[-1,-0.15,0.5],[-0.6,0,0],[0.6,0,0]]);
Su1 = CUBICHERMITE(S1)([[-1.8,0.15,1.8],[-1.8,-0.15,1.8],[-0.6,0,0],[0.6,0,0]]);
Sv0 = BEZIER(S2)([[-1.8,0.15,1.8],[-1.3,0.15,1.5],[-1.5,0.15,0.3],[-1,0.15,0.5]]);
Sv1 = BEZIER(S2)([[-1.8,-0.15,1.8],[-1.3,-0.15,1.5],[-1.5,-0.15,0.3],[-1,-0.15,0.5]]);
ext_h1 = COLOR([0,0,0])(MAP(COONSPATCH([Su1,Su0,Sv0,Sv1]))(dom2D));
Su2 = CUBICHERMITE(S1)([[-1.8,0.15,3.2],[-1.8,-0.15,3.2],[0,0,0.45],[0,0,-0.45]]);
Sv2 = BEZIER(S2)([[-1.8,0.15,3.2],[-2.2,0.15,3.1],[-2.6,0.15,2.4],[-1.8,0.15,1.8]]);
Sv3 = BEZIER(S2)([[-1.8,-0.15,3.2],[-2.2,-0.15,3.1],[-2.6,-0.15,2.4],[-1.8,-0.15,1.8]]);
ext_h2 = COLOR([0,0,0])(MAP(COONSPATCH([Su2,Su1,Sv2,Sv3]))(dom2D));

Su3 = CUBIC_HERMITE(S1)([[-1.2,0.15,2.5],[-1.2,-0.15,2.5],[0.6,0,0],[-0.6,0,0]]);
Sv4 = BEZIER(S2)([[-1.2,0.15,2.5],[-0.9,0.15,3],[-1.5,0.15,3.3],[-1.8,0.15,3.2]]);
Sv5 = BEZIER(S2)([[-1.2,-0.15,2.5],[-0.9,-0.15,3],[-1.5,-0.15,3.3],[-1.8,-0.15,3.2]]);
ext_h3 = COLOR([0,0,0])(MAP(COONSPATCH([Su3,Su2,Sv4,Sv5]))(dom2D));

#internal_handle
Su0_int = CUBICHERMITE(S1)([[-1,0.15,0.5],[-1,-0.15,0.5],[0.6,0,0],[-0.6,0,0]]);
Su1_int = CUBICHERMITE(S1)([[-1.8,0.15,1.8],[-1.8,-0.15,1.8],[0.6,0,0],[-0.6,0,0]]);
int_h1 = COLOR([1,1,1])(MAP(COONSPATCH([Su1_int,Su0_int,Sv0,Sv1]))(dom2D));
Su2_int = CUBICHERMITE(S1)([[-1.8,0.15,3.2],[-1.8,-0.15,3.2],[0,0,-0.45],[0,0,0.45]]);
int_h2 = COLOR([1,1,1])(MAP(COONSPATCH([Su2_int,Su1_int,Sv2,Sv3]))(dom2D));
Su3_int = CUBICHERMITE(S1)([[-1.2,0.15,2.5],[-1.2,-0.15,2.5],[-0.6,0,0],[0.6,0,0]]);
int_h3b = MAP(COONSPATCH([Su3_int,Su2_int,Sv4,Sv5]))(dom2D);
int_h3 = COLOR([1,1,1])(int_h3b);
ext_h3b = COLOR([0,0,0])(T([3])([0.1])(int_h3b));

#int_h4 = MAP(BEZIER(S2)([Su3,Su3_int]))(dom2D)
handle = S([1,2,3])([0.35,0.35,0.35])(STRUCT([ext_h1,ext_h2,ext_h3b,int_h1,int_h2,int_h3]))

#bottom
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

#spout
sv1 = BEZIER(S2)([[1,0,0.2],[3,0,0.2],[1.5,0,2.8],[2.5,0,2.8]])
sv2 = BEZIER(S2)([[1,0,1.85],[2,0,1.75],[1.5,0,3.1],[2.3,0,3.1]])
su1 = BEZIER(S1)([[1,0,0.2],[1.2,0.35,0.7],[1.2,0.35,1.2],[1,0,1.85]])
su2 = BEZIER(S1)([[2.5,0,2.8],[2.4,0.2,2.9],[2.1,0.2,3],[2.3,0,3.1]])
spout1 = MAP(COONSPATCH([su1,su2,sv1,sv2]))(dom2D)
su3 = BEZIER(S1)([[1,0,0],[1.5,-0.35,0.7],[1.5,-0.35,1.2],[1,0,1.85]])
su4 = BEZIER(S1)([[2.5,0,2.8],[2.4,-0.2,2.9],[2.1,-0.2,3],[2.3,0,3.1]])
spout2 = MAP(COONSPATCH([su3,su4,sv1,sv2]))(dom2D)
spout = COLOR([1,1,1])(S([1,2,3])([0.35,0.35,0.35])(STRUCT([spout1,spout2])))

teapot = STRUCT([teapot_form,handle,spout,cover,teapot1,sb])
VIEW(teapot)