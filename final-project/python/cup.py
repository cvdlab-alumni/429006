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



#saucer_form
su0 = CUBICHERMITE(S1)([[COS(PI/6),SIN(PI/12),0],[COS(PI/6),-SIN(PI/12),0],[TAN(PI/6),0,0],[0,TAN(PI/6),0]]);
su1 = CUBICHERMITE(S1)([[0.5+2.3*COS(PI/6),0.1+2.3*SIN(PI/12),1],[0.5+2.3*COS(PI/6),-0.1 +2.3*(-SIN(PI/12)),1],[2.3*TAN(PI/6),0,0],[0,2.3*TAN(PI/6),0]]);

saucer_f = MAP(BEZIER(S2)([su0,su1]))(dom2D)

su2 = CUBICHERMITE(S1)([[0.5+2.3*COS(PI/6),0.1+2.3*SIN(PI/12),1],[0.5+2.3*COS(PI/6),-0.1 +2.3*(-SIN(PI/12)),1],[2.3*TAN(PI/6),0,2],[0,2.3*TAN(PI/6),-2]]);
su12 = MAP(BEZIER(S2)([su1,su2]))(dom2D)

saucer_f_black = COLOR([0,0,0])(R([1,2])(PI/6)(STRUCT([saucer_f,su12])))
saucer_f_white = COLOR([1,1,1])(STRUCT([saucer_f,su12]))
saucer_form = STRUCT([saucer_f_black,saucer_f_white])
saucer_top = STRUCT([saucer_form, R([1,2])(PI/3)] * 6)

def arc(alpha,r,R):
	domain = PROD([INTERVALS(alpha)(36),INTERVALS(R-r)(1)])
	def mapping(v):
		a = v[0];
		r = v[1];
		return [r*COS(a), r*SIN(a)];
	model = MAP(mapping)(domain);
	return model;



arc1 = arc(PI/6,0,0.98);
arc_black = COLOR([0,0,0])(R([1,2])(PI/12)(arc1))
arc_white = COLOR([1,1,1])(R([1,2])(-PI/12)(arc1))
partial_saucer_bottom = STRUCT([arc_black,arc_white])
saucer_bottom = STRUCT([partial_saucer_bottom,R([1,2])(PI/3)]*6)
saucer = STRUCT([saucer_bottom,saucer_top])

#cup
s0 = CUBICHERMITE(S1)([[COS(PI/6),SIN(PI/6),0],[COS(PI/6),-SIN(PI/6),0],[COS(PI/3),-SIN(PI/3),0],[-COS(PI/3),-SIN(PI/3),0]]);
s1 = CUBICHERMITE(S1)([[2*COS(PI/6),2*SIN(PI/6),2],[2*COS(PI/6),2*(-SIN(PI/6)),2],[2*COS(PI/3),-2*SIN(PI/3),0],[-2*COS(PI/3),-2*SIN(PI/3),0]]);
cup_f_white = COLOR([1,1,1])(MAP(BEZIER(S2)([s1,s0]))(dom2D))
cup_f_black = COLOR([0,0,0])(R([1,2])(PI/3)(MAP(BEZIER(S2)([s1,s0]))(dom2D)))
partial_cup_form = STRUCT([cup_f_white,cup_f_black])
cup_form = T([3])([0.2])(STRUCT([partial_cup_form, R([1,2])(2*PI/3)] * 3))

arc2 = arc(PI/3,0,1);
arc2_black = COLOR([0,0,0])(R([1,2])(PI/6)(arc2))
arc2_white = COLOR([1,1,1])(R([1,2])(-PI/6)(arc2))
partial_cup_bottom = STRUCT([arc2_white,arc2_black])
cup_bottom = T([3])([0.2])(STRUCT([partial_cup_bottom,R([1,2])(2*PI/3)]*3))

#handle
Su0 = CUBICHERMITE(S1)([[-1,0.15,0.5],[-1,-0.15,0.5],[-0.6,0,0],[0.6,0,0]]);
Su1 = CUBICHERMITE(S1)([[-1.9,0.15,1],[-1.9,-0.15,1],[-0.6,0,0],[0.6,0,0]]);
Sv0 = BEZIER(S2)([[-1.9,0.15,1],[-1.3,0.15,0.8],[-1.7,0.15,0.4],[-1,0.15,0.5]]);
Sv1 = BEZIER(S2)([[-1.9,-0.15,1],[-1.3,-0.15,0.8],[-1.7,-0.15,0.4],[-1,-0.15,0.5]]);
ext_h1 = COLOR([0,0,0])(MAP(COONSPATCH([Su1,Su0,Sv0,Sv1]))(dom2D));
Su2 = CUBICHERMITE(S1)([[-2.5,0.15,2.05],[-2.5,-0.15,2.05],[0,0,0.45],[0,0,-0.45]]);
Sv2 = BEZIER(S2)([[-2.5,0.15,2.05],[-2.9,0.15,1.9],[-3.3,0.15,1.2],[-2.5,0.15,0.8]]);
Sv3 = BEZIER(S2)([[-2.5,-0.15,2.05],[-2.9,-0.15,1.9],[-3.3,-0.15,1.2],[-2.5,-0.15,0.8]]);
ext_h2 = COLOR([0,0,0])(MAP(COONSPATCH([Su2,Su1,Sv2,Sv3]))(dom2D));
Su3 = CUBICHERMITE(S1)([[-2,0.15,1.5],[-2,-0.15,1.5],[0.6,0,0],[-0.6,0,0]]);
Sv4 = BEZIER(S2)([[-2,0.15,1.5],[-1.7,0.15,2],[-2.2,0.15,2.3],[-2.5,0.15,2.05]]);
Sv5 = BEZIER(S2)([[-2,-0.15,1.5],[-1.7,-0.15,2],[-2.2,-0.15,2.3],[-2.5,-0.15,2.05]]);
ext_h3 = COLOR([0,0,0])(MAP(COONSPATCH([Su3,Su2,Sv4,Sv5]))(dom2D));

#internal_handle
Su0_int = CUBICHERMITE(S1)([[-1,0.15,0.5],[-1,-0.15,0.5],[0.6,0,0],[-0.6,0,0]]);
Su1_int = CUBICHERMITE(S1)([[-1.9,0.15,1],[-1.9,-0.15,1],[0.6,0,0],[-0.6,0,0]]);
int_h1 = COLOR([1,1,1])(MAP(COONSPATCH([Su1_int,Su0_int,Sv0,Sv1]))(dom2D));
Su2_int = CUBICHERMITE(S1)([[-2.5,0.15,2.05],[-2.5,-0.15,2.05],[0,0,-0.45],[0,0,0.45]]);
int_h2 = COLOR([1,1,1])(MAP(COONSPATCH([Su2_int,Su1_int,Sv2,Sv3]))(dom2D));
Su3_int = CUBICHERMITE(S1)([[-2,0.15,1.5],[-2,-0.15,1.5],[-0.6,0,0],[0.6,0,0]]);
int_h3 = COLOR([1,1,1])(MAP(COONSPATCH([Su3_int,Su2_int,Sv4,Sv5]))(dom2D));


upper = S([1,2,3])([0.33,0.33,0.33])(T([1,3])([1.25,0.38])(STRUCT([ext_h2,ext_h3,int_h2,int_h3])))
under = S([1,2,3])([0.33,0.33,0.33])(STRUCT([ext_h1,int_h1]))
handle = T([1,3])([-0.2,-0.1])(STRUCT([upper,under]))

cup = STRUCT([saucer,cup_bottom,cup_form,handle])
VIEW(cup)