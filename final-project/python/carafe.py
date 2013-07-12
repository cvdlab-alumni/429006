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

def move(p,h,s):
	p1 = p;
	p1[0] = [p[0][0]+s[0],p[0][1]+s[1],p[0][2]+h];
	p1[1] = [p[1][0]+s[0],p[1][1]+s[1],p[1][2]+h];
	return p1;


def side_form(p1,p2,s):
	c0 = CUBICHERMITE(S1)(p1)
	c1 = CUBICHERMITE(S1)(move(p1,0.5,[0,s*0.5,0]))
	c2 = CUBICHERMITE(S1)(move(p1,1.8,[0,s*0.3,0]))
	c3 = BEZIER(S1)(p2)
	b2 = MAP(BEZIER(S2)([c0,c1,c2,c3]))(dom2D)
	return b2;


#carafe's form

#body
#bezier point
b0 = [[1.5,0,2.9],[COS(PI/10)+0.1,-SIN(PI/10)+0.08,2.95],[COS(PI/10),-SIN(PI/10)+0.2,3]]
b1 = [[COS(PI/10),-SIN(PI/10),3],[COS(PI/10)-0.05,-SIN(PI/10)-0.15,3.04],[COS(3*PI/10)-0.05,-SIN(3*PI/10)-0.15,3.18],[COS(3*PI/10),-SIN(3*PI/10),3.1]]
b2 = [[COS(3*PI/10),-SIN(3*PI/10),3.1],[COS(3*PI/10)-0.05,-SIN(3*PI/10)-0.05,3.2],[COS(PI/2)-0.05,-SIN(PI/2)-0.05,3.05],[COS(PI/2),-SIN(PI/2),3]]
b3 = [[COS(PI/2),-SIN(PI/2),3],[COS(PI/2)-0.05,-SIN(PI/2)-0.05,3.05],[COS(7*PI/10)-0.05,-SIN(7*PI/10)-0.05,2.9],[COS(7*PI/10),-SIN(7*PI/10),2.82]]
b4 = [[COS(7*PI/10),-SIN(7*PI/10),2.82],[COS(7*PI/10)-0.3,-SIN(7*PI/10)+0.2,2.8],[COS(9*PI/10)-0.05,-SIN(9*PI/10)+0.2,2.7],[COS(9*PI/10),-SIN(9*PI/10),2.75]]
b5 = [[COS(9*PI/10),-SIN(9*PI/10),2.75],[COS(9*PI/10)-0.1,-SIN(9*PI/10)+0.2,2.8],[COS(PI),-SIN(PI),2.8]]
b6 = [[1,0,0],[1.9,0,1],[1,0,1.8],[1.5,0,2.9]]

#hermite point
p0 = [[1,0,0],[COS(PI/10),-SIN(PI/10),0],[COS(PI/10)/2,0,0],[-SIN(PI/10)/3,-2*COS(PI/10)/3,0]]
p01 = [[COS(PI/10),-SIN(PI/10),0],[COS(PI/10),-SIN(PI/10),2.95],[0,-1.95,0],[0,1.95,0]]
p1 = [[COS(PI/10),-SIN(PI/10),0],[COS(3*PI/10),-SIN(3*PI/10),0],[-SIN(PI/10)/3,-2*COS(PI/10)/3,0],[-2*SIN(PI/5)/3,-COS(PI/5)/3,0]]
p2 = [[COS(3*PI/10),-SIN(3*PI/10),0],[COS(PI/2),-SIN(PI/2),0],[-2*SIN(PI/5)/3,-COS(PI/5)/3,0],[-SIN(PI/5),0,0]]
p3 = [[COS(PI/2),-SIN(PI/2),0],[COS(7*PI/10),-SIN(7*PI/10),0],[-SIN(PI/5),0,0],[-2*SIN(PI/5)/3,COS(PI/5)/3,0]]
p4 = [[COS(7*PI/10),-SIN(7*PI/10),0],[COS(9*PI/10),-SIN(9*PI/10),0],[-2*SIN(PI/5)/3,COS(PI/5)/3,0],[-SIN(PI/10)/3,2*COS(PI/10)/3,0]]
p5 = [[COS(9*PI/10),-SIN(9*PI/10),0],[COS(PI),-SIN(PI),0],[-SIN(PI/10)/2,2*COS(PI/10)/3,0],[0,SIN(PI/10),0]]
p6 = [[-1,0,0],[-1,0,2.9],[-0.5,0,0],[0.5,0,0]]
p7 = [[COS(9*PI/10),-SIN(9*PI/10),0],[COS(9*PI/10),-SIN(9*PI/10),2.75],[0,-2,0],[0,2,0]]


#coons_patch
#front part
Su0 = CUBICHERMITE(S1)(p0)
Su1 = BEZIER(S1)(b0)
Sv0 = CUBICHERMITE(S2)(p01)
Sv1 = BEZIER(S2)(b6)
spout0 = S([1,2,3])([0.36,0.44,0.34])(MAP(COONSPATCH([Su0,Su1,Sv1,Sv0]))(dom2D))

#rotate point
p0r = [[1,0,0],[COS(PI/10),SIN(PI/10),0],[COS(PI/10)/2,0,0],[-SIN(PI/10)/3,2*COS(PI/10)/3,0]]
b0r = [[1.5,0,2.9],[COS(PI/10)+0.1,SIN(PI/10)+0.08,2.95],[COS(PI/10),SIN(PI/10)-0.2,3]]
p01r= [[COS(PI/10),SIN(PI/10),0],[COS(PI/10),SIN(PI/10),2.95],[0,1.95,0],[0,-1.95,0]]

Su00 = CUBICHERMITE(S1)(p0r)
Su11 = BEZIER(S1)(b0r)
Sv00 = CUBICHERMITE(S2)(p01r)
Sv11 = BEZIER(S2)(b6)

spout1 = S([1,2,3])([0.36,0.44,0.34])(MAP(COONSPATCH([Su00,Su11,Sv11,Sv00]))(dom2D))
spout = COLOR([1,1,1])(T([1])([-0.1])(STRUCT([spout0,spout1])))

#behind part
Su2 = CUBICHERMITE(S1)(p5)
Su3 = BEZIER(S1)(b5)
Sv22 = CUBICHERMITE(S2)(p6)
Sv33 = CUBICHERMITE(S2)(p7)
behindp0 = S([1,2,3])([0.34,0.42,0.33])(MAP(COONSPATCH([Su2,Su3,Sv33,Sv22]))(dom2D))

#rotate point
p5r = [[COS(9*PI/10),SIN(9*PI/10),0],[COS(PI),SIN(PI),0],[-SIN(PI/10)/2,-2*COS(PI/10)/3,0],[0,-SIN(PI/10),0]]
b5r = [[COS(9*PI/10),SIN(9*PI/10),2.75],[COS(9*PI/10)-0.1,SIN(9*PI/10)+0.2,2.8],[COS(PI),SIN(PI),2.8]]
p7r = [[COS(9*PI/10),SIN(9*PI/10),0],[COS(9*PI/10),SIN(9*PI/10),2.75],[0,2,0],[0,-2,0]]

Su2 = CUBIC_HERMITE(S1)(p5r)
Su3 = BEZIER(S1)(b5r)
Sv22 = CUBICHERMITE(S2)(p6)
Sv33 = CUBICHERMITE(S2)(p7r)
behindp1 = S([1,2,3])([0.34,0.48,0.33])(MAP(COONSPATCH([Su2,Su3,Sv33,Sv22]))(dom2D))
behind_part = COLOR([0,0,0])(STRUCT([behindp0,behindp1]))

#side

form1 = COLOR([0,0,0])(side_form(p1,b1,-1))
form2 = COLOR([1,1,1])(side_form(p2,b2,-1))
form3 = COLOR([0,0,0])(side_form(p3,b3,-1))
form4 = COLOR([1,1,1])(side_form(p4,b4,-1))

#rotate point
p1r = [[COS(PI/10),SIN(PI/10),0],[COS(3*PI/10),SIN(3*PI/10),0],[-SIN(PI/10)/3,2*COS(PI/10)/3,0],[-2*SIN(PI/5)/3,COS(PI/5)/3,0]]
p2r = [[COS(3*PI/10),SIN(3*PI/10),0],[COS(PI/2),SIN(PI/2),0],[-2*SIN(PI/5)/3,COS(PI/5)/3,0],[-SIN(PI/5),0,0]]
p3r = [[COS(PI/2),SIN(PI/2),0],[COS(7*PI/10),SIN(7*PI/10),0],[-SIN(PI/5),0,0],[-2*SIN(PI/5)/3,-COS(PI/5)/3,0]]
p4r = [[COS(7*PI/10),SIN(7*PI/10),0],[COS(9*PI/10),SIN(9*PI/10),0],[-2*SIN(PI/5)/3,-COS(PI/5)/3,0],[-SIN(PI/10)/3,-2*COS(PI/10)/3,0]]

b1r = [[COS(PI/10),SIN(PI/10),3],[COS(PI/10)-0.05,SIN(PI/10)-0.15,3.04],[COS(3*PI/10)-0.05,SIN(3*PI/10)-0.15,3.18],[COS(3*PI/10),SIN(3*PI/10),3.1]]
b2r = [[COS(3*PI/10),SIN(3*PI/10),3.1],[COS(3*PI/10)-0.05,SIN(3*PI/10)-0.05,3.2],[COS(PI/2)-0.05,SIN(PI/2)-0.05,3.05],[COS(PI/2),SIN(PI/2),3]]
b3r = [[COS(PI/2),SIN(PI/2),3],[COS(PI/2)-0.05,SIN(PI/2)-0.05,3.05],[COS(7*PI/10)-0.05,SIN(7*PI/10)-0.05,2.9],[COS(7*PI/10),SIN(7*PI/10),2.82]]
b4r = [[COS(7*PI/10),SIN(7*PI/10),2.82],[COS(7*PI/10)-0.3,SIN(7*PI/10)+0.2,2.8],[COS(9*PI/10)-0.05,SIN(9*PI/10)+0.2,2.7],[COS(9*PI/10),SIN(9*PI/10),2.75]]


form11 = COLOR([0,0,0])(side_form(p1r,b1r,1))
form22 = COLOR([1,1,1])(side_form(p2r,b2r,1))
form33 = COLOR([0,0,0])(side_form(p3r,b3r,1))
form44 = COLOR([1,1,1])(side_form(p4r,b4r,1))

carafe_form = STRUCT([spout,form1,form11,form2,form22,form3,form33,form4,form44,behind_part])

#handle
hSu0 = CUBICHERMITE(S1)([[-1.2,0.18,0.5],[-1.2,-0.15,0.5],[0,0.6,0],[0,-0.6,0]]);
hSu1 = CUBICHERMITE(S1)([[-1.8,0.15,1.8],[-1.8,-0.15,1.8],[-0.6,0,0],[0.6,0,0]]);
hSv0 = BEZIER(S2)([[-1.8,0.15,1.8],[-1.3,0.15,1.5],[-1.5,0.15,0.3],[-1.2,0.15,0.5]]);
hSv1 = BEZIER(S2)([[-1.8,-0.15,1.8],[-1.3,-0.15,1.5],[-1.5,-0.15,0.3],[-1.2,-0.15,0.5]]);
hext_h1 = COLOR([0,0,0])(MAP(COONSPATCH([hSu1,hSu0,hSv0,hSv1]))(dom2D));
hSu2 = CUBICHERMITE(S1)([[-1.8,0.15,3.2],[-1.8,-0.15,3.2],[0,0,0.45],[0,0,-0.45]]);
hSv2 = BEZIER(S2)([[-1.8,0.15,3.2],[-2.2,0.15,3.1],[-2.6,0.15,2.4],[-1.8,0.15,1.8]]);
hSv3 = BEZIER(S2)([[-1.8,-0.15,3.2],[-2.2,-0.15,3.1],[-2.6,-0.15,2.4],[-1.8,-0.15,1.8]]);
hext_h2 = COLOR([0,0,0])(MAP(COONSPATCH([hSu2,hSu1,hSv2,hSv3]))(dom2D));
hSu3 = CUBICHERMITE(S1)([[-1.23,0.15,2.5],[-1.23,-0.15,2.5],[0.6,0,0],[-0.6,0,0]]);
hSv4 = BEZIER(S2)([[-1.23,0.15,2.5],[-0.9,0.15,3],[-1.5,0.15,3.3],[-1.8,0.15,3.2]]);
hSv5 = BEZIER(S2)([[-1.23,-0.15,2.5],[-0.9,-0.15,3],[-1.5,-0.15,3.3],[-1.8,-0.15,3.2]]);
hext_h3 = COLOR([0,0,0])(MAP(COONSPATCH([hSu3,hSu2,hSv4,hSv5]))(dom2D));

#internal_handle
hSu0_int = CUBICHERMITE(S1)([[-1.2,0.15,0.5],[-1.2,-0.15,0.5],[0,0.6,0],[0,-0.6,0]]);
hSu1_int = CUBICHERMITE(S1)([[-1.8,0.15,1.8],[-1.8,-0.15,1.8],[0.6,0,0],[-0.6,0,0]]);
hint_h1 = COLOR([1,1,1])(MAP(COONSPATCH([hSu1_int,hSu0_int,hSv0,hSv1]))(dom2D));
hSu2_int = CUBICHERMITE(S1)([[-1.8,0.15,3.2],[-1.8,-0.15,3.2],[0,0,-0.45],[0,0,0.45]]);
hint_h2 = COLOR([1,1,1])(MAP(COONSPATCH([hSu2_int,hSu1_int,hSv2,hSv3]))(dom2D));
hSu3_int = CUBICHERMITE(S1)([[-1.23,0.15,2.5],[-1.23,-0.15,2.5],[-0.6,0,0],[0.6,0,0]]);
hint_h3 = COLOR([1,1,1])(MAP(COONSPATCH([hSu3_int,hSu2_int,hSv4,hSv5]))(dom2D));
hint_h4 = MAP(BEZIER(S2)([hSu3,hSu3_int]))(dom2D)
handle = T([1])([0.17])(S([1,2,3])([0.35,0.35,0.35])(STRUCT([hext_h1,hext_h2,hext_h3,hint_h1,hint_h2,hint_h3,hint_h4])))

Su0 = CUBICHERMITE(S1)([[COS(PI/10),-SIN(PI/10),0],[COS(3*PI/10),-SIN(3*PI/10),0],[-SIN(PI/10)/3,-2*COS(PI/10)/3,0],[-2*SIN(PI/5)/3,-COS(PI/5)/3,0]])
Su4 = CUBICHERMITE(S1)([[1.1*COS(PI/10),-1.1*SIN(PI/10),-0.5],[1.1*COS(3*PI/10),-1.1*SIN(3*PI/10),-0.5],[-2*SIN(PI/10)/6,-4*COS(PI/10)/6,0],[-4*SIN(PI/5)/6,-2*COS(PI/5)/6,0]])
tp4 = MAP(BEZIER(S2)([Su0,Su4]))(dom2D)

def arc(alpha,r,R):
    domain = PROD([INTERVALS(alpha)(36),INTERVALS(R-r)(1)])
    def mapping(v):
        a = v[0];
        r = v[1];
        return [r*COS(a), r*SIN(a)];
    model = MAP(mapping)(domain);
    return model;


arc1 = arc(PI/5,0,1.1);
arc_black = COLOR([0,0,0])(R([1,2])(PI/10)(arc1))
arc_white = COLOR([1,1,1])(R([1,2])(-PI/10)(arc1))
partial_sb_bottom = STRUCT([arc_black,arc_white])
sb_bottom = STRUCT([partial_sb_bottom,R([1,2])(2*PI/5)]*6)
sb = T([3])([-0.5])(sb_bottom)


partial_tp_black = COLOR([0,0,0])(tp4)
partial_tp_white = COLOR([1,1,1])(R([1,2])(PI/5)(tp4))
partial_tp = STRUCT([partial_tp_white,partial_tp_black])
carafebottom = STRUCT([partial_tp,R([1,2])(2*PI/5)]*5)

carafe = STRUCT([carafe_form,handle,carafebottom,sb])