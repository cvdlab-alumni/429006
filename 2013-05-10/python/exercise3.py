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
#funzioni

domain1D = INTERVALS(1)(32)
dom2D = GRID([20,20])

def image_bezierline(points):
	c = BEZIER(S1)(points)
	return MAP(c)(domain1D)


def image_hermiteline(points):
	c = CUBICHERMITE(S1)(points)
	return MAP(c)(domain1D)


def bezier_surface(points1,points2):
	c1 = BEZIER(S1)(points1)
	c2 = BEZIER(S1)(points2)
	c = BEZIER(S2)([c2,c1])
	return MAP(c)(dom2D)



#PROFILO SOPRA
profile_s0 = image_bezierline([[-9,0,0],[-8.8,2.3,0],[-8,3.3,0],[-5,3.8,0],[6.8,3.8,0],[8.4,3.8,0],[8.9,3.4,0],[9,0,0]])
profile_s00 = image_bezierline([[-9,0,0],[-8.8,-2.3,0],[-8,-3.3,0],[-5,-3.8,0],[6.8,-3.8,0],[8.4,-3.8,0],[8.9,-3.4,0],[9,0,0]])

profile_s = STRUCT([profile_s0,profile_s00])

#PROFILO DAVANTI
profile_d0 = image_bezierline([[0,0,2.3],[0,1.9,2.3],[0,2.5,1.6],[0,2.7,1.1]])
profile_d00 = image_bezierline([[0,0,2.3],[0,-1.9,2.3],[0,-2.5,1.6],[0,-2.7,1.1]])

profile_d1 = image_bezierline([[0,2.7,1.1],[0,3.3,1.1],[0,3.7,0.6],[0,3.7,-2]])
profile_d11 = image_bezierline([[0,-2.7,1.1],[0,-3.3,1.1],[0,-3.7,0.6],[0,-3.7,-2]])

profile_d2 = image_bezierline([[0,3.7,-2],[0,0,-2]])
profile_d22 = image_bezierline([[0,-3.7,-2],[0,0,-2]])

profile_d = STRUCT([profile_d0,profile_d1, profile_d2,profile_d00, profile_d11, profile_d22])

#PROFILO LATERALE
profile_l0 = image_bezierline([[-9,0,-1.1],[-7.1,0,0.5],[-4.3,0,0.6],[-3.5,0,0.7]])
profile_l1 = image_bezierline([[-3.5,0,0.7],[-1.4,0,2.8],[3.4,0,2.8],[7.2,0,1]])
profile_l2 = image_bezierline([[7.2,0,1],[7.9,0,1.1],[8.1,0,0.9],[8.6,0,1]])
profile_l3 = image_bezierline([[8.6,0,1],[9.4,0,-1.8],[8.9,0,-1.3],[8.6,0,-1.6]])
profile_l4 = image_bezierline([[8.6,0,-1.6],[8.4,0,-2],[7.6,0,-2]])
profile_l5 = image_hermiteline([[7.6,0,-2],[4.6,0,-2],[0.3,0,8.9],[-0.3,0,-8.9]])
profile_l6 = image_bezierline([[4.6,0,-2],[-2.6,0,-2]])
profile_l7 = image_hermiteline([[-2.6,0,-2],[-5.6,0,-2],[0.3,0,8.9],[-0.3,0,-8.9]])
profile_l8 = image_bezierline([[-5.6,0,-2],[-8.2,0,-1.9],[-8.1,0,-2.1],[-9,0,-1.1]])

profile_l = STRUCT([profile_l0,profile_l1,profile_l2,profile_l3,profile_l4,profile_l5,profile_l6,profile_l7,profile_l8])

profile = STRUCT([profile_d, profile_s,profile_l])
#EXERCISE 3: WHEEL
t0 = [[0,1.5,0],[1.5,0,0],[3,0,0],[0,-3,0]]
t00 = CUBICHERMITE(S1)(t0)
t1 = [[0,1,0],[1,0,0],[2,0,0],[0,-2,0]]
t11 = CUBICHERMITE(S1)(t1)
t10 = CUBICHERMITE(S2)([t11,t00,[0,0,2],[0,0,-2]])
tire1 = MAP(t10)(dom2D)
t01 = CUBICHERMITE(S2)([t00,t11,[0,0,-2],[0,0,2]])
tire2 = MAP(w01)(dom2D)
tire0 = STRUCT([tire1,tire2])
tire = COLOR(BLACK)(STRUCT([tire0, R([1,2])(PI/2)] * 4))

r0 = [[0,1,0],[1,0,0],[2,0,0],[0,-2,0]]
r00 = CUBICHERMITE(S1)(r0)
r1 = [[0,0.8,0],[0.8,0,0],[1.6,0,0],[0,-1.6,0]]
r11 = CUBICHERMITE(S1)(r1)
r10 = CUBICHERMITE(S2)([r11,r00,[0,0,1],[0,0,-1]])
rim1 = MAP(r10)(dom2D)
r01 = CUBICHERMITE(S2)([r00,r11,[0,0,-1],[0,0,1]])
rim2 = MAP(r01)(dom2D)
rim0 = STRUCT([rim1,rim2])
rim = COLOR(GRAY)(STRUCT([rim0, R([1,2])(PI/2)] * 4))


r0 = [[0,1.1,0.3],[1.1,0,0.3],[2,0,0],[0,-2,0]]
r00 = CUBICHERMITE(S1)(r0)
r1 = [[0,0.8,0.1],[0.8,0,0.1],[1.8,0,0],[0,-1.8,0]]
r11 = CUBICHERMITE(S1)(r1)
r = BEZIER(S2)([r11,r00])
rm = MAP(r)(dom2D)
rm_r = COLOR(GRAY)(STRUCT([rm, R([1,2])(PI/2)] * 4))
VIEW(STRUCT([rm_r,tire]))

p = POLYLINE([[0,0],[0.2,0.2],[0,0.85],[-0.2,0.2],[0,0]])
s = SOLIDIFY(p)
s1 = PROD([s,Q(0.1)])
t = POLYLINE([[0,0],[-0.19,0.32],[0.19,0.32],[0,0]])
st = SOLIDIFY(t)
st1 = R([2,3])(-PI/5)(st)

rm_s = COLOR(GRAY)(STRUCT([s1, R([1,2])(2*PI/5)] * 5))
rm_st = T([3])([0.3])(COLOR(GRAY)(STRUCT([st1, R([1,2])(2*PI/5)] * 5)))

disk = T(3)(-0.1)(COLOR([0.2,0.2,0.2])(CYLINDER([1,0.01])(36)))

wheel = S([1,2,3])([0.9,0.9,0.9])(STRUCT([disk,rim,tire,rm_s,rm_st]))

wheel_r = R([2,3])(PI/2)(wheel)

wheel1 = T([1,2,3])([-4.1,-2.5,-1.5])(wheel_r)
wheel2 = T([1,2,3])([6.15,-2.5,-1.5])(wheel_r)
wheel3 = R([1,2])(PI)(wheel1)
wheel4 = R([1,2])(PI)(wheel2)
wheel34 = T([1])([2.05])(STRUCT([wheel3,wheel4]))
wheels = STRUCT([wheel1,wheel2,wheel34])

VIEW(STRUCT([profile,wheels]))