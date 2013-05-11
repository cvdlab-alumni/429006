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


def trasla (p,v):
    q = []
    length=len(p)
    for i in range(length):
        q += [ADD([p[i],v])]
    return q;


def bezier_s2 (h0, spessore):
    h00 = trasla(h0,spessore)
    ch0 = BEZIER(S1)(h0)
    ch00 = BEZIER(S1)(h00)
    beziers2 = BEZIER(S2)([ch00,ch0,])
    return MAP(beziers2)(dom2D);


def bezier_s2_2punti (h0,v1, spessore):
    h00 = trasla(h0,spessore)
    v11 = trasla(v1, spessore)
    ch0 = BEZIER(S1)(h0)
    ch00 = BEZIER(S1)(h00)
    cv1 = BEZIER(S1)(v1)
    cv11 = BEZIER(S1)(v11)
    beziers2 = BEZIER(S2)([cv1,ch0,ch00,cv11])
    return MAP(beziers2)(dom2D);


def reverse (p):
    q = []
    length=len(p)
    for i in range(length):
        q += [p[length-1-i]]
    return q;


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
tire2 = MAP(t01)(dom2D)
tire0 = STRUCT([tire1,tire2])
tire = COLOR(BLACK)(STRUCT([tire0, R([1,2])(PI/2)] * 4))

r0 = [[0,1,0],[1,0,0],[2,0,0],[0,-2,0]]
r00 = CUBICHERMITE(S1)(r0)
r1 = [[0,0.8,0],[0.8,0,0],[1.6,0,0],[0,-1.6,0]]
r11 = CUBICHERMITE(S1)(r1)
r10 = CUBICHERMITE(S2)([r11,r00,[0,0,1.3],[0,0,-1.3]])
rim1 = MAP(r10)(dom2D)
r01 = CUBICHERMITE(S2)([r00,r11,[0,0,-1.3],[0,0,1.3]])
rim2 = MAP(r01)(dom2D)
rim0 = STRUCT([rim1,rim2])
rim = COLOR(GRAY)(STRUCT([rim0, R([1,2])(PI/2)] * 4))

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

wheel1 = T([1,2,3])([-4.1,-3.8,-1.5])(wheel_r)
wheel2 = T([1,2,3])([6.15,-3.8,-1.5])(wheel_r)
wheel3 = R([1,2])(PI)(wheel1)
wheel4 = R([1,2])(PI)(wheel2)
wheel34 = T([1])([2.05])(STRUCT([wheel3,wheel4]))
wheels = STRUCT([wheel1,wheel2,wheel34])


#EXERCISE4: steering_wheel
p1 = [[0.6,0,0],[0,0.6,0],[-0.2,1.5,0],[-1.2,0,0]]
p2 = [[0.45,0,0],[0,0.45,0],[-0.125,1.125,0],[-0.9,0,0]]

p11 = [[0,-0.6,0],[0.6,0,0],[1.2,0,0],[-0.2,1.5,0],]
p22 = [[0,-0.45,0],[0.45,0,0],[0.9,0,0],[-0.125,1.125,0]]

c1 = CUBICHERMITE(S1)(p1)
c2 = CUBICHERMITE(S1)(p2)
c11 = CUBICHERMITE(S1)(p11)
c22 = CUBICHERMITE(S1)(p22)

s_w1 = CUBICHERMITE(S2)([c1,c2,[0,0,0.3],[0,0,-0.3]])
s_w2 = CUBICHERMITE(S2)([c2,c1,[0,0,-0.3],[0,0,0.3]])

s_wheel_up1 = MAP(s_w1)(dom2D)
s_wheel_down2 = MAP(s_w2)(dom2D)

s_w11 = CUBICHERMITE(S2)([c11,c22,[0,0,0.3],[0,0,-0.3]])
s_w22 = CUBICHERMITE(S2)([c22,c11,[0,0,-0.3],[0,0,0.3]])

s_wheel_up11 = MAP(s_w11)(dom2D)
s_wheel_down22 = MAP(s_w22)(dom2D)

s_wheelsurface1 = STRUCT([s_wheel_up1,s_wheel_down2,s_wheel_up11,s_wheel_down22])
s_wheelsurface2 = R([1,3])(PI)(s_wheelsurface1)
s_wheelsurface = COLOR([0.1,0.1,0.1])(STRUCT([s_wheelsurface1,s_wheelsurface2]))

center0 = bezier_surface([[0.09,0.45],[0,0.25],[0.3,0.45],[0.3,0]],
    [[-0.09,0.45],[0,0.25],[-0.3,0.45],[-0.3,0]])
center0_ex = PROD([center0,Q(0.02)])

center1 = bezier_surface([[-0.45,-0.1],[-0.2,-0.2],[0.2,-0.2],[0.45,-0.1]],
    [[-0.45,0],[0.45,0]])
center1_ex = PROD([center1,Q(0.02)])

d = COLOR([0.1,0.1,0.1])(T([2,3])([0.04,0.02])(CYLINDER([0.17,0.05])(36)))
d1 = COLOR(RED)(T([1,2,3])([-0.2,0.2,0.02])(CYLINDER([0.05,0.03])(36)))
d2 = COLOR(RED)(T([1,2,3])([0.17,0.24,0.02])(CYLINDER([0.04,0.03])(36)))
d3 = COLOR([0.1,0.1,0.1])(T([1,2,3])([0.21,0.14,0.02])(CYLINDER([0.03,0.03])(36)))

center = COLOR([0.5,0.5,0.5])(STRUCT([center0_ex,center1_ex]))
steering_wheel0 = STRUCT([s_wheelsurface,center,d,d1,d2,d3])
steering_wheel_s = S([1,2,3])([1.6,1.6,1.6])(steering_wheel0)
steering_wheel_r0 = R([1,2])(-PI/2)(steering_wheel_s)
steering_wheel_r1 = R([1,3])(-PI/3)(steering_wheel_r0)
steering_wheel = T([1,2])([-1.6,-1.3])(steering_wheel_r1)

VIEW(STRUCT([profile,wheels,steering_wheel]))