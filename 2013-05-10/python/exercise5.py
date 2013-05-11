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



#VIEW(STRUCT([profile,wheels,steering_wheel]))


#EXERCISE 5

p0 = [[-5.6,-3.8,-2],[-8.2,-3.8,-1.9],[-10.1,-3.8,-1.1],[-9,-3.8,-1.1],[-7.1,-3.8,0.3],[-5.6,-3.8,0.5]]
p1 = [[-5.6,-3.8,0.5],[-5.6,-3.8,-2]]
c0 = BEZIER(S1)(p0)
c1 = BEZIER(S1)(p1)

c01 = BEZIER(S2)([c1,c0])
c01map = MAP(c01)(dom2D)

######HERMITE##############
p2 = [[-5.6,-3.8,0.5],[-4.3,-3.8,0.6], [-3.9,-3.8,0.7],[-2.6,-3.8,0.5]]
p3 = [[-5.6,-3.8,-2],[-2.6,-3.8,-2],[0.3,0,8.9],[-0.3,0,-8.9]]
c2 = BEZIER(S1)(p2)
c3 = CUBICHERMITE(S1)(p3)

c23 = BEZIER(S2)([c3,c2])
c32 = BEZIER(S2)([c2,c3])
c23map = MAP(c23)(dom2D)
c32map = T(2)([7.6])(MAP(c32)(dom2D))

p4 = [[-2.6,-3.8,0.5],[-2.2,-3.8,0.35],[-1.8,-3.8,0],[3.4,-3.8,0],[4.6,-3.8,1.9]]
p5 = [[-2.6,-3.8,-2],[4.6,-3.8,-2]]
c4 = BEZIER(S1)(p4)
c5 = BEZIER(S1)(p5)

c45 = BEZIER(S2)([c5,c4])
c45map = MAP(c45)(dom2D)

p6 = [[4.6,-3.8,1.9],[7.2,-3.8,1],[7.6,-3.8,1]] 
p7 = [[4.6,-3.8,-2],[7.6,-3.8,-2],[0.3,0,8.9],[-0.3,0,-8.9]]
c6 = BEZIER(S1)(p6)
c7 = CUBICHERMITE(S1)(p7)

c67 = BEZIER(S2)([c7,c6])
c76 = BEZIER(S2)([c6,c7])
c67map = MAP(c67)(dom2D)
c76map = T(2)([7.6])(MAP(c76)(dom2D))

p8 = [[7.6,-3.8,1],[7.9,-3.8,1.1],[8.1,-3.8,0.9],[8.6,-3.8,1]]
p9 = [[7.6,-3.8,-2],[8.4,-3.8,-2],[8.6,-3.8,-1.6]]

c8 = BEZIER(S1)(p8)
c9 = BEZIER(S1)(p9)

c89 = BEZIER(S2)([c9,c8])
c89map = MAP(c89)(dom2D)

p10 = [[8.6,-3.8,1],[9.4,-3.8,-1.8],[8.9,-3.8,-1.3],[8.6,-3.8,-1.6]]
p11 = [[8.6,-3.8,1],[8.6,-3.8,-1.6]]
c10 = BEZIER(S1)(p10)
c11 = BEZIER(S1)(p11)

c1011 = BEZIER(S2)([c11,c10])
c1011map = MAP(c1011)(dom2D)

facciata_laterale_destra = STRUCT([c01map,c23map,c45map, c67map,c89map,c1011map])

#FACCIATA LATERALE SINISTRA
c00 = BEZIER(S1)(trasla(reverse(p0),[0,7.6,0]))
c11 = BEZIER(S1)(trasla(reverse(p1),[0,7.6,0]))
c44 = BEZIER(S1)(trasla(reverse(p4),[0,7.6,0]))
c55 = BEZIER(S1)(trasla(reverse(p5),[0,7.6,0]))
c88 = BEZIER(S1)(trasla(reverse(p8),[0,7.6,0]))
c99 = BEZIER(S1)(trasla(reverse(p9),[0,7.6,0]))
c1010 = BEZIER(S1)(trasla(reverse(p10),[0,7.6,0]))
c1111 = BEZIER(S1)(trasla(reverse(p11),[0,7.6,0]))

c0011 = BEZIER(S2)([c11,c00])
c0011map = MAP(c0011)(dom2D)
c4455 = BEZIER(S2)([c55,c44])
c4455map = MAP(c4455)(dom2D)
c8899 = BEZIER(S2)([c99,c88])
c8899map = MAP(c8899)(dom2D)
c10101111 = BEZIER(S2)([c1111,c1010])
c10101111map = MAP(c10101111)(dom2D)

facciata_laterale_sinistra = STRUCT([c32map,c0011map, c4455map,c8899map, c76map,c10101111map])
facciate_laterali = STRUCT([facciata_laterale_sinistra,facciata_laterale_destra])

#PARTE SUPERIORE
a0 = BEZIER(S1)([[-0.5,-2.8,2],[-0.5,-3,2.8],[1.2,-2.5,2.9],[2.5,-3,2.9],[4.6,-3.8,1.9]])
a1 = BEZIER(S1)([[-0.5,-2.8,2],[3.5,-3,0.8],[3,-3.8,-0.2], [4.6,-3.8,1.9]])

a00 = BEZIER(S1)([[-0.5,2.8,2],[-0.5,3,2.8],[1.2,2.5,2.9],[2.5,3,2.9],[4.6,3.8,1.9]])
a11 = BEZIER(S1)([[-0.5,2.8,2],[3.5,3,0.8],[3,3.8,-0.2], [4.6,3.8,1.9]])

a01 = BEZIER(S2)([a1,a0,a00,a11])
a01map = MAP(a01)(dom2D)


a2 = reverse([[-5.6,-3.8,-2],[-8.2,-3.8,-1.9],[-10.8,-3.8,-1.1],[-9,-3.8,-1.1],[-7.1,-3.8,0.3],
    [-5.6,-3.8,0.5],[-4.3,-3.8,0.6],[-3.5,-3.8,0.7],[-2.6,-3.8,0.5]])
a22 = bezier_s2(a2,[0,7.6,0])

a3 = reverse([[4.6,-3.8,1.9],[7.2,-3.8,1],[7.6,-3.8,1],[7.9,-3.8,1.1],[8.1,-3.8,0.9],
[8.6,-3.8,1]])
a33 = bezier_s2(a3,[0,7.6,0])
a4 = reverse([[8.6,-3.8,1],[9.4,-3.8,-1.8],[8.9,-3.8,-1.3],[8.6,-3.8,-1.6],[7.6,-3.8,-2]])
a44 = bezier_s2(a4,[0,7.6,0])

a5 = [[4.6,-3.8,-2],[-2.6,-3.8,-2]]
a55 = bezier_s2(a5,[0,7.6,0])

top = STRUCT([COLOR([0.1,0.1,0.1])(a01map),a22,a33,a44,a55])
carrozzeria = COLOR(YELLOW)(STRUCT([facciate_laterali,top]))

#FINESTRINI
b0 = BEZIER(S1)([[-2.6,-3.6,0.6],[-0.6,-2.6,2]])
b1 = BEZIER(S1)([[-2.6,3.6,0.6],[-0.6,2.6,2]])
b2 = BEZIER(S2)([b0,b1])

glass = COLOR([230.0/255, 1.0, 1.0, 0.5])(MAP(b2)(dom2D))

car = STRUCT([glass,carrozzeria,profile,wheels,steering_wheel])
VIEW(car)