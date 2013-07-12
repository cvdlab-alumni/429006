//cornicioni
var domain = PROD1x1([INTERVALS(1)(14),INTERVALS(1)(14)]);
var dom = INTERVALS(1)(14)
var c1 = CUBIC_HERMITE(S0)([[5,5.5,5],[4.85,8,2.5],[0,4,0],[0,0,-4]]);
var c2 = CUBIC_HERMITE(S0)([[4.85,8,2.5],[4.7,5.5,0],[0,0,-4],[0,-4,0]]);
var c3 = CUBIC_HERMITE(S0)([[5,5.5,4.8],[4.85,7.8,2.5],[0,3.4,0],[0,0,-3.4]]);
var c4 = CUBIC_HERMITE(S0)([[4.85,7.8,2.5],[4.7,5.5,0.2],[0,0,-3.4],[0,-3.4,0]]);

var c11 = CUBIC_HERMITE(S0)([[4.5,5.5,5],[4.5,8,2.5],[0,4,0],[0,0,-4]]);
var c22 = CUBIC_HERMITE(S0)([[4.5,8,2.5],[4.5,5.5,0],[0,0,-4],[0,-4,0]]);
var c33 = CUBIC_HERMITE(S0)([[4.5,5.5,4.8],[4.5,7.8,2.5],[0,3.4,0],[0,0,-3.4]]);
var c44 = CUBIC_HERMITE(S0)([[4.5,7.8,2.5],[4.5,5.5,0.2],[0,0,-3.4],[0,-3.4,0]]);

var out1 = MAP(BEZIER(S1)([c1,c11]))(domain);
var out2 = MAP(BEZIER(S1)([c2,c22]))(domain);
var out3 = MAP(BEZIER(S1)([c3,c33]))(domain);
var out4 = MAP(BEZIER(S1)([c4,c44]))(domain);
var out5 = MAP(BEZIER(S1)([c1,c3]))(domain);
var out6 = MAP(BEZIER(S1)([c2,c4]))(domain);
var out7 = MAP(BEZIER(S1)([c11,c33]))(domain);
var out8 = MAP(BEZIER(S1)([c22,c44]))(domain);

s = STRUCT ([out1,out2,out3,out4,out5,out6,out7,out8])

var c1b = CUBIC_HERMITE(S0)([[-0.5,5.5,5],[-0.35,8,2.5],[0,4,0],[0,0,-4]]);
var c2b = CUBIC_HERMITE(S0)([[-0.35,8,2.5],[-0.2,5.5,0],[0,0,-4],[0,-4,0]]);
var c3b = CUBIC_HERMITE(S0)([[-0.5,5.5,4.8],[-0.35,7.8,2.5],[0,3.4,0],[0,0,-3.4]]);
var c4b = CUBIC_HERMITE(S0)([[-0.35,7.8,2.5],[-0.2,5.5,0.2],[0,0,-3.4],[0,-3.4,0]]);

var c11b = CUBIC_HERMITE(S0)([[0,5.5,5],[0,8,2.5],[0,4,0],[0,0,-4]]);
var c22b= CUBIC_HERMITE(S0)([[0,8,2.5],[0,5.5,0],[0,0,-4],[0,-4,0]]);
var c33b = CUBIC_HERMITE(S0)([[0,5.5,4.8],[0,7.8,2.5],[0,3.4,0],[0,0,-3.4]]);
var c44b = CUBIC_HERMITE(S0)([[0,7.8,2.5],[0,5.5,0.2],[0,0,-3.4],[0,-3.4,0]]);

var out1b = MAP(BEZIER(S1)([c1b,c11b]))(domain);
var out2b = MAP(BEZIER(S1)([c2b,c22b]))(domain);
var out3b = MAP(BEZIER(S1)([c3b,c33b]))(domain);
var out4b = MAP(BEZIER(S1)([c4b,c44b]))(domain);
var out5b = MAP(BEZIER(S1)([c1b,c3b]))(domain);
var out6b = MAP(BEZIER(S1)([c2b,c4b]))(domain);
var out7b = MAP(BEZIER(S1)([c11b,c33b]))(domain);
var out8b = MAP(BEZIER(S1)([c22b,c44b]))(domain);

b = STRUCT ([out1b,out2b,out3b,out4b,out5b,out6b,out7b,out8b])

cuboid = CUBOID([0.5,5.5,0.2])
ct = T([0,2])([4.5,4.8])(cuboid)
ctb = T([0,2])([-0.5,4.8])(cuboid)

cuboid1 = CUBOID([0.2,4.95,0.2])
ct1 = T([0,1])([4.5,0.55])(cuboid1)
ct1b = T([0,1])([-0.2,0.55])(cuboid1)

var c5 = CUBIC_HERMITE(S0)([[4.7,0.15,0.4],[4.7,0.55,0],[0,0,-0.6],[0,0.6,0]]);
var c6 = CUBIC_HERMITE(S0)([[4.7,0.35,0.4],[4.7,0.55,0.2],[0,0,-0.2],[0,0.2,0]]);
out9 = MAP(c5)(dom)
out10 = MAP(c6)(dom)

var out9 = MAP(CYLINDRICAL_SURFACE(c5)([-0.2,0,0]))(domain);
var out10 = MAP(CYLINDRICAL_SURFACE(c6)([-0.2,0,0]))(domain);
var out11 = MAP(BEZIER(S1)([c5,c6]))(domain);
var out12 = T([0])([-0.2])(out11);
s1 = STRUCT([out9,out10,out11,out12])

cuboid2 = CUBOID([0.2,0.2,4.4])
ct2 = T([0,1,2])([4.5,0.15,0.4])(cuboid2)

cs = STRUCT([ct2,s1,ct1])
csb = T([0])([-4.7])(cs)

//sfera
var domain_sfera = DOMAIN([[0,2*PI],[0,2*PI]])([24,36]);
var sfera = function(r){
	return function(v){
	return [r*SIN(v[0])*COS(v[1]),r*SIN(v[0])*SIN(v[1]),r*COS(v[0])]
	}
}
var mapping = sfera(0.2);
var model = MAP(mapping)(domain_sfera);
var mt = T([0,1,2])([4.85,0.23,4.6])(model)
var mt1 = T([0,1,2])([-0.35,0.23,4.6])(model)
var mt2 = T([0,1,2])([4.25,0.23,2.1])(model)
var mt3 = T([0,1,2])([0.2,0.23,2.1])(model)

struct = COLOR([184/255, 115/255, 51/255])(STRUCT([cs,csb,ctb,ct,s,b,mt,mt1,mt2,mt3]))


//rifinitura cornicione
function wall_y_with_hole(dx,dy,dz,hole_y,hole_z,hole_dy,hole_dz) {
	var g1 = SIMPLEX_GRID (([[dx],[hole_y,-hole_dy, dy-hole_y-hole_dy],[dz]]))
	var g2 = SIMPLEX_GRID ([[dx],[-hole_y, hole_dy],[hole_z, -hole_dz,dz-hole_z-hole_dz]])
	var g = STRUCT([g1,g2]);
	return g;
}


function wall_y_with_N_hole(dx,dy,dz,hole_y,hole_z,hole_dy,hole_dz,n){
	var f = wall_y_with_hole(dx,dy,dz,hole_y,hole_z,hole_dy,hole_dz);
	var t = f;
	for(i=1; i<n; i++){
		t = STRUCT([t,T([1])([i*dy])(f)])
	}
	return t;
}

f1 = wall_y_with_N_hole(0.1,0.4,3,0.1,1,0.2,1.5,5)
f2 = wall_y_with_N_hole(0.1,0.4,1.4,0.1,0.5,0.2,0.2,5)
f1t = T([0,1,2])([-0.15,2.2,1.7])(f1)
f2t = T([0,1,2])([-0.15,2.2,0.3])(f2)

var domain = PROD1x1([INTERVALS(1)(14),INTERVALS(1)(14)]);
var dom = INTERVALS(1)(14)


function semicerchio(p1,p2,vettore){
	var var1 = CUBIC_HERMITE(S0)(p1);
	var var2 = CUBIC_HERMITE(S0)(p2);
	var v12 = MAP(BEZIER(S1)([var1,var2]))(domain);
	var v21 = T([0,1,2])(vettore)(v12);
	var var1l = MAP(CYLINDRICAL_SURFACE(var1)(vettore))(domain);
	var var2l = MAP(CYLINDRICAL_SURFACE(var2)(vettore))(domain);
	var s = STRUCT([v12,v21,var1l,var2l]);
	return s;
}

var p0 = [[0,1.5,0.4],[0,1.7,0.2],[0,0,-0.3],[0,0.3,0]]
var p00 = [[0,1.6,0.4],[0,1.7,0.3],[0,0,-0.2],[0,0.2,0]]
var sm0 = semicerchio(p0,p00,[-0.2,0,0])

var p1 = [[0,1.5,4.6],[0,1.7,4.8],[0,0,0.3],[0,0.3,0]]
var p11 = [[0,1.6,4.6],[0,1.7,4.7],[0,0,0.2],[0,0.2,0]]
var sm1 = semicerchio(p1,p11,[-0.2,0,0])

var t1 = CUBIC_HERMITE(S0)([[-0.05,1.6,0.4],[-0.05,1.7,0.3],[0,0,-0.2],[0,0.2,0]]);
var t0 = CUBIC_HERMITE(S0)([[-0.05,1.6,4.6],[-0.05,1.7,4.7],[0,0,0.2],[0,0.2,0]]);

var t10 = MAP(BEZIER(S1)([t1,t0]))(domain)
var t01 = T([0])([-0.1])(t10)

var cx =  CUBOID([0.2,0.1,4.2])
var cxt = T([0,1,2])([-0.2,1.5,0.4])(cx)
var cx1 = CUBOID([0.1,0.5,4.4])
var cx1t = T([0,1,2])([-0.15,1.7,0.3])(cx1)

var rif = STRUCT([sm0,sm1,cxt,cx1t, t10,t01])
var sr = R([0,1])([PI])(rif)
var srt = T([0,1])([-0.2,6.4])(sr)

var cx2 = CUBOID([0.2,3,0.1])
var cx2t = T([0,1,2])([-0.2,1.7,0.2])(cx2)
var cx2tb = T([0,1,2])([-0.2,1.7,4.7])(cx2)

struct_rif = STRUCT([rif,srt,cx2t,cx2tb,f1t,f2t])
rif2 = T([0])([4.7])(struct_rif)
rifinitura = STRUCT([struct_rif,rif2])

sedia = COLOR([184/255, 115/255, 51/255])(STRUCT([rifinitura,struct]))

//schienale
var p2 = [[0,5.5,3],[0.4,5.5,2.6],[0,0,-0.8],[0.8,0,0]]
var p22 = [[0.2,5.5,3],[0.4,5.5,2.8],[0,0,-0.4],[0.4,0,0]]
var sm2 = semicerchio(p2,p22,[0,0.2,0])

var p3 = [[0,5.5,9.7],[0.4,5.5,10.1],[0,0,0.8],[0.8,0,0]]
var p33 = [[0.2,5.5,9.7],[0.4,5.5,9.9],[0,0,0.4],[0.4,0,0]]
var sm3 = semicerchio(p3,p33,[0,0.2,0])

var cil = T([1,2])([5.5,3])(CUBOID([0.2,0.2,6.7]))



var sc = STRUCT([sm2,sm3,cil])
var scr = T([0,2])([4.3,12.7])(R([0,2])([PI])(sc))
var cil1 = T([0,1,2])([0.3,5.5,2.6])(CUBOID([3.7,0.2,0.2]))
var cil1t = T([2])([7.3])(cil1)

var schienale_ext = STRUCT([sc,scr,cil1,cil1t])

//parte interna
var p3 = [[0.75,5.5,3.1],[1.05,5.5,2.8],[0,0,-0.6],[0.6,0,0]]
var p33 = [[0.9,5.5,3.1],[1.05,5.5,2.95],[0,0,-0.3],[0.3,0,0]]
var sm3 = semicerchio(p3,p33,[0,0.2,0])

var p4 = [[0.75,5.5,9.6],[1.05,5.5,9.9],[0,0,0.6],[0.6,0,0]]
var p44 = [[0.9,5.5,9.6],[1.05,5.5,9.75],[0,0,0.3],[0.3,0,0]]
var sm4 = semicerchio(p4,p44,[0,0.2,0])

var cil2 = T([0,1,2])([0.75,5.5,3])(CUBOID([0.15,0.2,6.7]))

var ti1 = CUBIC_HERMITE(S0)([[0.9,5.575,3.1],[1.05,5.575,2.95],[0,0,-0.3],[0.3,0,0]]);
var ti0 = CUBIC_HERMITE(S0)([[0.9,5.575,9.6],[1.05,5.575,9.75],[0,0,0.3],[0.3,0,0]]);

var ti10 = MAP(BEZIER(S1)([ti1,ti0]))(domain)
var ti01 = T([1])([0.05])(ti10)

var sc1 = STRUCT([sm3,sm4,cil2,ti10,ti01])
var scr1 = T([0,2])([4.3,12.7])(R([0,2])([PI])(sc1))
var cil3 = T([0,1,2])([1.05,5.5,2.8])(CUBOID([2.4,0.2,0.15]))
var cil3t = T([2])([6.95])(cil3)
var interno1 = SIMPLEX_GRID([[-1.05,2.4],[-5.575,0.05],
	[-2.8,0.375,-0.5,0.23,-0.5,0.23,-0.5,0.23,-0.5,0.23,-0.5,0.23,-0.5,0.23,-0.5,0.23,-0.5,0.23,-0.5,0.385]]) 
var interno2 = SIMPLEX_GRID([[-1.05,0.5,-0.5,0.4,-0.5,0.5],[-5.575,0.05],[-2.8,7.1]]) 
var schienale_int = STRUCT([sc1,scr1,cil3,cil3t,interno1,interno2])

var schienale = COLOR([184/255, 115/255, 51/255])(T([1,2])([-0.4,2.1])(R([1,2])([-PI/9])(STRUCT([schienale_int,schienale_ext]))))

//sedile
sedile = CUBOID([4.5,5.5,0.3])
sed_t = T([2])([2.3])(sedile)

var l1 = CUBIC_HERMITE(S0)([[0,0.8,2.3],[0.8,0,2.3],[0,-1.6,0],[1.6,0,0]]);
var l2 = CUBIC_HERMITE(S0)([[0,4.7,2.3],[0.8,5.5,2.3],[0,1.6,0],[1.6,0,0]]);
var l3 = BEZIER(S0)([[0,0.8,2.3],[0,4.7,2.3]]);

var v12 = MAP(BEZIER(S1)([l1,l2]))(domain);
var v21 = T([0,1,2])([0,0,0.3])(v12);
var var1l = MAP(CYLINDRICAL_SURFACE(l1)([0,0,0.3]))(domain);
var var2l = MAP(CYLINDRICAL_SURFACE(l2)([0,0,0.3]))(domain);
var var3l = MAP(CYLINDRICAL_SURFACE(l2)([0,0,0.3]))(domain);
var s = STRUCT([v12,v21,var1l,var2l,var3l]);

var srt = T([0,1])([4.5,5.5])(R([0,1])(PI)(s))
var centro = T([0,2])([0.8,2.3])(CUBOID([2.9,5.5,0.3]))

var sedile = COLOR([184/255, 115/255, 51/255])(STRUCT([centro, srt,s]))


//funghetti
function cylinder(r,h){
	d = DISK(r)(36)
	c = EXTRUDE([h])(d)
	return c
}

var domain_fung = DOMAIN([[-PI/2, PI/2],[-PI/2,PI/2]])([24,36]);
var mapping = sfera(0.15);
var model1 = T([2])([0.1])(MAP(mapping)(domain_fung));
var d = cylinder(0.1,0.1);
var d1 = T([2])([0.1])(DISK(0.15)(24,36));
var fung = STRUCT([model1,d,d1]);


dt1 = T([0,1,2])([-0.25,5.5,5])(fung)
dt2 = T([1,2])([0.1,0.75])(R([1,2])(-PI/20)(dt1))
dt3 = T([0,1,2])([0.02,0.2,1.4])(R([1,2])(-2*PI/20)(dt1))
dt4 = T([0,1,2])([0.03,0.285,2])(R([1,2])(-3*PI/20)(dt1))
dt5 = T([0,1,2])([0.04,0.45,2.45])(R([1,2])(-4*PI/20)(dt1))
funghetti_left = STRUCT([dt1,dt2,dt3,dt4,dt5])
funghetti_right = T([0])([4.95])(funghetti_left)

funghetti = COLOR([184/255, 115/255, 51/255])(STRUCT([funghetti_left,funghetti_right]))

cilindro = cylinder(0.07,6)
cilindro_r = R([0,2])(PI/2)(cilindro)
cilindro_t = COLOR([0,0,0])(T([0,1,2])([-0.7,6.8,4.6])(cilindro_r))
sfera1 = COLOR([184/255, 115/255, 51/255])(T([0,1,2])([-0.7,6.8,4.6])(model))
sfera2 = COLOR([184/255, 115/255, 51/255])(T([0,1,2])([5.3,6.8,4.6])(model))
cil_sf = STRUCT([cilindro_t,sfera2,sfera1])

//cerniera
cil4 = cylinder(0.1,0.6)
cil4r = R([0,2])([PI/2])(cil4)
cil4t1 = T([0,1,2])([0.8,5.6,2.65])(cil4r)
cil4t2 = T([0,1,2])([3.1,5.6,2.65])(cil4r)
cerniera = COLOR([0.3,0.3,0.3])(STRUCT([cil4t2,cil4t1]))

//ruote
var domain_ruota = DOMAIN([[0,PI/2],[3*PI/8,5*PI/8]])([24,36]);
var mapping1 = sfera(0.3);
var model1 = COLOR([184/255, 115/255, 51/255])(MAP(mapping1)(domain_ruota));

var disk = DISK(0.15)(24)
var disk1 = EXTRUDE([0.2])(disk)
var disk_r = R([0,2])(PI/2)(disk1)
var disk_t = T([0,1,2])([-0.045,0.08,0.08])(disk_r)
var ruota = STRUCT([COLOR([0.3,0.3,0.3])(disk_t),model1])
var ruota1 = T([0,1])([-0.1,6.5])(ruota)
var ruota2 = T([0,1])([4.6,6.5])(ruota)
var ruote = STRUCT([ruota2,ruota1])
var sedia = STRUCT([sedia,sedile,schienale,funghetti,cil_sf,cerniera,ruote])
DRAW(sedia)