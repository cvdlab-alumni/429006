var dom1D = INTERVALS(1)(14)
var dom2D = PROD1x1([INTERVALS(1)(14),INTERVALS(1)(14)]);
var light_brown = COLOR([194/255, 115/255, 51/255])
var brown = COLOR([150/255, 75/255, 0/255])

function cylinder(r,h){
	d = DISK(r)(36)
	c = EXTRUDE([h])(d)
	return c
}


function table_top(c0,c1,c00,c11,w,h,s){
	var cmap0 = MAP(CYLINDRICAL_SURFACE(c0)([0,0,h]))(dom2D);
	var cmap1 = MAP(CYLINDRICAL_SURFACE(c1)([0,0,h]))(dom2D);
	var cmap00 = MAP(CYLINDRICAL_SURFACE(c00)([0,0,h]))(dom2D);
	var cmap11 = MAP(CYLINDRICAL_SURFACE(c11)([0,0,h]))(dom2D);
	var cuboid = brown(T([0,1])([w-s,2*s])(CUBOID([s,6.1-1.2,h])))
	var cuboid1 = brown(T([0])([0.6])(CUBOID([w-1.2,s,h])))
	var cuboid2 = T([1])([6.1-s])(cuboid1)
	var cuboid3 = light_brown(T([0,1,2])([0.6,0.3,h-0.11])(CUBOID([w-1.2,6.1-0.6,0.01])))
	var bez1 = MAP(BEZIER(S1)([c0,c1]))(dom2D);
	var bez2 = MAP(BEZIER(S1)([c00,c11]))(dom2D);
	var bez3 = T([2])([h])(bez1);
	var bez4 = T([2])([h])(bez2);
	var bez5 = light_brown(T([2])([h-0.1])(MAP(BEZIER(S1)([c1,c11]))(dom2D)));
	var border_r0 = brown(STRUCT([bez1,bez2,bez3,bez4,cmap0,cmap1,cmap00,cmap11,cuboid]))
	var border_r = STRUCT([border_r0,bez5])
	var border_l = T([0,1])([4,6.1])(R([0,1])(PI)(border_r))
	var border = STRUCT([border_r,border_l,cuboid1,cuboid2,cuboid3])
	return border;
}

function table(punti0,punti1,punti2,n,w,d,s){
    var tavoli = []
	for (i=0;i<n;i++){
		var cyl = cylinder(0.1,6.8-i*s)
		var p0 = punti0.map(function(p){return[p[0]+i*s,p[1]+i*s,p[2]]});
		var p1 = punti1.map(function(p){return[p[0]-i*s,p[1]+i*s,p[2]]});
		var p2 = punti2.map(function(p){return[p[0]-i*s,p[1]+i*s,p[2]]});
		var c0 = CUBIC_HERMITE(S0)([p0[0],p0[1],[-s*2,0,0],[s*2,0,0]]);
		var c1 = CUBIC_HERMITE(S0)([p1[0],p1[1],[s*4,0,0],[0,s*4,0]]);
		var c2 = CUBIC_HERMITE(S0)([p2[0],p2[1],[s*2,0,0],[0,s*2,0]]);
		var p00 = punti0.map(function(p){return[p[0]+i*s,d-(p[1]+i*s),p[2]]});
		var p11 = punti1.map(function(p){return[p[0]-i*s,d-(p[1]+i*s),p[2]]});
		var p22 = punti2.map(function(p){return[p[0]-i*s,d-(p[1]+i*s),p[2]]});
		var p3 = p0.map(function(p){return[p[0]+i*s,p[1]+i*s,p[2]+s]});
		var p4 = p00.map(function(p){return[p[0]+i*s,p[1]+i*s,p[2]+s]});
		var c3 = CUBIC_HERMITE(S0)([p00[0],p00[1],[-s*2,0,0],[s*2,0,0]]);
		var c4 = CUBIC_HERMITE(S0)([p11[0],p11[1],[s*4,0,0],[0,-s*4,0]]);
		var c5 = CUBIC_HERMITE(S0)([p22[0],p22[1],[s*2,0,0],[0,-s*2,0]]);
		var c00 = BEZIER(S0)(p0);
		var c33 = BEZIER(S0)(p00);
		var cmap0 = MAP(CYLINDRICAL_SURFACE(c0)([0,0,s]))(dom2D);
		var cmap1 = MAP(CYLINDRICAL_SURFACE(c1)([0,0,s]))(dom2D);
		var cmap2 = MAP(CYLINDRICAL_SURFACE(c2)([0,0,s]))(dom2D);
		var cmap3 = MAP(CYLINDRICAL_SURFACE(c3)([0,0,s]))(dom2D);
		var cmap4 = MAP(CYLINDRICAL_SURFACE(c4)([0,0,s]))(dom2D);
		var cmap5 = MAP(CYLINDRICAL_SURFACE(c5)([0,0,s]))(dom2D);
		var bez1 = MAP(BEZIER(S1)([c1,c2]))(dom2D);
		var bez2 = MAP(BEZIER(S1)([c4,c5]))(dom2D);
		var bez3 = T([2])([s])(bez1);
		var bez4 = T([2])([s])(bez2);
		var bez5 = MAP(BEZIER(S1)([c0,c00]))(dom2D);
		var bez6 = MAP(BEZIER(S1)([c3,c33]))(dom2D);
		var bez7 = T([2])([s])(bez5);
		var bez8 = T([2])([s])(bez6);
		var cuboid0 = T([0,1])([i*s+s,i*s])(CUBOID([p1[0][0]-p0[0][0],s,s]))
		var cuboid1 = T([1])([p00[0][1]-i*s-s])(cuboid0)
		var cuboid2 = T([0,1])([w-i*s-s,i*s+2*s])(CUBOID([s,p22[1][1]-p2[1][1],s]))
		var basis = brown(STRUCT([cmap0,cmap1,cmap2,cmap3,cmap4,cmap5,bez1,bez2,bez3,bez4,bez5,bez6,bez7,bez8,cuboid0,cuboid1,cuboid2]))
		//table leg
		var leg1 = T([0,1])([i*s+s+0.1,i*s+0.15])(cyl)
		var leg2 = T([0])([p1[0][0]-p0[0][0]])(leg1)
		var leg12 = STRUCT([leg1,leg2])
		var leg4 = T([1])([p00[0][1]-i*s-s])(leg12)
		var table_leg = brown(STRUCT([leg12,leg4]))
		//
		if(i!=0){
			var cuboid3 = brown(T([0,1,2])([i*s+s,i*s,6.8-((i-1)*s)-s])(CUBOID([p1[0][0]-p0[0][0]+0.2,0.2,0.3])))
			var cuboid4 = brown(T([1])([p00[0][1]-i*s-s])(cuboid3))
			var cuboid45 = T([0,1,2])([p1[0][0]-p0[0][0]+i*s+s+0.1,i*s,6.8-((i-1)*s)-s])(CUBOID([0.2,p11[0][1]-p1[0][1]-0.1,0.3]))
			var cuboid55 = T([0,1,2])([i*s+s-0.1,i*s,6.8-((i-1)*s)-0.15])(CUBOID([0.2,p11[0][1]-p1[0][1]-0.1,0.15]))
			var cuboid5 = brown(STRUCT([cuboid45,cuboid55]))
			var cuboid6 = light_brown(T([0,1,2])([i*s+s+0.1,i*s+0.2,6.8-((i-1)*s)-0.1])(CUBOID([p1[0][0]-p0[0][0],p11[0][1]-p1[1][1]+0.1,0.05])))
		}
		else{
			var cuboid3 = brown(T([0,1,2])([i*s+s,i*s+0.1,5.9])(CUBOID([p1[0][0]-p0[0][0]+0.05,0.2,0.9])))
			var cuboid4 = brown(T([1])([p00[0][1]-i*s-s-0.05])(cuboid3))
			var cuboid5 = brown(T([0,1,2])([p1[0][0]-p0[0][0]+i*s+s,i*s+0.15,5.9])(CUBOID([0.2,p11[0][1]-p1[0][1]-0.2,0.9])))
			//da aggiustare la y
			var cuboid6 = T([2])([6.8])(table_top(c1,c2,c4,c5,w,0.6,s))
			}
		var cub = STRUCT([cuboid3,cuboid4,cuboid5,cuboid6])
		var table = STRUCT([basis,table_leg,cub])
		tavoli[i] = table
	}
	return tavoli;
}
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

