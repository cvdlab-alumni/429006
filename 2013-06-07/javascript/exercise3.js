//mountains
domain = DOMAIN([[0,20],[0,15]])([24,12])
//array dove memorizzo le posizioni randomiche delle montagne
var x = [];
var y = [];
var z = [];


var land = function(v){
		var a = v[0];
		var b = v[1];
		var u = a
		var v = b
		if (COS(a)*SIN(b)<0.1 && ((a>=0.02 && a<=7 && b>0.02 && b <14.98) 
			|| (a>7 && a<15 && b<8 && b>0.02) || (a > 16 && a<19.98 && b>10 && b<14.98)) ){
			    var w = 1 +2.5*Math.random();}
		else{
                var w =1};
        x.push(u)
        y.push(v)
        z.push(w)
		return [u,v,w];
}
depth = CUBOID([20,15,1])
mountain = MAP(land)(domain)
mountains = COLOR([205/255, 133/255, 63/255])(STRUCT([mountain,	depth]))

//lakes

l1 = T([0,1,2])([5,0,1.01])(CUBOID([4,3,0.1]))
l2 = T([0,1,2])([2.5,2.5,1.01])(CUBOID([3,3,0.1]))
l3 = T([0,1,2])([16.5,12.5,1.01])(CUBOID([3.5,2.5,0.1]))
lakes = COLOR([0,0.5,1])(STRUCT([l1,l2,l3]))


//trees

function cone(raggio,altezza){
	var domain = PROD1x1([INTERVALS(1)(20),INTERVALS(1)(6)]);
	var apex = [0,0,altezza];
	var funProfile = CUBIC_HERMITE(S0)([[raggio,0,0],[0,raggio,0],[0,2*raggio,0],[-2*raggio,0,0]]);
	var out1 = MAP(CONICAL_SURFACE(apex)(funProfile))(domain);
	var out2 = R([0,1])(PI)(out1)
	var out3 = STRUCT([out1,out2])
	var out4 = R([0,1])(PI/2)(out3)
	var out = STRUCT([out3,out4])
	return out;
}

function cylinder(r,h,discr){
	d = DISK(r)(discr)
	c = EXTRUDE([h])(d)
	return c
}

function tree(raggio,altezza,discr){
		base = COLOR([115/255, 74/255, 18/255])(cylinder(raggio/5,altezza,discr))
		cone1 = COLOR([0, 165/255, 80/255])(T([2])([2*altezza/3])(cone(raggio,2*altezza/3)))
		return STRUCT([base,cone1])
}
t = tree(0.25,0.8,20)
t1 = tree(0.3,0.9,20)


function forest(tree,x_in,x_fin,y_in,y_fin){
	s = ""
	for (i=0;i<x.length;i++){
		if((x[i]>=x_in && x[i]<x_fin) && (y[i]>=y_in && y[i]<y_fin))
					s = STRUCT([T([0,1,2])([x[i],y[i],z[i]])(tree),s])
			}
return s;
}
f1 = forest(t,0.2,3,0.2,14.98)
f2 = forest(t1,8,14,0.2,8)
f3 = forest(t,4,7.5,5.5,13)

forests = STRUCT([f1,f2,f3])

model = STRUCT([mountains,lakes,forests])
DRAW(model)
