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

model = STRUCT([mountains,lakes])
DRAW(model)
