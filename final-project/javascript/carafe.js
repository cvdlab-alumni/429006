var dom1D = INTERVALS(1)(32);
var dom2D = PROD1x1([INTERVALS(1)(16),INTERVALS(1)(16)]);

function move(p,h,s){
	var p1 = p;
	p1[0] = [p[0][0]+s[0],p[0][1]+s[1],p[0][2]+h];
	p1[1] = [p[1][0]+s[0],p[1][1]+s[1],p[1][2]+h];
	return p1;
}

function side_form(p1,p2,s){
	var c0 = CUBIC_HERMITE(S0)(p1)
	var c1 = CUBIC_HERMITE(S0)(move(p1,0.5,[0,s*0.5,0]))
	var c2 = CUBIC_HERMITE(S0)(move(p1,1.8,[0,s*0.3,0]))
	var c3 = BEZIER(S0)(p2)
	var b2 = MAP(BEZIER(S1)([c0,c1,c2,c3]))(dom2D)
	return b2;
}

function rotate(punti,v){
	var p1 = punti;
	var p1 = p1.map(function(p){return[p[0]*v[0],p[1]*v[1],p[2]*v[2]]});
	return p1;
}


//carafe's form

//body
//bezier point
var b0 = [[1.5,0,2.9],[COS(PI/10)+0.1,-SIN(PI/10)+0.08,2.95],[COS(PI/10),-SIN(PI/10),3]]
var b1 = [[COS(PI/10),-SIN(PI/10),3],[COS(PI/10)-0.05,-SIN(PI/10)-0.15,3.04],[COS(3*PI/10)-0.05,-SIN(3*PI/10)-0.15,3.18],[COS(3*PI/10),-SIN(3*PI/10),3.1]]
var b2 = [[COS(3*PI/10),-SIN(3*PI/10),3.1],[COS(3*PI/10)-0.05,-SIN(3*PI/10)-0.05,3.2],[COS(PI/2)-0.05,-SIN(PI/2)-0.05,3.05],[COS(PI/2),-SIN(PI/2),3]]
var b3 = [[COS(PI/2),-SIN(PI/2),3],[COS(PI/2)-0.05,-SIN(PI/2)-0.05,3.05],[COS(7*PI/10)-0.05,-SIN(7*PI/10)-0.05,2.9],[COS(7*PI/10),-SIN(7*PI/10),2.82]]
var b4 = [[COS(7*PI/10),-SIN(7*PI/10),2.82],[COS(7*PI/10)-0.3,-SIN(7*PI/10)+0.2,2.8],[COS(9*PI/10)-0.05,-SIN(9*PI/10)+0.2,2.7],[COS(9*PI/10),-SIN(9*PI/10),2.75]]
var b5 = [[COS(9*PI/10),-SIN(9*PI/10),2.75],[COS(9*PI/10)-0.1,-SIN(9*PI/10)+0.2,2.8],[COS(PI),-SIN(PI),2.8]]
var b6 = [[1,0,0],[1.9,0,1],[1,0,1.8],[1.5,0,2.9]]

//hermite point
var p0 = [[1,0,0],[COS(PI/10),-SIN(PI/10),0],[COS(PI/10)/2,0,0],[-SIN(PI/10)/3,-2*COS(PI/10)/3,0]]
var p01 = [[COS(PI/10),-SIN(PI/10),0],[COS(PI/10),-SIN(PI/10),2.95],[0,-1.95,0],[0,1.95,0]]
var p1 = [[COS(PI/10),-SIN(PI/10),0],[COS(3*PI/10),-SIN(3*PI/10),0],[-SIN(PI/10)/3,-2*COS(PI/10)/3,0],[-2*SIN(PI/5)/3,-COS(PI/5)/3,0]]
var p2 = [[COS(3*PI/10),-SIN(3*PI/10),0],[COS(PI/2),-SIN(PI/2),0],[-2*SIN(PI/5)/3,-COS(PI/5)/3,0],[-SIN(PI/5),0,0]]
var p3 = [[COS(PI/2),-SIN(PI/2),0],[COS(7*PI/10),-SIN(7*PI/10),0],[-SIN(PI/5),0,0],[-2*SIN(PI/5)/3,COS(PI/5)/3,0]]
var p4 = [[COS(7*PI/10),-SIN(7*PI/10),0],[COS(9*PI/10),-SIN(9*PI/10),0],[-2*SIN(PI/5)/3,COS(PI/5)/3,0],[-SIN(PI/10)/3,2*COS(PI/10)/3,0]]
var p5 = [[COS(9*PI/10),-SIN(9*PI/10),0],[COS(PI),-SIN(PI),0],[-SIN(PI/10)/2,2*COS(PI/10)/3,0],[0,SIN(PI/10),0]]
var p6 = [[-1,0,0],[-1,0,2.9],[-0.5,0,0],[0.5,0,0]]
var p7 = [[COS(9*PI/10),-SIN(9*PI/10),0],[COS(9*PI/10),-SIN(9*PI/10),2.75],[0,-2,0],[0,2,0]]

//rotate point
var p0r = rotate(p0,[1,-1,1])
var b0r = rotate(b0,[1,-1,1])
var p01r= rotate(p01,[1,-1,1])
var p1r = rotate(p1,[1,-1,1])
var b1r = rotate(b1,[1,-1,1])
var p2r = rotate(p2,[1,-1,1])
var b2r = rotate(b2,[1,-1,1])
var p3r = rotate(p3,[1,-1,1])
var b3r = rotate(b3,[1,-1,1])
var p4r = rotate(p4,[1,-1,1])
var b4r = rotate(b4,[1,-1,1])
var p5r = rotate(p5,[1,-1,1])
var b5r = rotate(b5,[1,-1,1])
var p7r = rotate(p7,[1,-1,1])

//coons_patch
//front part
var Su0 = CUBIC_HERMITE(S0)(p0)
var Su1 = BEZIER(S0)(b0)
var Sv0 = CUBIC_HERMITE(S1)(p01)
var Sv1 = BEZIER(S1)(b6)
var spout0 = MAP(COONS_PATCH([Su0,Su1,Sv1,Sv0]))(dom2D);

var Su00 = CUBIC_HERMITE(S0)(p0r)
var Su11 = BEZIER(S0)(b0r)
var Sv00 = CUBIC_HERMITE(S1)(p01r)
var Sv11 = BEZIER(S1)(b6)

var spout1 = MAP(COONS_PATCH([Su00,Su11,Sv11,Sv00]))(dom2D)
var spout = COLOR([1,1,1])(STRUCT([spout0,spout1]))

//behind part
var Su2 = CUBIC_HERMITE(S0)(p5)
var Su3 = BEZIER(S0)(b5)
var Sv22 = CUBIC_HERMITE(S1)(p6)
var Sv33 = CUBIC_HERMITE(S1)(p7)

var behindp0 = MAP(COONS_PATCH([Su2,Su3,Sv33,Sv22]))(dom2D)

var Su2 = CUBIC_HERMITE(S0)(p5r)
var Su3 = BEZIER(S0)(b5r)
var Sv22 = CUBIC_HERMITE(S1)(p6)
var Sv33 = CUBIC_HERMITE(S1)(p7r)

var behindp1 = MAP(COONS_PATCH([Su2,Su3,Sv33,Sv22]))(dom2D)
var behind_part = COLOR([0,0,0])(STRUCT([behindp0,behindp1]))

//side
var form1 = COLOR([0,0,0])(side_form(p1,b1,-1))
var form11 = COLOR([0,0,0])(side_form(p1r,b1r,1))
var form2 = COLOR([1,1,1])(side_form(p2,b2,-1))
var form22 = COLOR([1,1,1])(side_form(p2r,b2r,1))
var form3 = COLOR([0,0,0])(side_form(p3,b3,-1))
var form33 = COLOR([0,0,0])(side_form(p3r,b3r,1))
var form4 = COLOR([1,1,1])(side_form(p4,b4,-1))
var form44 = COLOR([1,1,1])(side_form(p4r,b4r,1))

var carafe_form = STRUCT([spout,form1,form11,form2,form22,form3,form33,form4,form44,behind_part])

//handle
var hSu0 = CUBIC_HERMITE(S0)([[-1,0.18,0.5],[-1,-0.15,0.5],[0,0.6,0],[0,-0.6,0]]);
var hSu1 = CUBIC_HERMITE(S0)([[-1.8,0.15,1.8],[-1.8,-0.15,1.8],[-0.6,0,0],[0.6,0,0]]);
var hSv0 = BEZIER(S1)([[-1.8,0.15,1.8],[-1.3,0.15,1.5],[-1.5,0.15,0.3],[-1,0.15,0.5]]);
var hSv1 = BEZIER(S1)([[-1.8,-0.15,1.8],[-1.3,-0.15,1.5],[-1.5,-0.15,0.3],[-1,-0.15,0.5]]);
var hext_h1 = COLOR([0,0,0])(MAP(COONS_PATCH([hSu1,hSu0,hSv0,hSv1]))(dom2D));

var hSu2 = CUBIC_HERMITE(S0)([[-1.8,0.15,3.2],[-1.8,-0.15,3.2],[0,0,0.45],[0,0,-0.45]]);
var hSv2 = BEZIER(S1)([[-1.8,0.15,3.2],[-2.2,0.15,3.1],[-2.6,0.15,2.4],[-1.8,0.15,1.8]]);
var hSv3 = BEZIER(S1)([[-1.8,-0.15,3.2],[-2.2,-0.15,3.1],[-2.6,-0.15,2.4],[-1.8,-0.15,1.8]]);
var hext_h2 = COLOR([0,0,0])(MAP(COONS_PATCH([hSu2,hSu1,hSv2,hSv3]))(dom2D));

var hSu3 = CUBIC_HERMITE(S0)([[-1.23,0.15,2.5],[-1.23,-0.15,2.5],[0.6,0,0],[-0.6,0,0]]);
var hSv4 = BEZIER(S1)([[-1.23,0.15,2.5],[-0.9,0.15,3],[-1.5,0.15,3.3],[-1.8,0.15,3.2]]);
var hSv5 = BEZIER(S1)([[-1.23,-0.15,2.5],[-0.9,-0.15,3],[-1.5,-0.15,3.3],[-1.8,-0.15,3.2]]);
var hext_h3 = COLOR([0,0,0])(MAP(COONS_PATCH([hSu3,hSu2,hSv4,hSv5]))(dom2D));

//internal_handle
var hSu0_int = CUBIC_HERMITE(S0)([[-1,0.15,0.5],[-1,-0.15,0.5],[0,0.6,0],[0,-0.6,0]]);
var hSu1_int = CUBIC_HERMITE(S0)([[-1.8,0.15,1.8],[-1.8,-0.15,1.8],[0.6,0,0],[-0.6,0,0]]);
var hint_h1 = COLOR([1,1,1])(MAP(COONS_PATCH([hSu1_int,hSu0_int,hSv0,hSv1]))(dom2D));

var hSu2_int = CUBIC_HERMITE(S0)([[-1.8,0.15,3.2],[-1.8,-0.15,3.2],[0,0,-0.45],[0,0,0.45]]);
var hint_h2 = COLOR([1,1,1])(MAP(COONS_PATCH([hSu2_int,hSu1_int,hSv2,hSv3]))(dom2D));

var hSu3_int = CUBIC_HERMITE(S0)([[-1.23,0.15,2.5],[-1.23,-0.15,2.5],[-0.6,0,0],[0.6,0,0]]);
var hint_h3 = COLOR([1,1,1])(MAP(COONS_PATCH([hSu3_int,hSu2_int,hSv4,hSv5]))(dom2D));

var hint_h4 = MAP(BEZIER(S1)([hSu3,hSu3_int]))(dom2D)
var hfondo = COLOR([0,0,0])(MAP(BEZIER(S1)([hSu0_int,hSu0]))(dom2D))
var handle = STRUCT([hext_h1,hext_h2,hext_h3,hint_h1,hint_h2,hint_h3,hint_h4,hfondo])

var Su0 = CUBIC_HERMITE(S0)([[COS(PI/10),-SIN(PI/10),0],[COS(3*PI/10),-SIN(3*PI/10),0],[-SIN(PI/10)/3,-2*COS(PI/10)/3,0],[-2*SIN(PI/5)/3,-COS(PI/5)/3,0]])
var Su4 = CUBIC_HERMITE(S0)([[1.1*COS(PI/10),-1.1*SIN(PI/10),-0.5],[1.1*COS(3*PI/10),-1.1*SIN(3*PI/10),-0.5],[-2*SIN(PI/10)/6,-4*COS(PI/10)/6,0],[-4*SIN(PI/5)/6,-2*COS(PI/5)/6,0]])
var Sv4 = BEZIER(S1)([[1.1*COS(PI/10),-1.1*SIN(PI/10),-0.5],[COS(PI/10),-SIN(PI/10),-0.3],[COS(PI/10),-SIN(PI/10),-0.1],[COS(PI/10),-SIN(PI/10),0]]);
var Sv5 = BEZIER(S1)([[1.1*COS(3*PI/10),-1.1*SIN(3*PI/10),-0.5],[COS(3*PI/10),-SIN(3*PI/10),-0.3],[COS(3*PI/10),-SIN(3*PI/10),-0.1],[COS(3*PI/10),-SIN(3*PI/10),0]]);
var tp4 = MAP(COONS_PATCH([Su4,Su0,Sv4,Sv5]))(dom2D);


function arc(alpha,r,R){
	var domain = DOMAIN([[0,alpha], [r,R]])([36,1]);
	var mapping = function(v){
		var a = v[0];
		var r = v[1];
		return [r*COS(a), r*SIN(a)];
	};
	var model = MAP(mapping)(domain);
	return model;
};
var arc1 = arc(PI/5,0,1.1);
var arc_black = COLOR([0,0,0])(R([0,1])(PI/10)(arc1))
var arc_white = COLOR([1,1,1])(R([0,1])(-PI/10)(arc1))
var partial_sb_bottom = STRUCT([arc_black,arc_white])
var sb_bottom = rotate_form(partial_sb_bottom,[0,1],2*PI/5,5)
var sb = T([2])([-0.5])(sb_bottom)


function rotate_form(partial,axes,angle,n){
	var form = partial
	var f = partial
	for (i=0; i<n; i++){
		f = R(axes)(angle)(f)
		form = STRUCT([form,f])
	}
	return form;
}

var partial_tp_black = COLOR([0,0,0])(tp4)
var partial_tp_white = COLOR([1,1,1])(R([0,1])(PI/5)(tp4))
var partial_tp = STRUCT([partial_tp_white,partial_tp_black])
var carafebottom = rotate_form(partial_tp,[0,1],2*PI/5,4)

var carafe = STRUCT([carafe_form,handle,carafebottom,sb])