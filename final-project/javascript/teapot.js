var dom1D = INTERVALS(1)(32);
var dom2D = PROD1x1([INTERVALS(1)(16),INTERVALS(1)(16)]);

//teapot's form
var Su0 = CUBIC_HERMITE(S0)([[COS(PI/6),SIN(PI/12),0],[COS(PI/6),-SIN(PI/12),0],[TAN(PI/6),0,0],[0,TAN(PI/6),0]]);
var Su1 = CUBIC_HERMITE(S0)([[0.2+COS(PI/6),SIN(PI/12),1.9],[0.2+COS(PI/6),(-SIN(PI/12)),1.9],[TAN(PI/6),0,0],[0,TAN(PI/6),0]]);
var Sv0 = BEZIER(S1)([[COS(PI/6),SIN(PI/12),0],[COS(PI/6)+0.4,SIN(PI/12)+0.2,0.6],[COS(PI/6)+0.4,SIN(PI/12)+0.2,1.3],[COS(PI/6),SIN(PI/12),1.9]]);
var Sv1 = BEZIER(S1)([[COS(PI/6),-SIN(PI/12),0],[COS(PI/6)+0.4,-SIN(PI/12)-0.2,0.6],[COS(PI/6)+0.4,-SIN(PI/12)-0.2,1.3],[COS(PI/6),-SIN(PI/12),1.9]]);
var tp1 = MAP(COONS_PATCH([Su0,Su1,Sv0,Sv1]))(dom2D);

//top
var Su2 = CUBIC_HERMITE(S0)([[COS(PI/6),SIN(PI/12),3.1],[COS(PI/6),-SIN(PI/12),3.1],[TAN(PI/6),0,0],[0,TAN(PI/6),0]]);
var Su3 = CUBIC_HERMITE(S0)([[0.1+COS(PI/6),SIN(PI/12),2.1],[0.1+COS(PI/6),(-SIN(PI/12)),2.1],[TAN(PI/6),0,0],[0,TAN(PI/6),0]]);
var Sv2 = BEZIER(S1)([[0.1+COS(PI/6),SIN(PI/12),2.1],[COS(PI/6)+0.2,SIN(PI/12)+0.1,2.4],[COS(PI/6)+0.2,SIN(PI/12)+0.1,2.8],[COS(PI/6),SIN(PI/12),3.1]]);
var Sv3 = BEZIER(S1)([[0.1+COS(PI/6),-SIN(PI/12),2.1],[COS(PI/6)+0.2,-SIN(PI/12)-0.1,2.4],[COS(PI/6)+0.2,-SIN(PI/12)-0.1,2.8],[COS(PI/6),-SIN(PI/12),3.1]]);
var tp2 = MAP(COONS_PATCH([Su1,Su2,Sv2,Sv3]))(dom2D);

//center part
var Sv4 = BEZIER(S1)([[0.2+COS(PI/6),SIN(PI/12),1.9],[COS(PI/6)+0.15,SIN(PI/12),1.95],[COS(PI/6)+0.15,SIN(PI/12),2],[0.1+COS(PI/6),SIN(PI/12),2.1]]);
var Sv5 = BEZIER(S1)([[0.2+COS(PI/6),-SIN(PI/12),1.9],[COS(PI/6)+0.15,-SIN(PI/12),1.95],[COS(PI/6)+0.15,-SIN(PI/12),2],[0.1+COS(PI/6),-SIN(PI/12),2.1]]);
var tp3 = MAP(COONS_PATCH([Su1,Su3,Sv4,Sv5]))(dom2D);

//bottom
var Su4 = CUBIC_HERMITE(S0)([[COS(PI/6),SIN(PI/12),-0.5],[COS(PI/6),-SIN(PI/12),-0.5],[TAN(PI/6),0,0],[0,TAN(PI/6),0]]);
var Sv6 = BEZIER(S1)([[0.1+COS(PI/6),SIN(PI/12),-0.5],[COS(PI/6)+0.15,SIN(PI/12),-0.3],[COS(PI/6)+0.15,SIN(PI/12),-0.1],[0.1+COS(PI/6),SIN(PI/12),0]]);
var Sv7 = BEZIER(S1)([[0.1+COS(PI/6),-SIN(PI/12),-0.5],[COS(PI/6)+0.15,-SIN(PI/12),-0.3],[COS(PI/6)+0.15,-SIN(PI/12),-0.1],[0.1+COS(PI/6),-SIN(PI/12),0]]);
var tp4 = MAP(COONS_PATCH([Su4,Su0,Sv6,Sv7]))(dom2D);

partial_teapot_form = STRUCT([tp1,tp2,tp3,tp4])

partial_tp_black = COLOR([0,0,0])(partial_teapot_form)
partial_tp_white = COLOR([1,1,1])(R([0,1])(PI/6)(partial_teapot_form))

partial_tp = STRUCT([partial_tp_white,partial_tp_black])
function rotate_form(partial,axes,angle,n){
	var form = partial
	var f = partial
	for (i=0; i<n; i++){
		f = R(axes)(angle)(f)
		form = STRUCT([form,f])
	}
	return form;
}
teapot_form = rotate_form(partial_tp,[0,1],PI/3,5)

//knob
k0 = CUBIC_HERMITE(S1)([[-0.15,0,0],[0.15,0,0.5],[0,0,0.5],[0.5,0,0]])
k1 = CUBIC_HERMITE(S1)([[0,0,0],[0.15,0,0.3],[0,0,0.3],[0.3,0,0]])
k2 = CUBIC_HERMITE(S0)([[-0.15,0,0],[0,0,0],[0,-0.15,0],[0,0.15,0]])
k3 = CUBIC_HERMITE(S0)([[0.15,0,0.5],[0.15,0,0.3],[0,-0.15,0],[0,0.15,0]])
knob_right = MAP(COONS_PATCH([k2,k3,k0,k1]))(dom2D);

k4 = CUBIC_HERMITE(S0)([[-0.15,0,0],[0,0,0],[0,0.15,0],[0,-0.15,0]])
k5 = CUBIC_HERMITE(S0)([[0.15,0,0.5],[0.15,0,0.3],[0,0.15,0],[0,-0.15,0]])
knob_left = MAP(COONS_PATCH([k4,k5,k0,k1]))(dom2D);

knob = COLOR([0,0,0])(T([2])([1])(STRUCT([knob_left,knob_right])))

var domain_sphere = DOMAIN([[0,PI/3],[-PI/12,PI/12]])([24,36]);
var sphere = function(r){
	return function(v){
	return [r*SIN(v[0])*COS(v[1]),r*SIN(v[0])*SIN(v[1]),r*COS(v[0])]
	}
}
var mapping = sphere(1.05);
var partial_cover = MAP(mapping)(domain_sphere);
var partial_cover_black = COLOR([0,0,0])(partial_cover)
var partial_cover_white = COLOR([1,1,1])(R([0,1])(PI/6)(partial_cover))
var partial_cover_bw = STRUCT([partial_cover_white,partial_cover_black])

cover = T([2])([2.6])(STRUCT([rotate_form(partial_cover_bw,[0,1],PI/3,5),knob]))

//handle
//external_handle

var Su0 = CUBIC_HERMITE(S0)([[-1,0.15,0.5],[-1,-0.15,0.5],[-0.6,0,0],[0.6,0,0]]);
var Su1 = CUBIC_HERMITE(S0)([[-1.8,0.15,1.8],[-1.8,-0.15,1.8],[-0.6,0,0],[0.6,0,0]]);
var Sv0 = BEZIER(S1)([[-1.8,0.15,1.8],[-1.3,0.15,1.5],[-1.5,0.15,0.3],[-1,0.15,0.5]]);
var Sv1 = BEZIER(S1)([[-1.8,-0.15,1.8],[-1.3,-0.15,1.5],[-1.5,-0.15,0.3],[-1,-0.15,0.5]]);
var ext_h1 = COLOR([0,0,0])(MAP(COONS_PATCH([Su1,Su0,Sv0,Sv1]))(dom2D));

var Su2 = CUBIC_HERMITE(S0)([[-1.8,0.15,3.2],[-1.8,-0.15,3.2],[0,0,0.45],[0,0,-0.45]]);
var Sv2 = BEZIER(S1)([[-1.8,0.15,3.2],[-2.2,0.15,3.1],[-2.6,0.15,2.4],[-1.8,0.15,1.8]]);
var Sv3 = BEZIER(S1)([[-1.8,-0.15,3.2],[-2.2,-0.15,3.1],[-2.6,-0.15,2.4],[-1.8,-0.15,1.8]]);
var ext_h2 = COLOR([0,0,0])(MAP(COONS_PATCH([Su2,Su1,Sv2,Sv3]))(dom2D));

var Su3 = CUBIC_HERMITE(S0)([[-1.2,0.15,2.5],[-1.2,-0.15,2.5],[0.6,0,0],[-0.6,0,0]]);
var Sv4 = BEZIER(S1)([[-1.2,0.15,2.5],[-0.9,0.15,3],[-1.5,0.15,3.3],[-1.8,0.15,3.2]]);
var Sv5 = BEZIER(S1)([[-1.2,-0.15,2.5],[-0.9,-0.15,3],[-1.5,-0.15,3.3],[-1.8,-0.15,3.2]]);
var ext_h3 = COLOR([0,0,0])(MAP(COONS_PATCH([Su3,Su2,Sv4,Sv5]))(dom2D));

//internal_handle
var Su0_int = CUBIC_HERMITE(S0)([[-1,0.15,0.5],[-1,-0.15,0.5],[0.6,0,0],[-0.6,0,0]]);
var Su1_int = CUBIC_HERMITE(S0)([[-1.8,0.15,1.8],[-1.8,-0.15,1.8],[0.6,0,0],[-0.6,0,0]]);
var int_h1 = COLOR([1,1,1])(MAP(COONS_PATCH([Su1_int,Su0_int,Sv0,Sv1]))(dom2D));

var Su2_int = CUBIC_HERMITE(S0)([[-1.8,0.15,3.2],[-1.8,-0.15,3.2],[0,0,-0.45],[0,0,0.45]]);
var int_h2 = COLOR([1,1,1])(MAP(COONS_PATCH([Su2_int,Su1_int,Sv2,Sv3]))(dom2D));

var Su3_int = CUBIC_HERMITE(S0)([[-1.2,0.15,2.5],[-1.2,-0.15,2.5],[-0.6,0,0],[0.6,0,0]]);
var int_h3 = COLOR([1,1,1])(MAP(COONS_PATCH([Su3_int,Su2_int,Sv4,Sv5]))(dom2D));

var int_h4 = MAP(BEZIER(S1)([Su3,Su3_int]))(dom2D)
var handle = STRUCT([ext_h1,ext_h2,ext_h3,int_h1,int_h2,int_h3,int_h4])

//spout
var sv1 = BEZIER(S1)([[1,0,0.2],[3,0,0.2],[1.5,0,2.8],[2.5,0,2.8]])
var sv2 = BEZIER(S1)([[1,0,1.85],[2,0,1.75],[1.5,0,3.1],[2.3,0,3.1]])
var su1 = BEZIER(S0)([[1,0,0.2],[1.2,0.35,0.7],[1.2,0.35,1.2],[1,0,1.85]])
var su2 = BEZIER(S0)([[2.5,0,2.8],[2.4,0.2,2.9],[2.1,0.2,3],[2.3,0,3.1]])
var spout1 = MAP(COONS_PATCH([su1,su2,sv1,sv2]))(dom2D)

var su3 = BEZIER(S0)([[1,0,0],[1.5,-0.35,0.7],[1.5,-0.35,1.2],[1,0,1.85]])
var su4 = BEZIER(S0)([[2.5,0,2.8],[2.4,-0.2,2.9],[2.1,-0.2,3],[2.3,0,3.1]])
var spout2 = MAP(COONS_PATCH([su3,su4,sv1,sv2]))(dom2D)

var spout = COLOR([1,1,1])(STRUCT([spout1,spout2]))

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
var arc1 = arc(PI/6,0,1.03);
var arc_white = COLOR([1,1,1])(R([0,1])(PI/12)(arc1))
var arc_black = COLOR([0,0,0])(R([0,1])(-PI/12)(arc1))
var partial_sb_bottom = STRUCT([arc_black,arc_white])
var sb_bottom = rotate_form(partial_sb_bottom,[0,1],PI/3,5)
var sb = T([2])([-0.5])(sb_bottom)


teapot = STRUCT([teapot_form,cover,handle,spout,sb])
DRAW(teapot)