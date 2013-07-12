var dom1D = INTERVALS(1)(32);
var dom2D = PROD1x1([INTERVALS(1)(16),INTERVALS(1)(16)]);

//saucer_form
var Su0 = CUBIC_HERMITE(S0)([[COS(PI/6),SIN(PI/12),0],[COS(PI/6),-SIN(PI/12),0],[TAN(PI/6),0,0],[0,TAN(PI/6),0]]);
var Su1 = CUBIC_HERMITE(S0)([[0.5+2.3*COS(PI/6),0.1+2.3*SIN(PI/12),1],[0.5+2.3*COS(PI/6),-0.1 +2.3*(-SIN(PI/12)),1],[2.3*TAN(PI/6),0,0],[0,2.3*TAN(PI/6),0]]);
var Sv0 = BEZIER(S1)([[COS(PI/6),SIN(PI/12),0],[+COS(PI/6),+0.005+SIN(PI/12),0.5],[2.3*COS(PI/6),0.07+SIN(PI/12),0.1],[0.5+2.3*COS(PI/6),2.3*SIN(PI/12),1]]);
var Sv1 = BEZIER(S1)([[COS(PI/6),-SIN(PI/12),0],[+COS(PI/6),-0.005-SIN(PI/12),0.5],[2.3*COS(PI/6),-0.07-SIN(PI/12),0.1],[0.5+2.3*COS(PI/6),2.3*(-SIN(PI/12)),1]]);

var saucer_f = MAP(COONS_PATCH([Su0,Su1,Sv0,Sv1]))(dom2D);

var Su2 = CUBIC_HERMITE(S0)([[0.5+2.3*COS(PI/6),0.1+2.3*SIN(PI/12),1],[0.5+2.3*COS(PI/6),-0.1 +2.3*(-SIN(PI/12)),1],[2.3*TAN(PI/6),0,2],[0,2.3*TAN(PI/6),-2]]);
var Su12 = MAP(BEZIER(S1)([Su1,Su2]))(dom2D)

var saucer_f_black = COLOR([0,0,0])(R([0,1])(PI/6)(STRUCT([saucer_f,Su12])))
var saucer_f_white = COLOR([1,1,1])(STRUCT([saucer_f,Su12]))

saucer_form = STRUCT([saucer_f_black,saucer_f_white])

function rotate_form(partial,axes,angle,n){
	var form = partial
	var f = partial
	for (i=0; i<n; i++){
		f = R(axes)(angle)(f)
		form = STRUCT([form,f])
	}
	return form;
}
saucer_top = rotate_form(saucer_form,[0,1],PI/3,5)

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
arc1 = arc(PI/6,0,0.98);
arc_black = COLOR([0,0,0])(R([0,1])(PI/12)(arc1))
arc_white = COLOR([1,1,1])(R([0,1])(-PI/12)(arc1))
partial_saucer_bottom = STRUCT([arc_black,arc_white])
saucer_bottom = rotate_form(partial_saucer_bottom,[0,1],PI/3,5)
saucer = STRUCT([saucer_bottom,saucer_top])

//cup
var s0 = CUBIC_HERMITE(S0)([[COS(PI/6),SIN(PI/6),0],[COS(PI/6),-SIN(PI/6),0],[COS(PI/3),-SIN(PI/3),0],[-COS(PI/3),-SIN(PI/3),0]]);
var s1 = CUBIC_HERMITE(S0)([[2*COS(PI/6),2*SIN(PI/6),2],[2*COS(PI/6),2*(-SIN(PI/6)),2],[2*COS(PI/3),-2*SIN(PI/3),0],[-2*COS(PI/3),-2*SIN(PI/3),0]]);

var cup_f_white = COLOR([1,1,1])(MAP(BEZIER(S1)([s1,s0]))(dom2D))
var cup_f_black = COLOR([0,0,0])(R([0,1])([PI/3])(MAP(BEZIER(S1)([s1,s0]))(dom2D)))

var partial_cup_form = STRUCT([cup_f_white,cup_f_black])
var cup_form = T([2])([0.2])(rotate_form(partial_cup_form,[0,1],2*PI/3,3))

var arc2 = arc(PI/3,0,1);
var arc2_black = COLOR([0,0,0])(R([0,1])(PI/6)(arc2))
var arc2_white = COLOR([1,1,1])(R([0,1])(-PI/6)(arc2))
var partial_cup_bottom = STRUCT([arc2_white,arc2_black])
var cup_bottom = T([2])([0.2])(rotate_form(partial_cup_bottom,[0,1],2*PI/3,3))

//handle
var Su0 = CUBIC_HERMITE(S0)([[-1.19,0.15,0.3],[-1.19,-0.15,0.3],[-0.6,0,0],[0.6,0,0]]);
var Su1 = CUBIC_HERMITE(S0)([[-1.9,0.15,1],[-1.9,-0.15,1],[-0.6,0,0],[0.6,0,0]]);
var Sv0 = BEZIER(S1)([[-1.9,0.15,1],[-1.3,0.15,0.8],[-1.7,0.15,0.4],[-1.19,0.15,0.3]]);
var Sv1 = BEZIER(S1)([[-1.9,-0.15,1],[-1.3,-0.15,0.8],[-1.7,-0.15,0.4],[-1.19,-0.15,0.3]]);
var ext_h1 = COLOR([0,0,0])(MAP(COONS_PATCH([Su1,Su0,Sv0,Sv1]))(dom2D));

var Su2 = CUBIC_HERMITE(S0)([[-2.5,0.15,2.05],[-2.5,-0.15,2.05],[0,0,0.45],[0,0,-0.45]]);
var Sv2 = BEZIER(S1)([[-2.5,0.15,2.05],[-2.9,0.15,1.9],[-3.3,0.15,1.2],[-2.5,0.15,0.8]]);
var Sv3 = BEZIER(S1)([[-2.5,-0.15,2.05],[-2.9,-0.15,1.9],[-3.3,-0.15,1.2],[-2.5,-0.15,0.8]]);
var ext_h2 = COLOR([0,0,0])(MAP(COONS_PATCH([Su2,Su1,Sv2,Sv3]))(dom2D));

var Su3 = CUBIC_HERMITE(S0)([[-2,0.15,1.5],[-2,-0.15,1.5],[0.6,0,0],[-0.6,0,0]]);
var Sv4 = BEZIER(S1)([[-2,0.15,1.5],[-1.7,0.15,2],[-2.2,0.15,2.3],[-2.5,0.15,2.05]]);
var Sv5 = BEZIER(S1)([[-2,-0.15,1.5],[-1.7,-0.15,2],[-2.2,-0.15,2.3],[-2.5,-0.15,2.05]]);
var ext_h3 = COLOR([0,0,0])(MAP(COONS_PATCH([Su3,Su2,Sv4,Sv5]))(dom2D));

//internal_handle
var Su0_int = CUBIC_HERMITE(S0)([[-1.19,0.15,0.3],[-1.19,-0.15,0.3],[0.6,0,0],[-0.6,0,0]]);
var Su1_int = CUBIC_HERMITE(S0)([[-1.9,0.15,1],[-1.9,-0.15,1],[0.6,0,0],[-0.6,0,0]]);
var int_h1 = COLOR([1,1,1])(MAP(COONS_PATCH([Su1_int,Su0_int,Sv0,Sv1]))(dom2D));

var Su2_int = CUBIC_HERMITE(S0)([[-2.5,0.15,2.05],[-2.5,-0.15,2.05],[0,0,-0.45],[0,0,0.45]]);
var int_h2 = COLOR([1,1,1])(MAP(COONS_PATCH([Su2_int,Su1_int,Sv2,Sv3]))(dom2D));

var Su3_int = CUBIC_HERMITE(S0)([[-2,0.15,1.5],[-2,-0.15,1.5],[-0.6,0,0],[0.6,0,0]]);
var int_h3 = COLOR([1,1,1])(MAP(COONS_PATCH([Su3_int,Su2_int,Sv4,Sv5]))(dom2D));

var int_h4 = MAP(BEZIER(S1)([Su3,Su3_int]))(dom2D)
var int_h5 = MAP(BEZIER(S1)([Su0,Su0_int]))(dom2D)

var handle = STRUCT([ext_h1,ext_h2,ext_h3,int_h1,int_h2,int_h3,int_h4,int_h5])

var cup = STRUCT([saucer,cup_bottom,cup_form,handle])