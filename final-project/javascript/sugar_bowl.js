var dom1D = INTERVALS(1)(32);
var dom2D = PROD1x1([INTERVALS(1)(16),INTERVALS(1)(16)]);

//sugar bowl
var Su0 = CUBIC_HERMITE(S0)([[COS(PI/6),SIN(PI/12),0],[COS(PI/6),-SIN(PI/12),0],[TAN(PI/6),0,0],[0,TAN(PI/6),0]]);
var Su1 = CUBIC_HERMITE(S0)([[0.2+COS(PI/6),SIN(PI/12),1.9],[0.2+COS(PI/6),(-SIN(PI/12)),1.9],[TAN(PI/6),0,0],[0,TAN(PI/6),0]]);
var Sv0 = CUBIC_HERMITE(S0)([[2*COS(PI/6),2*SIN(PI/12),1],[2*COS(PI/6),-2*SIN(PI/12),1],[2*TAN(PI/6),0,0],[0,2*TAN(PI/6),0]]);
var sb1 = MAP(BEZIER(S1)([Su1,Sv0,Su0]))(dom2D)

//sugar bowl bottom
var Su2 = CUBIC_HERMITE(S0)([[COS(PI/6),SIN(PI/12),-0.5],[COS(PI/6),-SIN(PI/12),-0.5],[TAN(PI/6),0,0],[0,TAN(PI/6),0]]);
var Sv2 = CUBIC_HERMITE(S0)([[1.4*COS(PI/6),1.4*SIN(PI/12),-0.5],[1.4*COS(PI/6),-1.4*SIN(PI/12),-0.5],[1.4*TAN(PI/6),0,0],[0,1.4*TAN(PI/6),0]]);
var sb2 = MAP(BEZIER(S1)([Su2,Sv2,Su0]))(dom2D)

var partial_sugarbowl_form = STRUCT([sb1,sb2])

var partial_sb_black = COLOR([0,0,0])(partial_sugarbowl_form)
var partial_sb_white = COLOR([1,1,1])(R([0,1])(PI/6)(partial_sugarbowl_form))

var partial_sb = STRUCT([partial_sb_white,partial_sb_black])
function rotate_form(partial,axes,angle,n){
	var form = partial
	var f = partial
	for (i=0; i<n; i++){
		f = R(axes)(angle)(f)
		form = STRUCT([form,f])
	}
	return form;
}
var sugarbowl_form = rotate_form(partial_sb,[0,1],PI/3,5)

//knob
var k0 = CUBIC_HERMITE(S1)([[-0.15,0,0],[0.15,0,0.5],[0,0,0.5],[0.5,0,0]])
var k1 = CUBIC_HERMITE(S1)([[0,0,0],[0.15,0,0.3],[0,0,0.3],[0.3,0,0]])
var k2 = CUBIC_HERMITE(S0)([[-0.15,0,0],[0,0,0],[0,-0.15,0],[0,0.15,0]])
var k3 = CUBIC_HERMITE(S0)([[0.15,0,0.5],[0.15,0,0.3],[0,-0.15,0],[0,0.15,0]])
var knob_right = MAP(COONS_PATCH([k2,k3,k0,k1]))(dom2D);

var k4 = CUBIC_HERMITE(S0)([[-0.15,0,0],[0,0,0],[0,0.15,0],[0,-0.15,0]])
var k5 = CUBIC_HERMITE(S0)([[0.15,0,0.5],[0.15,0,0.3],[0,0.15,0],[0,-0.15,0]])
var knob_left = MAP(COONS_PATCH([k4,k5,k0,k1]))(dom2D);
 
var knob = COLOR([0,0,0])(T([2])([1.15])(STRUCT([knob_left,knob_right])))

var domain_sphere = DOMAIN([[0,PI/3],[-PI/12,PI/12]])([24,36]);
var sphere = function(r){
	return function(v){
	return [r*SIN(v[0])*COS(v[1]),r*SIN(v[0])*SIN(v[1]),r*COS(v[0])]
	}
}
var mapping = sphere(1.3);
var partial_cover = MAP(mapping)(domain_sphere);
var partial_cover_black = COLOR([0,0,0])(partial_cover)
var partial_cover_white = COLOR([1,1,1])(R([0,1])(PI/6)(partial_cover))
var partial_cover_bw = STRUCT([partial_cover_white,partial_cover_black])

var cover = T([2])([1.26])(STRUCT([rotate_form(partial_cover_bw,[0,1],PI/3,5),knob]))

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



var sugarbowl = STRUCT([sugarbowl_form,cover,sb])