function lar_to_obj(model) {
 v = model[0]
 fv = model[1]
 s = "\n"
 for (i=0; i<v.length;i++){
 	s += "v "+v[i][0] +" " + v[i][1] +" " +0 + "\n"
 }
 s += "\n"
 for (i=0; i<fv.length;i++){
    s += "f "
 	for (j=0; j<fv[i].length;j++)
	 	s+= fv[i][j] +" "
 s += "\n"
 }
 return s;
}

//Lascio un codice d'esempio che mostra il risultato della funzione
V = [[0,6],[0,0],[3,0],[6,0],[0,3],[3,3],[6,3],[6,6],[3,6]]
FV = [[5,6,7,8],[0,5,8],[0,4,5],[1,2,4,5],[2,3,5,6],
//facce estene
[0,8,7],[3,6,7],[1,2,3],[0,4,1]]

l = lar_to_obj([V,FV])


