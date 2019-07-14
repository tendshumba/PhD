
M:=ZeroMatrix(R,Dimension(V),Dimension(V));
time
for i in [1..Dimension(V)] do 
	v:=BB[i];
	pd:=scaledProd(a,v);
	C:=Matrix(R,[[pd[k]:k in [1..Dimension(V)]]])*MatInv;
	for j in [1..Dimension(V)] do 
		M[i][j]:=C[1][j];
	end for;
end for;
