
 Mat:=ZeroMatrix(R,Dimension(V),Dimension(V));
for i in [1..Dimension(V)] do
	 for j in [1..Dimension(V)] do
 		Mat[i][j]:=BB[i][j];
 	end for;
 end for;
 MatInv:=Mat^-1;
