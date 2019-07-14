ToW:=function(u)
 u_1:=Eltseq(u);
 for i in [1..(Dimension(WW)-#u_1)] do 
 	Append(~u_1,0);
 end for;
 return WW!u_1;
end function;

ToWSeq:=function(l)
 for i in [1..Dimension(WW)-#l] do
 	Append(~l,0);
 end for;
 return l;
end function;
