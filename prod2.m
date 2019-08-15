prod2:=function(u,v)
if not ( u in V) or  not( v in V) then
  print "error, vector must be in V";
else
  return &+[u[Floor(i/1782)+1]*v[i-Floor(i/1782)*1782]*prods1[i]:i in [1..#prods1]];
end if;
end function;
