///////////////////////////////////////////////////////////////////////////////////////////////////////
// Algebra product 1 here,i.e.,pointwise multiplication and projection. We won't//sacle the product.
//
// ///////////////////////////////////////////////////////////////////////////////////////////////
prod1:=function(u,v)
        vec:=WW![u[i]*v[i]:i in [1..Degree(WW)]];/*algebra product is pointwise multiplication and projection */
return proj(vec);
end function;
///////////////////////////////////////////////////////////////////////////////
//Product 2 is presented here.
//
//////////////////////////////////////////////////////////////////////////////
times:=function(u,v)/*multiplication of standard basis vectors*/
  s_1:=Support(u);s_2:=Support(v);
S1:=SetToSequence(s_1);S2:=SetToSequence(s_2);
i1:=S1[1];i2:=S2[1];
if i1 eq i2 then 
return Zero(WW);
 else
   return u+v;
end if;
end function;
///////////////////////////////////////////////////////////////////////////////
//
//We present here the second algebra product here.Takes as input two vectors.
//Again we won't scale the function yet.
///////////////////////////////////////////////////////////////////////////////
prod2:=function(u,v)/*arguments are vectors u,v\in WW.*/
  //S_1:=Support(u);S_2:=Support(v);
  //S:=S_1 join S_2;
StBas:=[];/*we set up part of the standard basis for W (WW) involved in u&v.*/
for i in [1..Degree(WW)] do 
      w:=Zero(WW);
      w[i]:=1;
      Append(~StBas,w);
end for;
/*We rewrite u,v in terms of the standard basis now.*/
pd:=Zero(WW);
for i in [1..Degree(WW)] do 
      c:=u[i];
      for j in [1..Degree(WW)] do
	    d:=v[j];
            pd:=pd+(c*d)*times(StBas[i],StBas[j]);
      end for;
end for;
return Proj(pd);/*function returns the projection of the product to V=780.*/
end function;
/////////////////////////////////////////////////////////////////////////////
//In all likelihood, the product 1 and 2 are dependent.( based on evdience of wh//at been done so far. Conclusive test being made.
///////////////////////////////////////////////////////////////////////////////
OrbMat:=ZeroMatrix(R,1782,1782);
for i in [1..1782] do 
	for j in [1..1782] do 
		if i ne j then 
			if <i,j> in orbs[2] then
				OrbMat[i][j]:=1;
			end if;
		end if;
	end for;
end for;
//////////////////////////////////////////////////////////////////////////////
//we will set up the function times1 which is defined on the standard basis.
//a pair w_i,w_j multiplies to w_i+w_j if <i,j> is the orbit 1782*416 of Suz.2 
//acting on pairs
///////////////////////////////////////////////////////////////////////////////
times1:=function(u,v)/*multiplication of standard basis vectors*/
  s_1:=Support(u);s_2:=Support(v);
S1:=SetToSequence(s_1);S2:=SetToSequence(s_2);
i1:=S1[1];i2:=S2[1];
if not OrbMat[i1][i2] eq 1 then
return Zero(WW);
 else
   return u+v;
end if;
end function;

////////////////////////////////////////////////////////////////////////////////// we set up prod3 here. 
//
//////////////////////////////////////////////////////////////////////////////
prod3:=function(u,v)
	StBas:=[];/*we set up part of the standard basis for W (WW) involved in u&v.*/
for i in [1..Degree(WW)] do
      w:=Zero(WW);
      w[i]:=1;
      Append(~StBas,w);
end for;
/*We rewrite u,v in terms of the standard basis now.*/
pd:=Zero(WW);
for i in [1..Degree(WW)] do
      c:=u[i];
      for j in [1..Degree(WW)] do
            d:=v[j];
            pd:=pd+(c*d)*times1(StBas[i],StBas[j]);
      end for;
end for;
return Proj(pd);/*function returns the projection of the product to V=780.*/
end function;

times2:=function(u,v)/*this function takes arguments standard basis vectors.*/
  S1:=Support(u);S2:=Support(v);
S1:=SetToSequence(S1);S2:=SetToSequence(S2);
i:=S1[1];j:=S2[1];
if not OrbMat[i][j] eq 1 then /*check if <i,j>\in O.*/
   return Zero(WW);
 else
 x:={i,j};
KK:=Stabilizer(G,x);
XX:=GSet(G,{1..1782});
Orbs_KK:=Orbits(KK,XX);
pdct:=Zero(WW);
for l in Orbs_KK[2] do
	wl:=Zero(WW);wl[l]:=1; 
      pdct:=pdct+wl;
end for;
return pdct;
end if;
end function;
//////////////////////////////////////////////////////////////////////////////
//
//prod4 here.
/////////////////////////////////////////////////////////////////////////////
prod4:=function(u,v)
        StBas:=[];/*we set up part of the standard basis for W (WW) involved in u&v.*/
for i in [1..Degree(WW)] do
      w:=Zero(WW);
      w[i]:=1;
      Append(~StBas,w);
end for;
/*We rewrite u,v in terms of the standard basis now.*/
pd:=Zero(WW);
for i in [1..Degree(WW)] do
      c:=u[i];
      for j in [1..Degree(WW)] do
            d:=v[j];
            pd:=pd+(c*d)*times2(StBas[i],StBas[j]);
      end for;
end for;
return Proj(pd);/*function returns the projection of the product to V=780.*/
end function;
/////////////////////////////////////////////////////////////////////////////////////
//
//This should be faster than all the other point-wise multiplication followed by projection
//
///////////////////////////////////////////////////////////////////////////////////////////

Prod:=function(u,v)
vec:=WW!([u[i]*v[i]:i in [1..1782]]);
return vec*ProjV;
end function;
scaledProd:=function(u,v)
supp:=SetToSequence(Support(Prod(a,a)));
Pd:=Prod(u,v);
return (Pd)/(Prod(a,a)[supp[1]]);
end function;
