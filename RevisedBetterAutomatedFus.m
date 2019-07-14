

/*The grand aim here is that suppose we are given the ad_a matrix M and the algebra product.
can we find the fusion rules without human inteference?*/
/*Dependences:we will need ToW, to change the vectors of length rank M to length of vectors of V.*/
load "ToW.m";
/*Set Up Quat.*/
Quat:=ZeroMatrix(R,Degree(V),Degree(V));
for i in [1..Dimension(V)] do
	Quat[i]:=ToW(M[i]);
end for;
eigs:=Eigenvalues(Quat);
eigs:=SetToSequence(eigs);/*idea is to allow indexing. Can't index sets.*/
evals:=[x[1]:x in eigs|x[1] ne 0];/*extracting just the eigenvalues.*/
Evals_e:=[];/*initialising the set of eigenvalues whose coresponding spaces are even.*/
Evals_o:=[];/*same story as above, odd this time.*/
Even:=[];/*all the bases of even eigenspace will be placed here.*/
Odd:=[];
Dims_e:=[];Dims_o:=[];
for x in evals do 
  E:=Eigenspace(M,x);
  Bas:=[&+[(E.i)[j]*BB[j]:j in [1..Dimension(V)]]:i in [1..Dimension(E)]];
  if forall{v:v in Bas|v^p eq v} then /*everything in the space is even*/
   Append(~Evals_e,x);Append(~Dims_e,#Bas);Even:=Even cat Bas;/*book keeping exercises here.*/
   elif forall{v:v in Bas|v^p eq -1*v}  then /*everything is odd.*/ 
     Append(~Evals_o,x);Append(~Dims_o,#Bas);Odd:=Odd cat Bas;
   else
     E_e:=sub<WW|[(1/2)*(v+v^p):v in Bas]>; E_o:=sub<WW|[1/2*(v-v^p):v in Bas]>;
     Append(~Evals_e,x);Append(~Dims_e,Dimension(E_e));Append(~Evals_o,x);Append(~Dims_o,Dimension(E_o));
     Even:=Even cat [v:v in Basis(E_e)];Odd:=Odd cat [v:v in Basis(E_o)];
   end if;
 end for;
 /*At this stage we have two bases for the even and odd parts of V respectively.*/
VBas:=Even cat Odd; /*VV:=VectorSpaceWithBasis(VBas);*/
Evals:=Evals_e;
for x in Evals_o do 
if not x in Evals_e then 
Append(~Evals,x); printf "Eigenvalues %o is completely odd \n",x;
	else
		eig:=x+1;
		while eig/1 in Evals  or eig/1 in Evals_o do
			eig:=eig+1;
		end while;
			Append(~Evals,eig/1); printf "the %o-eigenvalue is using the dummy value %o for the odd part\n",x,eig;
		end if;
	end for;
Dims:=Dims_e cat Dims_o;
Inds:=[];
for i in [1..#Evals] do /*this piece sets up Inds, where each list in inds indicates positions occupied by each subspace.*/ 
  if i eq 1 then 
    Append(~Inds,[1..Dims[1]]);
  else
    Append(~Inds,[&+[Dims[c]:c in [1..i-1]]+1..&+[Dims[cc]:cc in [1..i]]]);
  end if;
end for;
//e1:=Zero(WW); e1[1]:=1;p1:=Proj(e1);
//C:=Coordinates(V,Proj(e1));
projs:=[];
time for i in [1..Degree(WW)] do 
	w:=Zero(V);
	w[i]:=1;
	Append(~projs,w*ProjV);
end for;

///////////////////////////////////////////////////////////////////
//New (better) projection function
///////////////////////////////////////////////////////////
proj:=function(u)
	return &+[u[i]*projs[i]: i in [1..Degree(u)]];
end function;

////////////////////////////////////////////////////////////////////////
//Prod is the new algebra product.
//
/////////////////////////////////////////////////////////////////////////
Prod:=function(u,v)
	return &+[u[i]*v[i]*projs[i]: i in [1..Degree(u)]];
end function;

////////////////////////////////////////////////////////////////////////
//scaledProd is the product above scaled so that a is an idempotent
//////////////////////////////////////////////////////////////////////
scaledProd:=function(u,v)
supp:=SetToSequence(Support(Prod(a,a)));
return (1/Prod(a,a)[supp[1]])*Prod(u,v);
end function; 
 //Mat:=Matrix(R,[Eltseq(v): v in VBas]);
//time sols:=Solution(Mat, projs);
/*Set Up TMInv*/
d:=Dimension(V);
TM:=ZeroMatrix(R,Dimension(V),d);
for j in [1..Dimension(V)] do
	vv:=Zero(VectorSpace(R,d));
	for i in [1..d] do
		vv[i]:=VBas[j][i];
	end for;
	TM[j]:=vv;
end for;
TMInv:=TM^-1;
//time sols:=[Matrix(R,[[projs[j][i]:i in [1..d]]])*TMInv: j in [1..#projs]];
//ProjSpaces:=[];/*here we get the projections of each standard basis vector to every eigenspace.*/ 
 //for i in [1..#Inds] do 
//	space:=[];
//	if #Inds[i] eq 1 then 
//		for j in [1..#sols] do 
//			Append(~space,sols[j][1][i]*VBas[i]);
//		end for;
//	else
//		for k in [1..#sols] do 
//			zz:=Zero(V);
//			for j in [Inds[i][1]..Inds[i][2]] do
//				zz:=zz+sols[k][1][j]*VBas[j];
//			end for;
//			Append(~space,zz);
//		end for;
//	end if;
//	Append(~ProjSpaces,space);
 //end for;
/*next order of business is to find out how eigenvalues multiply.*/
time 
for s in [1..#Evals] do 
  for t in [1..#Evals] do
    if s le t then 
      prods:=[];
      if Dims[s] le Dims[t] then 
	for i in [1..#Inds[s]] do 
	  for j in [1..#Inds[t]] do 
	    if i le j then 
	      u:=VBas[Inds[s][i]];v:=VBas[Inds[t][j]];
	      pd:=scaledProd(u,v);
		c:=Matrix(R,[[pd[u]:u in [1..Dimension(V)]]])
	      for l in [1..#Inds] do 
		if &+[c[1][m]*VBas[l][m]:m in Inds[l]] ne Zero(V) then
		 if not Evals[l] in prods then
		  Append(~prods,Evals[l]);
		end if;
	      end if;
	      end for;
	    end if;
	  end for;
	end for;
	printf "E%o times E%o is %o\n",s,t,prods;
      else /*dim Es greater that dim Et*/
	prods:=[];
	for i in [1..#Inds[s]] do
	  for j in [1..#Inds[t]] do
	    if i ge j then 
	      u:=VBas[Inds[s][i]];v:=VBas[Inds[t][j]];
	      pd:=scaledProd(u,v);
		c:=Matrix(R,[[pd[u]]:u in [1..Dimension(V)]]);
	      for l in [1..#Inds] do 
		if &+[c[1][l]*VBas[l][m]:m in Inds[l]] ne Zero(V) then 
		  if not Evals[l] in prods then 
		    Append(~prods,Evals[l]/1);
		  end if;
		end if;
	      end for;
	    end if;
	  end for;
	end for;
	printf "E%o times E%o is %o\n",s,t,prods;
      end if;
    end if;
  end for;
end for;


      



