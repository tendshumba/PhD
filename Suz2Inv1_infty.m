SetLogFile("Suz2Inv1.txt");
SetSeed(1);
G:=PrimitiveGroups(1782)[2];
CompositionFactors(G);
CClas:=ConjugacyClasses(G);
invs:=[cc: cc in CClas|cc[1] eq 2];/* filtering to get only the classes of involutions*/
#invs;

X:=GSet(G,{t:t in CartesianPower(GSet(G),2)});
orbs:=Orbits(G,X);
lorbs:=[#l:l in orbs];
R:=Rationals();
W:=PermutationModule(G,R);
A:=ZeroMatrix(R,lorbs[1],lorbs[1]);/*setting up the adjacency matrix of the Suz graph*/
for x in orbs[2] do
        A[x[1]][x[2]]:=Identity(R);
 end for;
ee:=Eigenvalues(A);
ee; /*verifying that we have the right thing here.*/
BB:=Basis(Eigenspace(A,20)) cat Basis(Eigenspace(A,416)) cat Basis(Eigenspace(A,-16));
ee1:=SetToSequence(ee);
ind:=[v:v in ee1|v[2] eq 780];
ind;
inde:=Index(ee1,ind[1]);/*establishing the index of the 780 dim module in the list of eigenvalues*/
ranges:=[];
for i in [1..#ee1] do
    if i eq 1 then
                 range:=[1,ee1[i][2]];
                 Append(~ranges,range);
         else
                 range:=[ranges[i-1][2]+1,ranges[i-1][2]+ee1[i][2]];
                 Append(~ranges,range);

         end if;
 end for;
ranges;/*finding the positions each eigenspace occupies in BB*/
ranges[inde];
/*V:=VectorSpaceWithBasis([BB[i]: i in [ranges[inde][1]..ranges[inde][2]]]); Space on which we want to define an algebra product*/
WW:=VectorSpaceWithBasis(BB);
V:=sub<WW|[WW!BB[i]:i in [ranges[inde][1]..ranges[inde][2]]]>;
Dimension(V);

 p:=invs[1][3];/*choosing an involution to use.*/
 K:=Centraliser(G,p);
CompositionFactors(K);
 printf "The order of the centraliser is %o\n",Order(K);
 oo:=Orbits(K);
 invbas:=[];
 for o in oo do
    vec:=Zero(WW);
    for i in o do
                vec[i]:=Identity(R);
        end for;
     Append(~invbas,vec);
 end for;/*basis for the fixed space of K*/
Invspace:=sub<WW|invbas>;
Intspace:=V meet Invspace;
                                                              
Dimension(Intspace);
if Dimension(Intspace) ne 0 then
 a:=Basis(Intspace)[1];
end if;
//end for;
////////////////////////////////////////////////////////////////////////////////////////////////
// The projection function here.
//
// ///////////////////////////////////////////////////////////////////////////////////////
Proj:=function(u)
  proj:=Zero(WW);
   for i in [ranges[inde][1]..ranges[inde][2]] do
   c:=Coordinates(WW,u)[i];
        proj:=proj+c*BB[i];
    end for;
 return proj;
 end function;
Proj(a) eq a; /*checking if it works OK.*/
////////////////////////////////////////////////////////////////////////////////////////////////////////
                                                             
prod:=function(u,v)
        vec:=WW![u[i]*v[i]:i in [1..Degree(WW)]];/*algebra product is pointwise multiplication and projection */
return Proj(vec);
end function;
scaledprod:=function(u,v)
       Sup:=SetToSequence(Support(a));
i_1:=Sup[1];
epsi:=prod(a,a)[i_1]/a[i_1];
	return (1/(epsi))*prod(u,v);// this is to make a an idempotent.
end function;
Proj:=function(u)
  proj:=Zero(WW);
C:=Coordinates(WW,u);
   for i in [ranges[inde][1]..ranges[inde][2]] do
   c:=C[i];
        proj:=proj+c*BB[i];
    end for;
 return proj;
 end function;
Proj(a) eq a; /*checking if it works OK.*/
load "BetterProj.m";
load "Suz2Prods.m";
time scaledprod(a,a) eq a;
