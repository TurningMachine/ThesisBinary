###OtherFunctions
###Consists of Frobenius Form, Radical, Ideals, Quotients and Miyamoto invariant checks.

### Frobenius Form Check

FFCheck:=function(V)
local frob, A, B, n, prod, i, j, k, ret;

ret:=[];
frob:=V[5];
A:=V[1];
B:=Basis(A);
n:=Dimension(A);
prod:=V[4];

for i in [1..n] do
for j in [1..n] do
for k in [1..n] do

if prod(B[i],B[j])*frob*B[k]<>B[i]*frob*prod(B[j],B[k]) then Add(ret, [B[i],B[j],B[k]]);
fi;

od;
od;
od;

if ret=[] then Print("Frobenius Form");
else Print("Not Frobenius Form  ");
if ret<>0 then return ret;
fi;
fi;
end;

### Radical Related

FFDeterminant:=function(V,gens)
local frob, n, i, z, det;

frob:=V[5];
n:=Size(gens);

z:= 1;
for i in [1..n] do
z:= (gens[i])^0*z;
od;

det:=(Determinant(z*frob));

return det;
end;

FFDetInfo:=function()
Print("If Algebra has any indeterminants (al, bt, dt etc.), please put them in [] for the second component of the function");
end;


Rad:=function(V)
local list, i, n, ret, FF, frob;
ret:=[];
FF:=V[2];
frob:=V[5];

if FFDeterminant(V,[1])<>0 then Print("Simple Algebra ");
else list:=Eigenvectors(FF,frob);
n:=Size(list);
for i in [1..n] do
if frob*list[i]=Zero(V[1]) then Add(ret,list[i]);
fi;
od;
fi;
return ret;
end;

RadInfo:=function()
Print("This function only works with no indeterminants!");
end;

###Quotient Related

IdealCheck:=function(V,I)
local A, B, BI, n, m, prod, ret, i, Ad, j;
A:=V[1];
B:=Basis(A);
BI:=Basis(I);
n:=Size(B);
m:=Size(BI);
prod:=V[4];
ret:=[];

for i in [1..n] do
for j in [1..m] do
if not  prod(B[i],BI[j]) in I then
Add(ret,[i,j]);
fi;
od;od;

if ret=[] then Print("Ideal");
else Print("Not Ideal");
if ret<>0 then return ret;
fi;
fi;
end;

QuotientAxial:=function(V,I)
local A, F, Q, B, h, n, mult, ret, x, prod, Qprod, frob, new, j, k;
A:=V[1];
F:=V[2];
Q:=A/I;
B:=Basis(Q);
n:=Size(B);
mult:=V[3];

h:= NaturalHomomorphismBySubspace(A,I);
x:=NullMat(1,n);
for i in [1..n] do
x[i]:=PreImagesRepresentative( h, B[i] );
od;

ret:=NullMat(1,n);

for i in [1..n] do
ret[i]:=NullMat(n,n);
for j in [1..n] do
ret[i][j]:=Image(h,x[i]*(x[j]*mult));
od;od;

prod:=V[4];
Qprod := function(u,v)

local x,y;
x:=PreImagesRepresentative( h, u );
y:=PreImagesRepresentative( h, v );
return Image(h, prod(x,y));
end;


frob:=V[5];

new:=NullMat(n,n);
for i in [1..n] do
for j in [1..n] do
new[i][j]:=(x[i]*frob)*x[j];
od;od;


return [Image(h,A),F,ret,Qprod,new];;
end;

###Invariant under Miyamoto involutions check for Monster Fusion Law

MiyamotoInvariant:= function(V,bt,S)
local Ad, n, f10al, A, B, A10al, s, I, a, i,j, mas, Bas,b, new, prod;

s:=Size(S);
I:=V[1];
A:=V[1];
n:=Dimension(I);
B:=Basis(A);
for i in [1..s] do
a:=S[i];
Ad:=V[3]*a;
f10al:=(Ad*B-bt*B);
A10al:=Subspace(A, f10al);
I:=Intersection(A10al,I);
od;
return [I,V[2],V[3],V[4],V[5]];;
end;
