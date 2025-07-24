### Fusion Law Checks
### Only coded for Jordan and Monster fusion laws.

### Jordan Fusion Law Check

JordanCheck:= function(V,a,eta)
local Ad, n, f1, f0, feta, f10, A, B, V1, V0, Veta, V10, i, j, ret,prod;

ret:=[];
Ad:=V[3]*a;
n:=Dimension(V[1]);

A:=V[1];
B:=Basis(V[1]);

f1:=(B*Ad-eta*B)*Ad;
f0:=(B*Ad-B)*Ad-eta*(B*Ad-B);
feta:=(B*Ad-B)*Ad;
f10:=B*Ad-eta*B;

V1:=Subspace(A,f1);
V0:=Subspace(A, f0);
Veta:=Subspace(A, feta);
V10:=Subspace(A, f10);

prod:=V[4];

if not prod(a,a)-a in Subspace(A,[]) then 
Add(ret,[1]);
fi;

if not Dimension(V1)+Dimension(V0)+Dimension(Veta)=Dimension(A) then
Add(ret,[2]);
fi;

if not Dimension(V1)=1 then
   Add(ret,[4]);
fi;

for i in f0 do
for j in f0 do
if not  prod(i,j) in V0 then
Add(ret,[0,0,i,j]);
fi;
od;od;

for i in f0 do
for j in feta do
if not  prod(i,j) in Veta then
Add(ret,[0,eta,i,j]);
fi;
od;od;


for i in feta do
for j in feta do
if not  prod(i,j) in V10 then
Add(ret,[eta,eta,i,j]);
fi;
od;od;
if ret=[] then Print("Primitive Jordan axis");
else Print("Not Primitive Jordan axis");
if ret<>0 then return ret;
fi;
fi;
end;

### Monster Fusion Law Check

MonsterCheck:= function(V,a,al,bt)
local Ad, n, f1, f0, fal, fbt, f10, f10al, A, B, V1, V0, Val, Vbt, V10, V10al, i, j, ret, prod;

ret:=[];
Ad:=V[3]*a;
n:=Dimension(V[1]);

A:=V[1];
B:=Basis(A);

f1:=((B*Ad-bt*B)*Ad)*Ad-al*(B*Ad-bt*B)*Ad;
f0:=((B*Ad-B)*Ad-bt*(B*Ad-B))*Ad-al*((B*Ad-B)*Ad-bt*(B*Ad-B));
fal:=((B*Ad-B)*Ad)*Ad-bt*((B*Ad-B)*Ad);
fbt:=((B*Ad-B)*Ad)*Ad-al*((B*Ad-B)*Ad);
f10:=(B*Ad-bt*B)*Ad-al*(B*Ad-bt*B);
f10al:=(B*Ad-bt*B);

V1:=Subspace(A,f1);
V0:=Subspace(A, f0);
Val:=Subspace(A, fal);
Vbt:=Subspace(A, fbt);
V10:=Subspace(A, f10);
V10al:=Subspace(A, f10al);

prod:=V[4];

if not prod(a,a)-a in Subspace(A,[]) then 
Add(ret,[1]);
fi;

if not Dimension(V1)+Dimension(V0)+Dimension(Val)+Dimension(Vbt)=Dimension(A) then
Add(ret,[2]);
fi;

if not Dimension(V1)=1 then
   Add(ret,[4]);
fi;

for i in f0 do
for j in f0 do
if not  prod(i,j) in V0 then
Add(ret,[0,0,i,j]);
fi;
od;od;

for i in f0 do
for j in fal do
if not  prod(i,j) in Val then
Add(ret,[0,al,i,j]);
fi;
od;od;

for i in f0 do
for j in fbt do
if not  prod(i,j) in Vbt then
Add(ret,[0,bt,i,j]);
fi;
od;od;

for i in fal do
for j in fal do
if not  prod(i,j) in V10 then
Add(ret,[al,al,i,j]);
fi;
od;od;

for i in fal do
for j in fbt do
if not  prod(i,j) in Vbt then
Add(ret,[al,bt,i,j]);
fi;
od;od;

for i in fbt do
for j in fbt do
if not  prod(i,j) in V10al then
Add(ret,[bt,bt,i,j]);
fi;
od;od;

if ret=[] then Print("Primitive Monster axis");
else Print("Not Primitive Monster axis  ");
if ret<>0 then return ret;
fi;
fi;
end;
