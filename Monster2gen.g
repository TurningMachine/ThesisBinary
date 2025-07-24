### Monster 2gen
### Here we present the 2-generated axial algebras of Monster type that are used in the thesis (or the ones needed)

###Highwater quotients

###H4
HW4:=function(F)
local mas, prod, a0, a1, a2, a3, h, c, A, frob;


## order of the basis: a0, a1, a2, a3, h, c

mas :=[
[[1,0,0,0,0,0], (1/4)*[3,3,1,1,-3,-4], (1/2)*[1,0,1,0,3,0], (1/4)*[3,1,1,3,-3,-4], (1/2)*[-1,0,1,0,3,0], (1/2)*[3,-1,-1,-1,0,2]],
[(1/4)*[3,3,1,1,-3,-4], [0,1,0,0,0,0], (1/4)*[1,3,3,1,-3,-4], (1/2)*[0,1,0,1,3,0], (1/2)*[0,-1,0,1,3,0], (1/2)*[-1,3,-1,-1,0,2]], 
[(1/2)*[1,0,1,0,3,0], (1/4)*[1,3,3,1,-3,-4], [0,0,1,0,0,0], (1/4)*[1,1,3,3,-3,-4], (1/2)*[1,0,-1,0,3,0], (1/2)*[-1,-1,3,-1,0,2]],
[(1/4)*[3,1,1,3,-3,-4], (1/2)*[0,1,0,1,3,0], (1/4)*[1,1,3,3,-3,-4], [0,0,0,1,0,0], (1/2)*[0,1,0,-1,3,0], (1/2)*[-1,-1,-1,3,0,2]],
[(1/2)*[-1,0,1,0,3,0], (1/2)*[0,-1,0,1,3,0], (1/2)*[1,0,-1,0,3,0], (1/2)*[0,1,0,-1,3,0], [0,0,0,0,1,0], [0,0,0,0,0,0]],
[(1/2)*[3,-1,-1,-1,0,2], (1/2)*[-1,3,-1,-1,0,2], (1/2)*[-1,-1,3,-1,0,2], (1/2)*[-1,-1,-1,3,0,2], [0,0,0,0,0,0], [0,0,0,0,0,1]]
];

prod := function(u,v)
local i,j,k, ans;
ans:=[0,0,0,0,0,0];
for i in [1..6] do
for j in [1..6] do
for k in [1..6] do
ans[k] := ans[k] + u[i]*v[j]*mas[i][j][k];
od;
od;
od;
return ans;
end;

a0:=[1,0,0,0,0,0];;
a1:=[0,1,0,0,0,0];;
a2:=[0,0,1,0,0,0];;
a3:=[0,0,0,1,0,0];;
h:=[0,0,0,0,1,0];;
c:=[0,0,0,0,0,1];;

A:=VectorSpace(F,[a0,a1,a2,a3,h,c]);;

frob:=[[1,1,1,1,0,1],[1,1,1,1,0,1],[1,1,1,1,0,1],[1,1,1,1,0,1],[0,0,0,0,0,0],[1,1,1,1,0,1]];

return [A,F,mas,prod,frob];
end;

HW4Details:=function();
Print("1:Vector Space, 2:Field, 3:Mult Matrix, 4:Mult Function, 5:Frobenius Matrix, Basis is a0,a1,a2,a3,h,c"); 
end;

###Algebras with fixed axet

###4A(1/4,beta)

4A:=function(bt,FF)
local F, d, mas, prod, a0, a1, a2, a3, q, A, frob;

if bt in FF then F:=FF;
else F:=FunctionField(FF,["bt"]);
fi;

## order of the basis: a0, a1, a2, a3, q

d:=bt^0;

mas :=d * [
[[1,0,0,0,0], [(1+4*bt)/8, (1+4*bt)/8, (1-4*bt)/8,(1-4*bt)/8,1], [0,0,0,0,0], [(1+4*bt)/8, (1-4*bt)/8, (1-4*bt)/8,(1+4*bt)/8,1], [(2*bt-1)/8,0,0,0,0]],
[[(1+4*bt)/8, (1+4*bt)/8, (1-4*bt)/8,(1-4*bt)/8,1], [0,1,0,0,0], [(1-4*bt)/8, (1+4*bt)/8, (1+4*bt)/8,(1-4*bt)/8,1], [0,0,0,0,0], [0,(2*bt-1)/8,0,0,0]], 
[[0,0,0,0,0], [(1-4*bt)/8, (1+4*bt)/8, (1+4*bt)/8,(1-4*bt)/8,1], [0,0,1,0,0], [(1-4*bt)/8, (1-4*bt)/8, (1+4*bt)/8,(1+4*bt)/8,1], [0,0,(2*bt-1)/8,0,0]],
[[(1+4*bt)/8, (1-4*bt)/8, (1-4*bt)/8,(1+4*bt)/8,1], [0,0,0,0,0], [(1-4*bt)/8, (1-4*bt)/8, (1+4*bt)/8,(1+4*bt)/8,1], [0,0,0,1,0], [0,0,0,(2*bt-1)/8,0]],
[[(2*bt-1)/8,0,0,0,0], [0,(2*bt-1)/8,0,0,0], [0,0,(2*bt-1)/8,0,0], [0,0,0,(2*bt-1)/8,0], [0,0,0,0,(2*bt-1)/8]]
];

prod := function(u,v)
local i,j,k, ans;
ans:=d * [0,0,0,0,0];
for i in [1..5] do
for j in [1..5] do
for k in [1..5] do
ans[k] := ans[k] + u[i]*v[j]*mas[i][j][k];
od;
od;
od;
return ans;
end;

a0:=d*[1,0,0,0,0];;
a1:=d*[0,1,0,0,0];;
a2:=d*[0,0,1,0,0];;
a3:=d*[0,0,0,1,0];;
q:=d*[0,0,0,0,1];;

A:=VectorSpace(F,[a0,a1,a2,a3,q]);;

frob:=[[1,bt,0,bt,(2*bt-1)/8],[bt,1,bt,0,(2*bt-1)/8],[0,bt,1,bt,(2*bt-1)/8],[bt,0,bt,1,(2*bt-1)/8],((2*bt-1)/8)*[1,1,1,1,(2*bt-1)/2]];

return [A,F,mas,prod,frob];;
end;

4ADetails:=function();
Print("1:Vector Space, 2:Field, 3:Mult Matrix, 4:Mult Function, 5:Frobenius Matrix, Basis is a0,a1,a2,a3,q"); 
end;

###4A(1/4,1/2)x

4AX:=function(F)
local mas, prod, a0, a1, a2, a3, A, frob;

## order of the basis: a0, a1, a2, a3,

mas :=[
[[1,0,0,0], [3/8, 3/8, -1/8, -1/8], [0,0,0,0], [3/8, -1/8, -1/8, 3/8]],
[[3/8, 3/8, -1/8, -1/8], [0,1,0,0], [-1/8, 3/8, 3/8,-1/8], [0,0,0,0]], 
[[0,0,0,0], [-1/8, 3/8, 3/8,-1/8], [0,0,1,0], [-1/8, -1/8, 3/8,3/8]],
[[3/8, -1/8, -1/8, 3/8], [0,0,0,0], [-1/8, -1/8, 3/8,3/8], [0,0,0,1]],
];

prod := function(u,v)
local i,j,k, ans;
ans:=[0,0,0,0];
for i in [1..4] do
for j in [1..4] do
for k in [1..4] do
ans[k] := ans[k] + u[i]*v[j]*mas[i][j][k];
od;
od;
od;
return ans;
end;

a0:=[1,0,0,0];;
a1:=[0,1,0,0];;
a2:=[0,0,1,0];;
a3:=[0,0,0,1];;

A:=VectorSpace(F,[a0,a1,a2,a3]);;

frob:=[[1,1/2,0,1/2],[1/2,1,1/2,0],[0,1/2,1,1/2,0],[1/2,0,1/2,1,]];

return [A,F,mas,prod,frob];;
end;

4AXDetails:=function();
Print("1:Vector Space, 2:Field, 3:Mult Matrix, 4:Mult Function, 5:Frobenius Matrix, Basis is a0,a1,a2,a3"); 
end;

###4B(al,al^2/2)

4B:=function(al,FF)
local F, bt, z, mas, prod, a0, a1, a2, a3, d, A, frob;

if al in [1,0,2] then Print("Not Defined!"); 
else

if al in FF then F:=FF;
else F:=FunctionField(FF,["al"]);
fi;
bt:=al^2/2;

## order of the basis: a0, a1, a2, a3, d

z:=al^0;

mas :=z * [
[[1,0,0,0,0], (bt/2)*[1,1,-1,-1,1], (al/2)*[1,0,1,0,-1], (bt/2)*[1,-1,-1,1,1], (al/2)*[1,0,-1,0,1]],
[(bt/2)*[1,1,-1,-1,1], [0,1,0,0,0], (bt/2)*[-1,1,1,-1,1], (al/2)*[0,1,0,1,-1], (al/2)*[0,1,0,-1,1]], 
[(al/2)*[1,0,1,0,-1], (bt/2)*[-1,1,1,-1,1], [0,0,1,0,0], (bt/2)*[-1,-1,1,1,1], (al/2)*[-1,0,1,0,1]],
[(bt/2)*[1,-1,-1,1,1], (al/2)*[0,1,0,1,-1], (bt/2)*[-1,-1,1,1,1], [0,0,0,1,0], (al/2)*[0,-1,0,1,1]],
[(al/2)*[1,0,-1,0,1], (al/2)*[0,1,0,-1,1], (al/2)*[-1,0,1,0,1], (al/2)*[0,-1,0,1,1], [0,0,0,0,1]]
];

prod := function(u,v)
local i,j,k, ans;
ans:=z * [0,0,0,0,0];
for i in [1..5] do
for j in [1..5] do
for k in [1..5] do
ans[k] := ans[k] + u[i]*v[j]*mas[i][j][k];
od;
od;
od;
return ans;
end;

a0:=z*[1,0,0,0,0];;
a1:=z*[0,1,0,0,0];;
a2:=z*[0,0,1,0,0];;
a3:=z*[0,0,0,1,0];;
d:=z*[0,0,0,0,1];;

A:=VectorSpace(F,[a0,a1,a2,a3,d]);;

frob:=[[1,bt/2,al/2,bt/2,al/2],[bt/2,1,bt/2,al/2,al/2],[al/2,bt/2,1,bt/2,al/2],[bt/2,al/2,bt/2,1,al/2],[al/2,al/2,al/2,al/2,1]];

return [A,F,mas,prod,frob];;
fi;
end;

4BDetails:=function();
Print("1:Vector Space, 2:Field, 3:Mult Matrix, 4:Mult Function, 5:Frobenius Matrix, Basis is a0,a1,a2,a3,d"); 
end;

###4B(-1,1/2)x

4BX:=function(F)
local mas, prod, a0, a1, a2, a3, A;

## order of the basis: a0, a1, a2, a3

mas :=[
[[1,0,0,0], (1/8)*[1,1,-3,-3], (-1/4)*[3,1,3,1], (1/8)*[1,-3,-3,1]],
[(1/8)*[1,1,-3,-3], [0,1,0,0], (1/8)*[-3,1,1,-3], (-1/4)*[1,3,1,3]], 
[(-1/4)*[3,1,3,1], (1/8)*[-3,1,1,-3], [0,0,1,0], (1/8)*[-3,-3,1,1]],
[(1/8)*[1,-3,-3,1], (-1/4)*[1,3,1,3], (1/8)*[-3,-3,1,1], [0,0,0,1]],
];

prod := function(u,v)
local i,j,k, ans;
ans:=[0,0,0,0];
for i in [1..4] do
for j in [1..4] do
for k in [1..4] do
ans[k] := ans[k] + u[i]*v[j]*mas[i][j][k];
od;
od;
od;
return ans;
end;

a0:=[1,0,0,0];;
a1:=[0,1,0,0];;
a2:=[0,0,1,0];;
a3:=[0,0,0,1];;

A:=VectorSpace(F,[a0,a1,a2,a3]);;

frob:=[[1,1/4,-1/2,1/4],[1/4,1,1/4,-1/2],[-1/2,1/4,1,1/4],[1/4,-1/2,1/4,1]];

return [A,F,mas,prod,frob];;
end;

4BXDetails:=function();
Print("1:Vector Space, 2:Field, 3:Mult Matrix, 4:Mult Function, 5:Frobenius Matrix, Basis is a0,a1,a2,a3"); 
end;

###4J(2bt,bt)

4J:=function(bt,FF)
local F, d, mas, prod, a0, a1, a2, a3, w, A, frob;

if bt in FF then F:=FF;
else F:=FunctionField(FF,["bt"]);
fi;

## order of the basis: a0, a1, a2, a3, w

d:=bt^0;

mas :=d * [
[[1,0,0,0,0], (bt/2)*[2,2,0,0,-1], [0,0,0,0,0], (bt/2)*[2,0,0,2,-1], bt*[2,-1,0,-1,1]],
[(bt/2)*[2,2,0,0,-1], [0,1,0,0,0], (bt/2)*[0,2,2,0,-1], [0,0,0,0,0], bt*[-1,2,-1,0,1]], 
[[0,0,0,0,0], (bt/2)*[0,2,2,0,-1], [0,0,1,0,0], (bt/2)*[0,0,2,2,-1], bt*[0,-1,2,-1,1]],
[(bt/2)*[2,0,0,2,-1], [0,0,0,0,0], (bt/2)*[0,0,2,2,-1], [0,0,0,1,0], bt*[-1,0,-1,2,1]],
[bt*[2,-1,0,-1,1], bt*[-1,2,-1,0,1], bt*[0,-1,2,-1,1], bt*[-1,0,-1,2,1], [0,0,0,0,1]]
];

prod := function(u,v)
local i,j,k, ans;
ans:=d * [0,0,0,0,0];
for i in [1..5] do
for j in [1..5] do
for k in [1..5] do
ans[k] := ans[k] + u[i]*v[j]*mas[i][j][k];
od;
od;
od;
return ans;
end;

a0:=d*[1,0,0,0,0];;
a1:=d*[0,1,0,0,0];;
a2:=d*[0,0,1,0,0];;
a3:=d*[0,0,0,1,0];;
w:=d*[0,0,0,0,1];;

A:=VectorSpace(F,[a0,a1,a2,a3,w]);;

frob:=[[1,bt,0,bt,2*bt],[bt,1,bt,0,2*bt],[0,bt,1,bt,2*bt],[bt,0,bt,1,2*bt],2*[bt,bt,bt,bt,1]];

return [A,F,mas,prod,frob];;
end;

4JDetails:=function();
Print("1:Vector Space, 2:Field, 3:Mult Matrix, 4:Mult Function, 5:Frobenius Matrix, Basis is a0,a1,a2,a3,w"); 
end;

###4J(-1/2,-1/4)x

4JX:=function(F)
local mas, prod, a0, a1, a2, a3, A, frob;

## order of the basis: a0, a1, a2, a3,

mas :=[
[[1,0,0,0], (-1/8)*[3,3,1,1], [0,0,0,0], (-1/8)*[3,1,1,3]],
[(-1/8)*[3,3,1,1], [0,1,0,0], (-1/8)*[1,3,3,1], [0,0,0,0]], 
[[0,0,0,0], (-1/8)*[1,3,3,1], [0,0,1,0], (-1/8)*[1,1,3,3]],
[(-1/8)*[3,1,1,3], [0,0,0,0], (-1/8)*[1,1,3,3], [0,0,0,1]],
];

prod := function(u,v)
local i,j,k, ans;
ans:=[0,0,0,0];
for i in [1..4] do
for j in [1..4] do
for k in [1..4] do
ans[k] := ans[k] + u[i]*v[j]*mas[i][j][k];
od;
od;
od;
return ans;
end;

a0:=[1,0,0,0];;
a1:=[0,1,0,0];;
a2:=[0,0,1,0];;
a3:=[0,0,0,1];;

A:=VectorSpace(F,[a0,a1,a2,a3]);;

frob:=[[1,-1/4,0,-1/4],[-1/4,1,-1/4,0],[0,-1/4,1,-1/4],[-1/4,0,-1/4,1,]];

return [A,F,mas,prod,frob];;
end;

4JXDetails:=function();
Print("1:Vector Space, 2:Field, 3:Mult Matrix, 4:Mult Function, 5:Frobenius Matrix, Basis is a0,a1,a2,a3"); 
end;

###4Y(al,1-al^2/2)

4Y:=function(al,FF)
local F, bt, z, mas, prod, a0, a1, a2, a3, e, A, gm, frob;

if al in FF then F:=FF;
else F:=FunctionField(FF,["al"]);
fi;
bt:=(1-al^2)/2;

## order of the basis: a0, a1, a2, a3, e

z:=al^0;

mas :=z * [
[[1,0,0,0,0], (1/2)*[bt,bt,-bt,-bt,(al+1)^2/2], (1/2)*[al-1,0,al-1,0,al+1], (1/2)*[bt,-bt,-bt,bt,(al+1)^2/2], (1/2)*[1-al,0,al-1,0,al+1]],
[(1/2)*[bt,bt,-bt,-bt,(al+1)^2/2], [0,1,0,0,0], (1/2)*[-bt,bt,bt,-bt,(al+1)^2/2], (1/2)*[0,al-1,0,al-1,al+1], (1/2)*[0,1-al,0,al-1,al+1]], 
[(1/2)*[al-1,0,al-1,0,al+1], (1/2)*[-bt,bt,bt,-bt,(al+1)^2/2], [0,0,1,0,0], (1/2)*[-bt,-bt,bt,bt,(al+1)^2/2], (1/2)*[al-1,0,1-al,0,al+1]],
[(1/2)*[bt,-bt,-bt,bt,(al+1)^2/2], (1/2)*[0,al-1,0,al-1,al+1], (1/2)*[-bt,-bt,bt,bt,(al+1)^2/2], [0,0,0,1,0], (1/2)*[0,al-1,0,1-al,al+1]],
[(1/2)*[1-al,0,al-1,0,al+1], (1/2)*[0,1-al,0,al-1,al+1], (1/2)*[al-1,0,1-al,0,al+1], (1/2)*[0,al-1,0,1-al,al+1], [0,0,0,0,1]]
];

prod := function(u,v)
local i,j,k, ans;
ans:=z * [0,0,0,0,0];
for i in [1..5] do
for j in [1..5] do
for k in [1..5] do
ans[k] := ans[k] + u[i]*v[j]*mas[i][j][k];
od;
od;
od;
return ans;
end;

a0:=z*[1,0,0,0,0];;
a1:=z*[0,1,0,0,0];;
a2:=z*[0,0,1,0,0];;
a3:=z*[0,0,0,1,0];;
e:=z*[0,0,0,0,1];;

A:=VectorSpace(F,[a0,a1,a2,a3,e]);;

gm:=(2-al)*(al+1)/4;

frob:=[[1,gm,al/2,gm,1-al/2],[gm,1,gm,al/2,1-al/2],[al/2,gm,1,gm,1-al/2],[gm,al/2,gm,1,1-al/2],((2-al)/2)*[1,1,1,1,2/(al+1)]];

return [A,F,mas,prod,frob];
end;

4YDetails:=function();
Print("1:Vector Space, 2:Field, 3:Mult Matrix, 4:Mult Function, 5:Frobenius Matrix, Basis is a0,a1,a2,a3,e"); 
end;

###4Y(1/2,bt)

4Yh:=function(bt,FF)
local F, z, mas, prod, a0, a1, a2, a3, f, A, frob;

if bt in FF then F:=FF;
else F:=FunctionField(FF,["bt"]);
fi;

## order of the basis: a0, a1, a2, a3, f

z:=bt^0;

mas :=z * [
[[1,0,0,0,0], (bt/2)*[1,1,-1,-1,8*bt], ((1-4*bt)/2)*[1,0,1,0,-8*bt], (bt/2)*[1,-1,-1,1,8*bt], (1/4)*[1,0,-1,0,8*bt]],
[(bt/2)*[1,1,-1,-1,8*bt], [0,1,0,0,0], (bt/2)*[-1,1,1,-1,8*bt], ((1-4*bt)/2)*[0,1,0,1,-8*bt], (1/4)*[0,1,0,-1,8*bt]], 
[((1-4*bt)/2)*[1,0,1,0,-8*bt], (bt/2)*[-1,1,1,-1,8*bt], [0,0,1,0,0], (bt/2)*[-1,-1,1,1,8*bt], (1/4)*[-1,0,1,0,8*bt]],
[(bt/2)*[1,-1,-1,1,8*bt], ((1-4*bt)/2)*[0,1,0,1,-8*bt], (bt/2)*[-1,-1,1,1,8*bt], [0,0,0,1,0], (1/4)*[0,-1,0,1,8*bt]],
[(1/4)*[1,0,-1,0,8*bt], (1/4)*[0,1,0,-1,8*bt], (1/4)*[-1,0,1,0,8*bt], (1/4)*[0,-1,0,1,8*bt], [0,0,0,0,1]]
];

prod := function(u,v)
local i,j,k, ans;
ans:=z * [0,0,0,0,0];
for i in [1..5] do
for j in [1..5] do
for k in [1..5] do
ans[k] := ans[k] + u[i]*v[j]*mas[i][j][k];
od;
od;
od;
return ans;
end;

a0:=z*[1,0,0,0,0];;
a1:=z*[0,1,0,0,0];;
a2:=z*[0,0,1,0,0];;
a3:=z*[0,0,0,1,0];;
f:=z*[0,0,0,0,1];;

A:=VectorSpace(F,[a0,a1,a2,a3,f]);;

frob:=[[1,4*bt^2,(4*bt-1)^2,4*bt^2,2*bt],[4*bt^2,1,4*bt^2,(4*bt-1)^2,2*bt],[(4*bt-1)^2,4*bt^2,1,4*bt^2,2*bt],[4*bt^2,(4*bt-1)^2,4*bt^2,1,2*bt],[2*bt,2*bt,2*bt,2*bt,1]];

return [A,F,mas,prod,frob];
end;

4YhDetails:=function();
Print("1:Vector Space, 2:Field, 3:Mult Matrix, 4:Mult Function, 5:Frobenius Matrix, Basis is a0,a1,a2,a3,f"); 
end;

###Algebras with variable axets

###IY3(al,1/2,1)

IY3:=function(al,FF)
local F, d, mas, prod, a0, a1, z, n, A, frob;

if al in FF then F:=FF;
else F:=FunctionField(FF,["al"]);
fi;

## order of the basis: a0, a1, z, n

d:=dt^0;

mas :=d * [
[[1,0,0,0],(1/2)*[1,1,2*al-1,2],[0,0,al,0],[0,0,0,0]],
[(1/2)*[1,1,2*al-1,2],[0,1,0,0],[0,0,al,0],[0,0,0,0]], 
[[0,0,al,0],[0,0,al,0], [0,0,0,0],[0,0,0,0]],
[[0,0,0,0],[0,0,0,0], [0,0,0,0],[0,0,0,0]]];

prod := function(u,v)
local i,j,k, ans;
ans:=d * [0,0,0,0];
for i in [1..4] do
for j in [1..4] do
for k in [1..4] do
ans[k] := ans[k] + u[i]*v[j]*mas[i][j][k];
od;
od;
od;
return ans;
end;

a0:=d*[1,0,0,0];;
a1:=d*[0,1,0,0];;
z:=d*[0,0,1,0];;
n:=d*[0,0,0,1];;

A:=VectorSpace(F,[a0,a1,z,n]);;

frob:=[[1,1,0,0],[1,1,0,0],[0,0,0,0],[0,0,0,0]];

return [A,F,mas,prod,frob];;
end;

IY3Details:=function();
Print("1:Vector Space, 2:Field, 3:Mult Matrix, 4:Mult Function, 5:Frobenius Matrix, Basis is a0,a1,z,n"); 
end;

###IY(al,1/2)

IY5:=function(al,FF)
local F, d, mas, prod, a0, a1, b0, b1, p, n, A, frob;

if al in FF then F:=FF;
else F:=FunctionField(FF,["al"]);
fi;

## order of the basis: a0, a1, b0, b1, p, n

d:=al^0;

mas :=d*[
[[1,0,0,0,0,0], (1/2)*[1,1,0,0,2,0], (1/2)*[0,0,1,0,-2,0], (1/4)*[0,0,0,2,-12,1], ((2*al-1)/8)*[-1,1,-1,0,4,0], [0,0,0,0,0,0]],
[(1/2)*[1,1,0,0,2,0], [0,1,0,0,0,0], (1/4)*[0,0,2,0,-12,1], (1/2)*[0,0,0,1,-2,0], ((2*al-1)/8)*[1,-1,0,-1,4,0], [0,0,0,0,0,0]], 
[(1/2)*[0,0,1,0,-2,0], (1/4)*[0,0,2,0,-12,1], [0,0,0,0,-2,0], [0,0,0,0,2,-1], ((2*al-1)/8)*[2,-2,1,-1,0,-1], [0,0,0,0,0,0]],
[(1/4)*[0,0,0,2,-12,1], (1/2)*[0,0,0,1,-2,0], [0,0,0,0,2,-1], [0,0,0,0,-2,0], ((2*al-1)/8)*[-2,2,-1,1,0,-1], [0,0,0,0,0,0]],
[((2*al-1)/8)*[-1,1,-1,0,4,0], ((2*al-1)/8)*[1,-1,0,-1,4,0], ((2*al-1)/8)*[2,-2,1,-1,0,-1], ((2*al-1)/8)*[-2,2,-1,1,0,-1], 
(((2*al-1)*(2*al-3))/32)*[0,0,0,0,0,1], [0,0,0,0,0,0]],
[[0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0]]
];


prod := function(u,v)
local i,j,k, ans;
ans:=d * [0,0,0,0,0,0];
for i in [1..6] do
for j in [1..6] do
for k in [1..6] do
ans[k] := ans[k] + u[i]*v[j]*mas[i][j][k];
od;
od;
od;
return ans;
end;

a0:=d*[1,0,0,0,0,0];;
a1:=d*[0,1,0,0,0,0];;
b0:=d*[0,0,1,0,0,0];;
b1:=d*[0,0,0,1,0,0];;
p:=d*[0,0,0,0,1,0];;
n:=d*[0,0,0,0,0,1];;

A:=VectorSpace(F,[a0,a1,b0,b1,p,n]);;

frob:=[[1,1,0,0,0,0],[1,1,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0]];

return [A,F,mas,prod,frob];;
end;

IY5Details:=function();
Print("1:Vector Space, 2:Field, 3:Mult Matrix, 4:Mult Function, 5:Frobenius Matrix, Basis is a0,a1,b0,b1,p,n"); 
end;

###Split Spin Factor Algebra

Split:=function(dt,al,FF)
local F, d, mas, prod, e0, e1, z1, z2, A, frob;

if al in FF then 
if dt in FF then F:=FF;
else F:=FunctionField(FF,["dt"]);
fi;
else if dt in FF then F:=FunctionField(FF,["al"]);
else F:=FunctionField(FF,["al","dt"]);
fi;
fi;

## order of the basis: e0, e1, z1, z2

d:=dt^0*al^0;

mas :=d * [
[[0,0,al*(2-al),(1-al)*(al+1)],[0,0,dt*al*(2-al),dt*(1-al)*(al+1)],[al,0,0,0],[1-al,0,0,0]],
[[0,0,dt*al*(2-al),dt*(1-al)*(al+1)],[0,0,al*(2-al),(1-al)*(al+1)],[0,al,0,0],[0,1-al,0,0]], 
[[al,0,0,0],[0,al,0,0], [0,0,1,0],[0,0,0,0]],
[[1-al,0,0,0],[0,1-al,0,0], [0,0,0,0],[0,0,0,1]]];

prod := function(u,v)
local i,j,k, ans;
ans:=d * [0,0,0,0];
for i in [1..4] do
for j in [1..4] do
for k in [1..4] do
ans[k] := ans[k] + u[i]*v[j]*mas[i][j][k];
od;
od;
od;
return ans;
end;

e0:=d*[1,0,0,0];;
e1:=d*[0,1,0,0];;
z1:=d*[0,0,1,0];;
z2:=d*[0,0,0,1];;

A:=VectorSpace(F,[e0,e1,z1,z2]);;

frob:=[(al+1)*(2-al)*[1,dt,0,0],(al+1)*(2-al)*[dt,1,0,0],[0,0,al+1,0],[0,0,0,2-al]];

return [A,F,mas,prod,frob];;
end;

SplitDetails:=function();
Print("1:Vector Space, 2:Field, 3:Mult Matrix, 4:Mult Function, 5:Frobenius Matrix, Basis is e0, e1, z1, z2"); 
end;

###Exceptional Split Spin Algebra

SpEx:=function(dt,FF)
local F, d, mas, prod, e0, e1, z1, A, frob;

if dt in FF then F:=FF;
else F:=FunctionField(FF,["dt"]);
fi;

## order of the basis: e0, e1, z1,

d:=dt^0;

mas :=d * [
[[0,0,-3],[0,0,-3*dt],[-1,0,0]],
[[0,0,-3*dt],[0,0,-3],[0,-1,0]], 
[[-1,0,0],[0,-1,0],[0,0,1]]];

prod := function(u,v)
local i,j,k, ans;
ans:=d * [0,0,0];
for i in [1..3] do
for j in [1..3] do
for k in [1..3] do
ans[k] := ans[k] + u[i]*v[j]*mas[i][j][k];
od;
od;
od;
return ans;
end;

e0:=d*[1,0,0];;
e1:=d*[0,1,0];;
z1:=d*[0,0,1];;

A:=VectorSpace(F,[e0,e1,z1]);;

frob:=[[3,3*dt,0],[3*dt,3,0],[0,0,1]];

return [A,F,mas,prod,frob];;
end;

SpExDetails:=function();
Print("1:Vector Space, 2:Field, 3:Mult Matrix, 4:Mult Function, 5:Frobenius Matrix, Basis is e0,e1,z1"); 
end;

###Exceptional Split Spin Algebra Cover

SpExC:=function(dt,FF)
local F, d, mas, prod, e0, e1, z1, n, A, frob;

if dt in FF then F:=FF;
else F:=FunctionField(FF,["dt"]);
fi;

## order of the basis: e0, e1, z1, n

d:=dt^0;

mas :=d * [
[[0,0,-3,2],[0,0,-3*dt,2*dt],[-1,0,0,0],[0,0,0,0]],
[[0,0,-3*dt,2*dt],[0,0,-3,2],[0,-1,0,0],[0,0,0,0]], 
[[-1,0,0,0],[0,-1,0,0], [0,0,1,0],[0,0,0,0]],
[[0,0,0,0],[0,0,0,0], [0,0,0,0],[0,0,0,0]]];

prod := function(u,v)
local i,j,k, ans;
ans:=d * [0,0,0,0];
for i in [1..4] do
for j in [1..4] do
for k in [1..4] do
ans[k] := ans[k] + u[i]*v[j]*mas[i][j][k];
od;
od;
od;
return ans;
end;

e0:=d*[1,0,0,0];;
e1:=d*[0,1,0,0];;
z1:=d*[0,0,1,0];;
n:=d*[0,0,0,1];;

A:=VectorSpace(F,[e0,e1,z1,n]);;

frob:=[[3,3*dt,0,0],[3*dt,3,0,0],[0,0,1,0],[0,0,0,0]];

return [A,F,mas,prod,frob];;
end;

SpExCDetails:=function();
Print("1:Vector Space, 2:Field, 3:Mult Matrix, 4:Mult Function, 5:Frobenius Matrix, Basis is e0,e1,z1,n"); 
end;

###Non-Symmetric Algebras

###Q2

Q2:=function(eta,FF)
local F, z, mas, prod, s1, s2, d1, d2, A, frob;

F:=FunctionField(FF,["eta"]);;

## order of the basis: s1, s2, d1, d2

z:=eta^0;

mas :=z * [
[[1,0,0,0],[0,0,0,0],[eta,0,eta/2,-eta/2],[eta,0,-eta/2,eta/2]], [[0,0,0,0],[0,1,0,0],[0,eta,eta/2,-eta/2],[0,eta,-eta/2,eta/2]], 
[[eta,0,eta/2,-eta/2],[0,eta,eta/2,-eta/2], [0,0,1,0],[-eta,-eta,eta,eta]],[[eta,0,-eta/2,eta/2],[0,eta,-eta/2,eta/2],[-eta,-eta,eta,eta],[0,0,0,1]]];

prod := function(u,v)
local i,j,k, ans;
ans:=z * [0,0,0,0];
for i in [1..4] do
for j in [1..4] do
for k in [1..4] do
ans[k] := ans[k] + u[i]*v[j]*mas[i][j][k];
od;
od;
od;
return ans;
end;

s1:=z*[1,0,0,0];;
s2:=z*[0,1,0,0];;
d1:=z*[0,0,1,0];;
d2:=z*[0,0,0,1];;

A:=VectorSpace(F,[s1,s2,d1,d2]);;

frob:=[[1,0,eta,eta],[0,1,eta,eta],[eta,eta,2,2*eta],[eta,eta,2*eta,2]];

return [A,F,mas,prod,frob];
end;
