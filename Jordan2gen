###Here we provide the 2-generated axial algebras of Jordan type. 
###These were classified in 'rimitive axial algebras of Jordan type' by Hall, Rehren and Shpectorov.

###2B
2B:=function(F)
local mas, prod, a, b, A, frob;
## order of the basis: a, b
mas :=[[[1,0],[0,0]],[[0,0],[0,1]]];

prod := function(u,v)
local i,j,k, ans;
ans:=[0,0];
for i in [1..2] do
for j in [1..2] do
for k in [1..2] do
ans[k] := ans[k] + u[i]*v[j]*mas[i][j][k];
od;
od;
od;
return ans;
end;

a:=[1,0];;
b:=[0,1];;

A:=VectorSpace(F,[a,b]);;

frob:=[[1,0],[0,1]];

return [A,F,mas,prod,frob];
end;

2BDetails:=function();
Print("1:Vector Space, 2:Field, 3:Mult Matrix, 4:Mult Function, 5:Frobenius Matrix, Basis is a,b"); 
end;

##3C(-1)^x

3CX:=function(F)
local mas, prod, x, y, A, frob;

## order of the basis: x, y

mas :=[[[1,0],[-1,-1]],[[-1,-1],[0,1]]];

prod := function(u,v)
local i,j,k, ans;
ans:=[0,0];
for i in [1..2] do
for j in [1..2] do
for k in [1..2] do
ans[k] := ans[k] + u[i]*v[j]*mas[i][j][k];
od;
od;
od;
return ans;
end;

x:=[1,0];;
y:=[0,1];;

A:=VectorSpace(F,[x,y]);;

frob:=[[1,-1/2],[-1/2,1]];;

return [A,F,mas,prod, frob];
end;

3CXDetails:=function();
Print("1:Vector Space, 2:Field, 3:Mult Matrix, 4:Mult Function, 5:Frobenius Matrix, Basis is x,y"); 
end;

###3C(eta)

3C:=function(eta,FF)
local d, mas, prod, x, y, z, A, F, frob;

if eta in FF then F:=FF;
else F:=FunctionField(FF,["eta"]);
fi;

## order of the basis: x, y, z

d:=eta^0;

mas :=d*[[[1,0,0],[eta/2,eta/2,-eta/2],[eta/2,-eta/2,eta/2]],[[eta/2,eta/2,-eta/2],[0,1,0],[-eta/2,eta/2,eta/2]],[[eta/2,-eta/2,eta/2],[-eta/2,eta/2,eta/2],[0,0,1]]];

prod := function(u,v)
local i,j,k, ans;
ans:=d*[0,0,0];
for i in [1..3] do
for j in [1..3] do
for k in [1..3] do
ans[k] := ans[k] + u[i]*v[j]*mas[i][j][k];
od;
od;
od;
return ans;
end;

x:=d*[1,0,0];;
y:=d*[0,1,0];;
z:=d*[0,0,1];;

A:=VectorSpace(F,[x,y,z]);;

frob:=[[1,eta/2,eta/2],[eta/2,1,eta/2],[eta/2,eta/2,1]];

return [A,F,mas,prod,frob];;
end;

3CDetails:=function();
Print("1:Vector Space, 2:Field, 3:Mult Matrix, 4:Mult Function, 5:Frobenius Matrix, Basis is x,y,z"); 
end;

###SpinFactor

Spin:=function(dt,FF)
local z, mas, prod, a, b, id, A, F, frob;

if dt in FF then F:=FF;
else F:=FunctionField(FF,["dt"]);
fi;

## order of the basis: a, b, id

z:=dt^0;

mas :=z*[[[1,0,0],[1/2,1/2,(dt-2)/8],[1,0,0]],[[1/2,1/2,(dt-2)/8],[0,1,0],[0,1,0]],[[1,0,0],[0,1,0],[0,0,1]]];

prod := function(u,v)
local i,j,k, ans;
ans:=z*[0,0,0];
for i in [1..3] do
for j in [1..3] do
for k in [1..3] do
ans[k] := ans[k] + u[i]*v[j]*mas[i][j][k];
od;
od;
od;
return ans;
end;

a:=z*[1,0,0];;
b:=z*[0,1,0];;
id:=z*[0,0,1];;

A:=VectorSpace(F,[a,b,id]);;

frob:=[[1,(dt+2)/4,1],[(dt+2)/4,1,1],[1,1,2]];

return [A,F,mas,prod,frob];;
end;

SpinDetails:=function();
Print("1:Vector Space, 2:Field, 3:Mult Matrix, 4:Mult Function, 5:Frobenius Matrix, Basis is a,b,id"); 
end;

###S(2)^circ
S2c:=function(F)
local mas, prod, a, b, A, frob;

## order of the basis: a, b

mas :=[[[1,0],[1/2,1/2]],[[1/2,1/2],[0,1]]];

prod := function(u,v)

local i,j,k, ans;

ans:=[0,0];

for i in [1..2] do
for j in [1..2] do
for k in [1..2] do
ans[k] := ans[k] + u[i]*v[j]*mas[i][j][k];
od;
od;
od;
return ans;
end;

a:=[1,0];;
b:=[0,1];;

A:=VectorSpace(F,[a,b]);;

frob:=[[1,1],[1,1]];

return [A,F,mas,prod,frob];;
end;

###S(hat)(2)^circ

S2cDetails:=function();
Print("1:Vector Space, 2:Field, 3:Mult Matrix, 4:Mult Function, 5:Frobenius Matrix, Basis is a,b"); 
end;

S2:=function(F)
local mas, prod, a, b, n, A, frob;

## order of the basis: a, b, n

mas :=[[[1,0,0],[1/2,1/2,1],[0,0,0]],[[1/2,1/2,1],[0,1,0],[0,0,0]],[[0,0,0],[0,0,0],[0,0,0]]];

prod := function(u,v)
local i,j,k, ans;
ans:=[0,0,0];
for i in [1..3] do
for j in [1..3] do
for k in [1..3] do
ans[k] := ans[k] + u[i]*v[j]*mas[i][j][k];
od;
od;
od;
return ans;
end;

a:=[1,0,0];;
b:=[0,1,0];;
n:=[0,0,1];;

A:=VectorSpace(F,[a,b,n]);;

frob:=[[1,1,0],[1,1,0],[0,0,0]];

return [A,F,mas,prod,frob];;
end;

S2Details:=function();
Print("1:Vector Space, 2:Field, 3:Mult Matrix, 4:Mult Function, 5:Frobenius Matrix, Basis is a,b,n"); 
end;












