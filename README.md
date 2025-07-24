# ThesisBinary
This is the code used in my thesis 'Binary Axial Algebras of Monster Type.' The main goal of this code is to check certain algebras are axial. It also includes stuff relating to the radical and ideals.


The algebras constructed are accompanied with lots of information. For an algebra, we output a vector space, the field, multiplication in the termos of an array of matricies, multiplication in terms of a binary operator, and a Frobenius form. We also make a Information function for each algebra to state the basis and any more information through text. 

We also provide functions to check certain elements have Jordan or Monster fusion laws. 

For example, we check that the algebra M21 (See Bin(M21) in the thesis) has a primitive M(2,1/2)-axis. 
gap> F:=Rationals;;
gap> A:=M21(F);;
gap> a0:=[1,0,0,0,0,0,0];;
gap> a1:=[0,0,1,0,0,0,0];;
gap> a2:=[0,0,0,0,1,0,0];;
gap> MonsterCheck(A,a0,2,1/2);
Primitive Monster axis

We also add some functions on calculating the orthogonal complement with respect to the Frobenius form. 
