/*Shallue-van de Woestijne map for BW13-P310*/
z:=-2224;
p:=1749234309176102157657582860550885176950582224007184238236721873530271444092780\
387026731606667;
r:=2143085360734996117913472445644488914854141302995428209978212956147876055499508\
01;
t:=-72424742659885778123097924206425573989051009199;
Fp:=GF(p);
F13<u>:=ExtensionField<Fp,u|u^13+2>;
E:=EllipticCurve([F13|0,-17]);
b:=-17;
t:=Random(Fp);
t:=t*u;
sqrt3:=Sqrt(Fp!(-3));
j:=(-1+sqrt3)/2;
w:=sqrt3*(t/(1+b+t^2));
x1:=j-t*w;
x2:=-1-x1;
x3:=1+1/w^2;
a1:=-1;a2:=-1;
if IsSquare(x1^3+b) eq true  then 
a1:=1;
end if;
if IsSquare(x2^3+b) eq true  then 
a2:=1;
end if;
x:=(a1-1)*a2 mod 3 +1;
if x eq 1 then 
  Sqrt(x1^3+b);
elif  x eq 2 then
 Sqrt(x2^3+b);
elif  x eq 3 then
 Sqrt(x3^3+b);
end if;
x;


   
