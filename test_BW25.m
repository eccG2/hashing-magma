/*The code is used to implement hashing to G2 on BW19-P286
using Method II.
*/
z:=6995;
p:=(z^52-z^51+z^50-z^27-2*z^26-z^25+z^2 + 2*z + 1)div 3;
r:=z^40 + z^35 - z^25 - z^20 - z^15 + z^5 + 1;
t:=-z^26 + z + 1;
Fp:=GF(p);
F25:=GF(p^25);
/*w1,w2 are two primitive cube root over Fp*/
w1:=1735097437082519757632288525845952039444587736336558963475975342275380694134291\
2274367063789025545039268507465847101947098962640750317042393482806722969352351\
84585374973743655181055516217618729;
w2:=w1^2 mod p;

E:=EllipticCurve([F25|0,31]);

/* endo is the efficiently computable enomorphism \psi^i=(\pi*\phi)^2, w
where \phi:(x,y)->(w1*x,y)*/
endo:=function(Q,i);
      if i mod 3 eq 0 then 
         R:=E![Frobenius(Q[1],Fp,i), Frobenius(Q[2],Fp,i),1];
      elif i mod 3 eq 1 then 
         R:=E![Frobenius(Q[1],Fp,i)*w1, Frobenius(Q[2],Fp,i),1];
     else
         R:=E![Frobenius(Q[1],Fp,i)*w2, Frobenius(Q[2],Fp,i),1];
    end if;
    return R;
end function;

/*computing Tr(Q)-k*Q*/



Q:=Random(E);

/*hasing to G2 using Method I*/

    P:=E![Frobenius(Q[1],Fp,5),Frobenius(Q[2],Fp,5),1]-Q;
    P1:=z*P;P2:=P1-P;P3:=z*P1;
    R26:=-2*P;R24:=P2-P3;
    R23:=-z*R24; R22:=-z*R23;R21:=-z*R22;R20:=-z*R21;
    R19:=-z*R20;R18:=-z*R19;R17:=-z*R18;R16:=-z*R17;
    R15:=-z*R16;R14:=-z*R15;R13:=-z*R14;R12:=-z*R13;
    R11:=-z*R12;R10:=-z*R11;R9:=-z*R10;R8:=-z*R9;
    R7:=-z*R8;R6:=-z*R7;R5:=-z*R6;R4:=-z*R5;
    R3:=-z*R4;R2:=-z*R3;R1:=-z*R2-P;
    U0:=(R1+R4+R7+R10+R13+R16+R19+R22+P1);
    U1:=(R3+R6+R9+R12+R15+R18+R21+R24);
    U2:=(R2+R5+R8+R11+R14+R17+R20+R23);
    R0:=U0-U1;R25:=U0-U2+P2;
    R:=R0+endo(R1,1)+endo(R2,2)+endo(R3,3)+endo(R4,4)+endo(R5,5)+endo(R6,6)+endo(R7,7)+endo(R8,8)+endo(R9,9)
    +endo(R10,10)+endo(R11,11)+endo(R12,12)+endo(R13,13)+endo(R14,14)+endo(R15,15)+endo(R16,16) +endo(R17,17)
    +endo(R18,18)+endo(R19,19)+endo(R20,20)+endo(R21,21)+endo(R22,22)+endo(R23,23)+endo(R24,24)+endo(R25,25)+endo(R26,26);
T0:=R;T1:=r*R;T2:=E![Frobenius(R[1],Fp,1),Frobenius(R[2],Fp,1),1]-p*R;
if T0[3] ne 0 and T1[3] eq 0 and T2[3] eq 0 then
    printf"CORRECT-- Method II!\n";
elif  T0[3] eq 0 then 
    printf"Select random point again!\n" ;
else
    printf"ERROR --Method II! ";
end if;

