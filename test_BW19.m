/*The code is used to implement hashing to G2 on BW19-P286
using Method I and II.
*/
z:=-145;
p:=(z+1)^2*(z^38-z^19+1)div 3-z^39;
t:=-z^20+z+1;
r:=Evaluate(CyclotomicPolynomial(114),z);
Fp:=GF(p);
F19<u>:=ExtensionField<Fp,u|u^19+u+31>;
E:=EllipticCurve([F19|0,31]);
/*w1,w2 are two primitive cube root over Fp*/
w1:=1364614900770326371786571581106544735789185132426333114702258894978661275889775\
5244;
w2:=w1^2 mod p;

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

/*hasing to G2 using Method I*/
Q:=Random(E);
P:=E![Frobenius(Q[1],Fp,1),Frobenius(Q[2],Fp,1),1]-Q;
   S1:=z*P;S2:=z*S1; S3:=z*S2; S4:=z*S3;S5:=z*S4;S6:=z*S5;S7:=z*S6;
    S8:=z*S7;S9:=z*S8;S10:=z*S9;S11:=z*S10;S12:=z*S11;S13:=z*S12;S14:=z*S13;
S15:=z*S14;S16:=z*S15;S17:=z*S16;S18:=z*S17;S19:=z*S18;S20:=z*S19;
 M0:=(2*S2+S1);M1:=(2*S5+S4);M2:=(2*S8+S7);M3:=(2*S11+S10);M4:=(2*S14+S13);M5:=(2*S17+S16);
 R0:=M0-P;H0:=2*S3-M0;R1:=M1-S3;H1:=2*S6-M1;
 R2:=M2-S6;H2:=2*S9-M2;R3:=M3-S9;H3:=2*S12-M3;
 R4:=M4-S12;H4:=2*S15-M4;R5:=M5-S15;H5:=2*S18-M5;

  L0:=(S20-S19+S18);L1:=z*L0;L2:=z*L1;L3:=z*L2;L4:=z*L3;
    L5:=z*L4;L6:=z*L5; L7:=z*L6;L8:=z*L7;L9:=z*L8;
    L10:=z*L9;L11:=z*L10;L12:=z*L11;L13:=z*L12;L14:=z*L13;
L15:=z*L14;L16:=z*L15;L17:=z*L16;L18:=z*L17;
K0:=L14+H4;K3:=-(L11+H3);K6:=(L8+H2);K9:=-(L5+H1);K12:=(L2+H0);
K1:=R4;K4:=-R3;K7:=R2;K10:=-R1;K13:=R0;
K2:=-(L12-3*S12); K5:=(L9-3*S9);K8:=-(L6-3*S6);K11:=(L3-3*S3);K14:=-(L0-3*P);
K15:=-L18+3*(S18-P);
K16:=-L17-H5;
K17:=-R5;
K18:=L15-3*S15;

  /*R in the target point*/
R:=K0+E![Frobenius(K1[1],Fp,1),Frobenius(K1[2],Fp,1)]+E![Frobenius(K2[1],Fp,2),Frobenius(K2[2],Fp,2)]
+E![Frobenius(K3[1],Fp,3),Frobenius(K3[2],Fp,3)]+E![Frobenius(K4[1],Fp,4),Frobenius(K4[2],Fp,4)]
+E![Frobenius(K5[1],Fp,5),Frobenius(K5[2],Fp,5)]+E![Frobenius(K6[1],Fp,6),Frobenius(K6[2],Fp,6)]
+E![Frobenius(K7[1],Fp,7),Frobenius(K7[2],Fp,7)]+E![Frobenius(K8[1],Fp,8),Frobenius(K8[2],Fp,8)]
+E![Frobenius(K9[1],Fp,9),Frobenius(K9[2],Fp,9)]+E![Frobenius(K10[1],Fp,10),Frobenius(K10[2],Fp,10)]
+E![Frobenius(K11[1],Fp,11),Frobenius(K11[2],Fp,11)]+E![Frobenius(K12[1],Fp,12),Frobenius(K12[2],Fp,12)]
+E![Frobenius(K13[1],Fp,13),Frobenius(K13[2],Fp,13)]+E![Frobenius(K14[1],Fp,14),Frobenius(K14[2],Fp,14)]
+E![Frobenius(K15[1],Fp,15),Frobenius(K15[2],Fp,15)]+E![Frobenius(K16[1],Fp,16),Frobenius(K16[2],Fp,16)]
+E![Frobenius(K17[1],Fp,17),Frobenius(K17[2],Fp,17)]+E![Frobenius(K18[1],Fp,18),Frobenius(K18[2],Fp,18)];
/*Check method II is correct*/
T0:=R;T1:=r*R;T2:=E![Frobenius(R[1],Fp,1),Frobenius(R[2],Fp,1),1]-p*R;
if T0[3] ne 0 and T1[3] eq 0 and T2[3] eq 0 then
    printf"CORRECT--Method I!\n";
elif  T0[3] eq 0 then 
    printf"Reselect random point!" ;
else
    printf"ERROR -- Method I! ";
end if;


/*hasing to G2 using Method II*/

    P:=E![Frobenius(Q[1],Fp,1),Frobenius(Q[2],Fp,1),1]-Q;
    P1:=z*P;P2:=P1-P;P3:=z*P1;
    R20:=2*P; R18:=P3-P2;
 R17:=-z*R18;R16:=-z*R17;R15:=-z*R16;R14:=-z*R15;
 R13:=-z*R14;R12:=-z*R13;
    R11:=-z*R12;R10:=-z*R11;R9:=-z*R10;R8:=-z*R9;
    R7:=-z*R8;R6:=-z*R7;R5:=-z*R6;R4:=-z*R5;
    R3:=-z*R4;R2:=-z*R3;R1:=-z*R2+P;
    U0:=(R1+R4+R7+R10+R13+R16-P1);
    U1:=(R3+R6+R9+R12+R15+R18);
    U2:=(R2+R5+R8+R11+R14+R17);
    R0:=U0-U1;R19:=U0-U2-P2;
   /*R in the target point*/
    R:=R0+endo(R1,1)+endo(R2,2)+endo(R3,3)+endo(R4,4)+endo(R5,5)+endo(R6,6)+endo(R7,7)+endo(R8,8)
    +endo(R9,9)+endo(R10,10)+endo(R11,11)+endo(R12,12)+endo(R13,13)+endo(R14,14)+endo(R15,15)
    +endo(R16,16)+endo(R17,17)+endo(R18,18)+endo(R19,19)+endo(R20,20);
 
/*Check method II is correct*/
T0:=R;T1:=r*R;T2:=E![Frobenius(R[1],Fp,1),Frobenius(R[2],Fp,1),1]-p*R;
if T0[3] ne 0 and T1[3] eq 0 and T2[3] eq 0 then
    printf"CORRECT--Method II!\n";
elif  T0[3] eq 0 then 
    printf"Reselect random point!" ;
else
    printf"ERROR -- Method II! ";
end if;
