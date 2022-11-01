/*The code is used to find the target vector h for hashing to G2 on BW13-P310
with Method I and II.
*/
/*Curve paramters*/
P<u>:=PolynomialRing(RationalField());
r:=P!CyclotomicPolynomial(78);
p:=1/3*(u+1)^2*(u^26-u^13+1)-u^27;
z:=-2224;
p:=1749234309176102157657582860550885176950582224007184238236721873530271444092780\
387026731606667;
r:=2143085360734996117913472445644488914854141302995428209978212956147876055499508\
01;
t:=-72424742659885778123097924206425573989051009199;
Fp:=GF(p); 
/*w is a cube root in F^*_p*/
w:=1060961839710140962051029043989073482206300122999381892347557062395216868763690\
225024943;

/*f is the conductor, i.e., t^2-4*p=-3*f^2*/
f:=43*1069*1543*340677927255779144657220930205829502307;
F13:=ext<Fp|u^13+2>;
E:=EllipticCurve([F13|0,-17]);
/*n is the large cofacor of the cyclotomic subgroup G0. In other words, #G_0=r*n.
Given a random point of Q in G_0, then $n*Q\in G_2$*/

n:=3829449736163031399356110014646477253722949499224952771796808311526882284721622\
6525362977813692040700124543771134399033985394310987087382177047993102558030633\
4118770998203112469254700908857130114479496362519790167986675682039976464992393\
7192261075795207044450734446337600975254534173311469296388258219480949553532981\
3584712597113506252391725991204165885267872301249352188839268343758113438924152\
7012172338101075047568719807680411669875073643468745720987082622769647119644685\
3518974478239955997008167587101106670263648202647771435149397680378857936916475\
6444764371373435794214471623071568027118476511569297479889166678583016683344933\
8508974458721686999528912068342165095442833670236992302066157114412659248677177\
7509353628154364568961757033235826243320590035374411793068434381606088516745870\
3771547355668175073589144273557533492991476545917600348712011388282544095223066\
7335388118061271134463808566549744570756793861197125672593923227809450685046369\
4968224256300346209546584470105129662738139446243928338507535589194507819849529\
767732848457;

/*f2_inv=(1/(2f) mod nr)i. The value can be obtained using the command XGCD(2*f, r*cof);*/

f2_inv:=-2548094967119561086559523596620129350097834973229266954410007570927106340899\
4320438535070878868471204378167403917493251442254820882119792194980813455457940\
3152148218039490001937523674937724894139028900179684368693649700475573204631789\
3208842644215505934451306537817708313370139476656283466760422659704201847424170\
5552852995954934642542026862427439734321664852234008030965846025730121800018461\
2798224456223381119551817754582539594278930056020665951351108763864896071294097\
2201595624983519854856502979006326097260423751405952796442761170466585653964161\
0962509533794195558115407508113469792231691898603166504769502643205160782204866\
5392723810495490320763560273181319497136793312406035244842539260549258484989427\
8169573776322178476119778961846540401639283210113052377074960614385530677695997\
0258365212230232484024125363755527479908083838241567909897349051025801503814887\
6773606662910257077847493278812544224281728226607585850134210233430013978984451\
0867685711051503220766489548144014427858765424598443620671432021723009856468531\
6690705281616062491964051866316000132873806956686074733523822698495627453884809\
5784247484489013;

/*Computing the parameter a  such that  \pi(P)=[a]P for all P\in G*/
P<x>:=PolynomialRing(Integers());
g:=x^2-t*x+p; 
q:=P!CyclotomicPolynomial(13);
R:=quo<P|g>;
h:=R!q;
H:=Coefficients(h);
a:=Solution(H[2],-H[1],r*n);

/*computing b such that (w*x,y)=b(x,y) for (x,y)\in G0*/
b:=((-f+2*a-t)*f2_inv) mod (n*r);

/*computing  s such that (w*x^p,y^p)=s(x,y) for (x,y)\in G0*/
s:=a*b mod (n*r);


/*computing the vector h of Method I */
B:=RMatrixSpace(Integers(),12,12)!0;
B[1][1]:=n;
for i:=2 to 12 do
 B[i][1]:=-a^(i-1) mod (r*n);B[i][i]:=1;
end for;
L := LatticeWithBasis(B);
h1:=ShortestVector(L);

/*check that (h1[1]+h1[2]*a+...h1[12]*a^11) mod r ne 0*/
sum:=0;
for i:=1 to 12 do
    sum:=sum+ h1[i]*a^(i-1);
end for;
if sum mod r eq 0 then 
    printf"Method I:The vector is trival, please seclect a new one \n";
else
    printf"The vector h of Method I:\n";
    h1;
end if;


/*computing the vector h of Method II */
B:=RMatrixSpace(Integers(),23,23)!0;
B[1][1]:=n;
for i:=2 to 23 do
 B[i][1]:=-s^(i-1) mod (r*n);B[i][i]:=1;
end for;
L := LatticeWithBasis(B);
h2:=ShortestVector(L);
/*check that (h2[1]+h2[2]*a+...h[23]*a^22) mod r ne 0*/
sum:=0;
for i:=1 to 23 do
    sum:=sum+ h2[i]*s^(i-1);
end for;
if sum mod r eq 0 then 
    printf"Method II:The vector is trival, please seclect a new one \n";
else
    printf"The vector h of Method II:\n";
    h2;
end if;
