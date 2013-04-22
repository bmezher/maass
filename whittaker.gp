\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ define some useful constants                                       \\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
r=vector(2);		\\ eigenvalues

\\ vector of gamma factors
\\ the gamma factors are \Gamma(\frac{s+gammaV[i]}{2})
gammaV=vector(3);   

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\  Series at 0                                                       \\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ here we store factors for the recursions of the coefficients
\\ coefficients for each of the 6 series are of the form
\\                                   (1-mnfactor)_(m+n)
\\ c_{m,n} = -----------------------------------------------------------------
\\           (1-mnfactor)_(m) (1-mfactor)_m (1-mnfactor)_n (1-nfactor)_n n! m!
\\
global(mnfactor);
global(mfactor);
global(nfactor);

cMN = matrix(600,600); \\ store coefficients 
global(cMNn); \\ what's the largest n
global(cMNm); \\ what's the largest m

\\ powers of y1 and y2 that premultiply each of 6 series
global(y1power);
global(y2power);
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ set precision defaults                                             \\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
answerprecision = 25;
default(realprecision,answerprecision);
termprecision = 35;
taylorprecision = 50;
integralprecision = 30;
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Polynomial recursion with bessel function                          \\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ won't need more than 100 coefficients
polyP = vector(300);
polyQ = vector(300);
polyC = vector(300);

\\what's the largest n?
global(maxN);

\\factors a_n and b_n used in recursion
global(recA,recB);

\\ vector of mu^2's used in the recursion
global(bigMu);

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Stade's integral and bessel functions                              \\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
global(mmu); global(nu);

h=0.05; \\ integration step
global(exponentVector); \\ pre-compute exponents
global(sqrtPlusVector); \\ pre-compute square roots inside bessel
global(sqrtMinusVector);

\\ besselC; \\coeffs

debug=0;

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\whittaker function is virtually zero after this 
whCutOff = 10^10;

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ some helpful functions                                             \\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ given a vector of eigenvalues, initialize the appropriate constants
whittakerInit(v)=
{
   local(gArgs,convFactor,cutoff);

   default(realprecision,taylorprecision);

   \\ v is the vector (r1,r2)
   \\ we have (\nu_1,\nu_2) = (1/3+3/2*i*r_1,1/3+3/2*i*r_2)
   r[1] = 3*v[1]*I/2;
   r[2] = 3*v[2]*I/2;

   \\ initialize coefficient factors
   mnfactor = [r[1]+r[2],r[2],r[1],-r[1]-r[2],-r[2],-r[1]];
   mfactor = [r[1],-r[1],r[1]+r[2],-r[2],-r[1]-r[2],r[2]];
   nfactor = [r[2],r[1]+r[2],-r[2],-r[1],r[1],-r[1]-r[2]];
   
   \\ initialize powers
   y1power = (2/3)*[-r[1]-2*r[2], -r[1]-2*r[2], -r[1]+r[2], 2*r[1]+r[2], -r[1]+r[2], 2*r[1]+r[2]];
   y2power = (2/3)*[-r[2]-2*r[1], r[1]-r[2], -2*r[1]-r[2], r[1]+2*r[2], r[1]+2*r[2], r[1]-r[2]];

   \\initialize coefficients
   cMNn = 0;
   cMNm = 0;

   \\ factor out gamma(r[1])*gamma(r[2])*gamma(r[1]+r[2])
   \\ this is what remains of each gamma factor
   gArgs = -2*I*[arg(gamma(r[1])), arg(gamma(r[2])), arg(gamma(r[1]+r[2]))];
   
   cMN[1,1] = [1, exp(gArgs[1]), exp(gArgs[2]), exp(gArgs*[1,1,1]~),
                exp(gArgs*[0,1,1]~), exp(gArgs*[1,0,1]~)];

   \\ we want this to coincide with the other method
   \\ so we multiply it by gammas and e
   convFactor = lngamma(r[1])+lngamma(r[2])+lngamma(r[1]+r[2])-(r[1]+r[2])*Pi*I;
   cMN *= exp(convFactor);


   \\let's get this party started
   recurseN();  recurseN(); recurseN(); recurseN(); recurseN();
   recurseM();  recurseM(); recurseM(); recurseM(); recurseM();

   \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
   \\ initialize stuff for polynomial recursion

   \\ this is the (\alpha,\beta,\gamma) from the notes
   gammaV = 2*[-(r[1]+2*r[2]),r[2]-r[1],2*r[1]+r[2]]/3;
   
   \\ factors a_n and b_n used in recursion
   recA = 1;
   recB = 3*gammaV/2;

   \\ P_0=4, Q_0 = 0
   polyP[1] = [2,2,2]; \\2 is the right normalization for some reason
   polyQ[1] = [0,0,0];

   \\ initialize mu^2 vector
   \\ bigMu = [(gammaV[2]-gammaV[3])^2,(gammaV[1]-gammaV[3])^2,(gammaV[1]-gammaV[2])^2]/4;
   bigMu = [r[1]^2,(r[1]+r[2])^2,r[2]^2];

   \\ what's the largest n?
   maxN = 1;

   \\ the constant term is gamma(gammav[i]-gammav[j]) for all possible i!=j
   \\ then we scale everything by exp(Pi*(r[1]+r[2]))
   \\polyC[1]=exp([lngamma(gammaV[2]-gammaV[1])+lngamma(gammaV[3]-gammaV[1]),\\
   \\              lngamma(gammaV[1]-gammaV[2])+lngamma(gammaV[3]-gammaV[2]),\\
   \\              lngamma(gammaV[1]-gammaV[3])+lngamma(gammaV[2]-gammaV[3])]\\
   \\              -[1,1,1]*((r[1]+r[2])*Pi*I))
   polyC[1]=exp([lngamma(r[2])+lngamma(r[1]+r[2]),lngamma(-r[2])+lngamma(r[1]),lngamma(-(r[1]+r[2]))+lngamma(-r[1])]-[1,1,1]*((r[1]+r[2])*Pi*I));

   for(i=1,99,recursePoly());


   \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
   \\ initialize stuff we'll need for the integral
   mmu = r[1]-r[2];
   nu = r[1]+r[2];

   \\pre-compute exponents we'll need
   exponentVector = vector(ceil(15/h),n,exp((n-1)*h*mmu));
   sqrtPlusVector = vector(ceil(15/h),n,sqrt(1+exp(2*h*(n-1))));
   sqrtMinusVector = vector(ceil(15/h),n,sqrt(1+exp(-2*h*(n-1))));
   
   \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
   \\ determine where to cut whittaker function off
   debug=1;
   cutoff = abs(nu/2);
   while(abs(whittaker([cutoff, 0.5]))>10^(-answerprecision),cutoff++);
   whCutOff = cutoff;

   default(realprecision,answerprecision);

}


\\ compute one more coefficient for N
recurseN()=
{     
   local(tempcN);

   \\ we are inserting a new column
   \\ take the last column and use recursion
   \\                         (m+n+1 - mnfactor)
   \\   c_{m,n+1} = -------------------------------------- * c_{m,n}
   \\               (n+1) (n+1 - mnfactor) (n+1 - nfactor)
   \\
   cMNn++;

   for(m=1,cMNm+1,cMN[m,cMNn+1]= vector(6,i,
      cMN[m,cMNn][i]*((m-1)+cMNn - mnfactor[i])/
      ((cMNn-mnfactor[i])*(cMNn-nfactor[i])) )/cMNn);

}

\\ compute one more coefficient for M
recurseM()= {     
   local(tempcM);

   \\ we are inserting a new row
   \\ take the last row and use recursions
   \\                         (m+n+1 - mnfactor)
   \\   c_{m+1,n} = -------------------------------------- * c_{m,n}
   \\               (m+1) (m+1 - mnfactor) (m+1 - mfactor)
   \\
   cMNm++;
   
   for(n=1,cMNn+1,cMN[cMNm+1,n]= vector(6,i,
      cMN[cMNm,n][i]*(cMNm+(n-1) - mnfactor[i])/
      ((cMNm-mnfactor[i])*(cMNm-mfactor[i])) )/cMNm);
}

\\ compute the next polynomial
recursePoly() = {
   \\increment recB
   recB+= [2,2,2];

   \\compute the next polynomials
   polyP[maxN+1] = y*deriv(polyP[maxN],y)+4*(Pi^2)*y^2*polyQ[maxN]+ptMult(bigMu,polyQ[maxN])+ptMult(recB,polyP[maxN]);
   polyQ[maxN+1] = polyP[maxN]+y*deriv(polyQ[maxN],y)+ptMult(recB,polyQ[maxN]);

   \\compute the next coeff   
   polyC[maxN+1] = ptMult(polyC[maxN],[1/((maxN-r[2])*(maxN-r[1]-r[2])),1/((maxN+r[2])*(maxN-r[1])),1/((maxN+r[1]+r[2])*(maxN+r[1]))]);   
 
   \\increment the count
   maxN++;
}

\\point-wise multiplication of two vectors
ptMult(vec1,vec2) = vector(length(vec1),i,vec1[i]*vec2[i]);

\\derivative of the K-Bessel function
besselkPrime(myNu,x) = -(besselk(myNu-1,x)+besselk(myNu+1,x))/2;

testCoeff(m,n) = {
   local(result);
   result = vector(6,i,gamma(nfactor[i]-n)*gamma(mnfactor[i]-n)*
                gamma(mfactor[i]-m)*gamma(mnfactor[i]-m)/
           (gamma(mnfactor[i]-m-n) * factorial(m) *factorial(n))  );
           
   if((m+n) %2,result*=-1);
   result;
}

\\ d = y1^2 y2, y2
myW(d,y) = whittaker([sqrt(d/y),y]);

\\ whittaker(y1,y2) using series expansion at y=0
whittakerSeries(y) =
{
   local(y1,y2,YPower,YPowerOld,res,nn,mm,change,flag);

   if(debug,print("Using Series"));

   \\ upgrade precision
   default(realprecision,taylorprecision);
   y = precision(y,taylorprecision);


   \\ y always comes with Pi
   y *= Pi;

   res=0; nn=0; change=1E100;

   \\ we'll need squares of y
   y2 = y[2]^2;
   y1 = y[1]^2;
      
   \\ initial values 
   YPower = vector(6,j,y[1]^(1+y1power[j])*y[2]^(1+y2power[j]));
   YPowerOld = YPower;


   flag = 0;

   until (flag,
     nn++;
     if(nn>cMNn,recurseN());

     \\ back up powers of y1 and y2; will need to kick powers of y2 back to 1
        
     mm = 0;
     until (abs(change)<10^-(taylorprecision),
        mm++;
        if(mm>cMNm,recurseM());

        change = YPower*cMN[mm,nn]~;

        res += change;
        YPower *= y2;
     );
     
     \\ we reached our precision goal
     if(mm<2,flag=1);

     \\ revert ypower, and increment power of y1
     YPowerOld *= y1;
     YPower=YPowerOld;
   );
   default(realprecision,answerprecision);
   res
}

\\ whittaker(y1,y2) using series expansion for y_1 (assumed small)
\\ and polynomial recursion in y_2 (assumed large)
whittakerPoly(yvec) =
{
   local(y1,y,YPower,Yterm,res,nn,change);

   if(debug,print("Using Poly"));

   \\ upgrade precision
   default(realprecision,taylorprecision);
   yvec = precision(yvec,taylorprecision);

   \\y1 always comes with Pi
   y1=Pi*yvec[1];
   y=yvec[2];


   res=0; nn=0; change=1E100;
   
   \\ YPower is going to contain the powers of Y times the K-bessel function
   \\ that we multiply by the polynomials P and Q
   YPower = vector(2);

   \\ y-factors in front of the sum 
   YPower[1] = vector(3,i,y1^(1+gammaV[i])*(Pi*y)^(1+gammaV[i]/2));

   \\multiply by the correct version of K-Bessel
   YPower[2] = 2*Pi*y*ptMult(YPower[1],[besselkPrime(r[1],2*Pi*y),besselkPrime(r[1]+r[2],2*Pi*y),besselkPrime(r[2],2*Pi*y)]);
   YPower[1] = ptMult(YPower[1],[besselk(r[1],2*Pi*y),besselk(r[1]+r[2],2*Pi*y),besselk(r[2],2*Pi*y)]);

   \\ we'll need half-squares of y1
   y1 = y1^2/2;
   

   until (abs(change)<10^-(taylorprecision),
      nn++;

      Yterm = ptMult(eval(polyP[nn]),YPower[1])+ ptMult(eval(polyQ[nn]),YPower[2]);

      change = Yterm*polyC[nn]~;

      res += change;
      YPower *= y1;
      YPower /= nn;
   );
     
   default(realprecision,answerprecision);
   res
}


\\ will return exp(-i*pi*r/2)*besselk(r,x)/sqrt(pi)  
\\ that is, exp(-i*pi*r/2-x)*(2x)^(r)*U(r+1/2,2r+1,2x)
myBesselK(r,x) = {
   local(result);
   \\ might want to take care of precision issues here...
   
\\   print("besslk(",r,",",x,")");
   
   x *= 2;


   result = hyperu(r+1/2,2*r+1,x);
   result *= x^r;
   \\ this thing has to be real
   result = real(result);

   result *= exp((-I*Pi*r-x)/2);
   result
}

\\whittaker(y1,y2) using stade's integral formula
whittakerStade(y) = {
   local(yfactor,result,leftBound,rightBound,term,base);

   if(debug,print("Using Stade"));

   \\ let's precompute some stuff
   y *= Pi;
   y *= 2;
\\   yfactor = 2*Pi*(y[1]^(1+mmu/3))*(y[2]^(1+2*mmu/3));
   yfactor = 2*Pi*(y[1]^(1+mmu/3))*(y[2]^(1-mmu/3));



   leftBound=1; rightBound=1;

   result = 0;

   \\ use trapezoid rule to integrate the left part
   term = 0.5*myBesselK(nu,y[2]*sqrt(2))*myBesselK(nu,y[1]*sqrt(2));
   base = term;
   while(abs(term/base)>10^(-termprecision), 
      \\ keep going further left if the integrand is large
      result += h*term/exponentVector[leftBound];
      leftBound++;
      term = myBesselK(nu,y[2]*sqrtPlusVector[leftBound])*myBesselK(nu,y[1]*sqrtMinusVector[leftBound]);
   );
   result += 0.5*h*term/exponentVector[leftBound]; \\ the endpoint weighs only half as much

   \\ use trapezoid rule to integrate the right part
   term = base;
   while(abs(term/base)>10^(-termprecision), 
        \\ keep going further right if the integrand is large
      result += h*term*exponentVector[rightBound];
      rightBound++;
      term = myBesselK(nu,y[2]*sqrtMinusVector[rightBound])*myBesselK(nu,y[1]*sqrtPlusVector[rightBound]);

   );
   result += 0.5*h*term*exponentVector[rightBound]; \\ the endpoint weighs only half as much


   if(debug,print("on [-",(leftBound-1)*h,",",(rightBound-1)*h,"]")); 

   result *= yfactor;
}

\\whittaker
whittaker(y) = {
   if(debug,print("y1=",y[1]," y2=",y[2]));

   \\ don't bother computing if it's gonna be too small anyway
   if(max(y[1],y[2])>whCutOff,return(0));

   if(y[1]*y[2]>40,
      whittakerStade(y), \\ both are large
      if(y[1]<y[2],         \\ if y_1 is smaller...
          whittakerPoly(y), \\ ...use whittakerPoly..
          conj(whittakerPoly([y[2],y[1]])) \\..else flip the y's and conjugate
      );
   );
}

\\dual function
whittakerDual(y) = whittaker([y[2],y[1]]);

makerfromlift(R) = [R,R]*2/3;
