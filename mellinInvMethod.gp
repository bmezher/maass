\\ compute GL(3) Whittaker function by taking the inverse
\\ Mellin transform of the product of gamma functions

\\ abscissas and weights for integration
\\ abscissas are k*hh, weights are weightsuplist[k+1] for k>=0
\\ weightsdownlist[-k] for k<0
weightsuplist;
weightsdownlist;

\\ store the values of gamma-factors
\\ one entry for each quadrant, 3rd,4th,2nd,1st
gammalist = [vector(600),vector(600),vector(600),vector(600)];


\\ see notes, case 4 of whittaker function
\\ ssigma[1]=\frac{\sigma_1}{2}
\\ ssigma[2]=\sigma_2 - \frac{\sigma_1}{2}
\\ might want to play with this later
ssigma = [6,2];
\\ssigma = [3,3];

\\ step size
hh = 1/2;

\\ parameters coming from sym square lifts
makerfromlift(R) = [R,R]*2/3;

\\sets up the parameters for integrals given a vector of eigenvalues
mellinSetup(r,prec = default(realprecision)) = {
   local(scale,maxcoeff,tempprec,gg,i,eps,small=vector(2),kk,jj,s1only,num,denom,numindex,denomindex,gammaV,ss=vector(2),up,right);
   
   \\ compute vector of gamma parameters
   gammaV = vecsort([-r[1]-2*r[2],-r[1]+r[2],r[2]+2*r[1]])*I/2;
   print("Setting up integration for ", gammaV,": ");
   
   \\ increase precision
   tempprec = default(realprecision);
   default(realprecision,20+prec);

   gammaV = precision(gammaV,default(realprecision));

   \\ integration step   
   \\ hh=6/prec;
   
   \\ what error can we tolerate
   eps = 10^(-prec-5);
   
   \\ scaling factor
   scale = (abs(gammaV[1])+abs(gammaV[3]))*Pi;

   \\ how many things do we think we need to precompute
   maxcoeff = ceil(200/hh);
   
   \\ stuff that depends on ss[1] only   
   print1("Computing gamma factors that depend on s1 only: ");
   s1only = vector(2);
   s1only[1] = vector(maxcoeff,i,ss[1]=ssigma[1]+(i-1)*hh*I;scale+lngamma(ss[1]+gammaV[1])+lngamma(ss[1]+gammaV[2])+lngamma(ss[1]+gammaV[3]));
   s1only[2] = vector(maxcoeff,i,ss[1]=ssigma[1]-i*hh*I;scale+lngamma(ss[1]+gammaV[1])+lngamma(ss[1]+gammaV[2])+lngamma(ss[1]+gammaV[3]));
   print("OK");

   \\ gamma factors in the numerator
   print1("Computing other gamma factors in the numerator: ");
   num = vector(maxcoeff*4,i,ss = (ssigma[1]+ssigma[2]+(i-2*maxcoeff)*hh*I)/2;lngamma(ss-gammaV[1])+lngamma(ss-gammaV[2])+lngamma(ss-gammaV[3]));
   print("OK");
   
   \\ gamma factor in the denominator
   print1("Computing gamma factors in the denominator: ");
   denom = vector(maxcoeff*8,i,ss = (3*ssigma[1]+ssigma[2]+(i-4*maxcoeff)*hh*I)/2;lngamma(ss));
   print("OK");
   
   print1("Precomputing all integration weights: ");
   for(up=0,1, 
      small[2] = 0;
      kk=0;

      numindex = 2*maxcoeff-up;\\      if(up==1,numindex--)
      denomindex = 4*maxcoeff-up;
         
      while(small[2]<3,
         kk++;
         
         if(up==1,numindex++;denomindex++,numindex--;denomindex--);
      
         for(right=0,1,
            \\ what quadrant are we in now?
            i = 1 + 2*up + right;
           
            small[1]=0;

            jj = 0;
            tempvector = vector(maxcoeff);
            while(small[1]<2,
               jj++;
 \\              print("numindex:", numindex+if(right==1,jj-1,-jj));
               gg = exp(s1only[2-right][jj]+num[numindex+if(right==1,jj-1,-jj)]-denom[denomindex+3*if(right==1,jj-1,-jj)]);
               
               tempvector[jj] = gg;
               if(abs(gg)<eps,small[1]++,small[1]=0);
            );
            gammalist[i][kk]=vector(jj,j,tempvector[j]);
         );
      
         if(length(gammalist[i][kk])<3,small[2]++,small[2]=0);

      );
      if(up==0,weightsdownlist = vector(kk,i,0),weightsuplist = vector(kk,i,0));
   );
   
   gammalist *= hh^2/(8*Pi^2); \\ extra factor of 4 in the denominator for consistency with other versions
   
   \\ reset precision
   default(realprecision,tempprec);

   print("OK");

   return;
}

computeWeights(d) = {
   local(prec,dmult,dd);
   prec = precision(gammalist);
   d = precision(d,prec);
   d *= Pi^3;
   
   for(i=1,length(weightsdownlist),
      weightsdownlist[i]=0;
      \\right
      dd = d^(-I*hh); dmult = 1; for(j=1,length(gammalist[2][i]),weightsdownlist[i]+= dmult*gammalist[2][i][j];dmult *= dd);
      \\left
      dd = 1/dd; dmult = dd; for(j=1,length(gammalist[1][i]),weightsdownlist[i]+= dmult*gammalist[1][i][j];dmult *= dd);
   );

   for(i=1,length(weightsuplist),
      weightsuplist[i]=0;
      \\right
      dd = d^(-I*hh); dmult = 1; for(j=1,length(gammalist[4][i]),weightsuplist[i]+= dmult*gammalist[4][i][j];dmult *= dd);
      \\left
      dd = 1/dd; dmult = dd; for(j=1,length(gammalist[3][i]),weightsuplist[i]+= dmult*gammalist[3][i][j];dmult *= dd);
   );


   dmult = d^(1/2-ssigma[1]);
   weightsdownlist *= dmult;
   weightsuplist *= dmult;

}

\\ evaluates mellin transform at \pi y
whittakerEvaluate(y) = {
   local(result,prec,ymult,yy);
   prec = precision(weightsuplist);
   y = precision(y,prec);
   y *= Pi;
   
   result = 0;
   yy = y^(I*hh);
   ymult = yy;

   for(i=1,length(weightsdownlist),
      result += ymult*weightsdownlist[i];
      ymult *= yy
   );
   yy = 1/yy;
   ymult = 1;
   for(i=1,length(weightsuplist),
      result += ymult*weightsuplist[i];
      ymult *= yy
   );

   result *= y^(1/2-ssigma[2]);

\\   print("y2: ",y/Pi," -> ", real(result));
   
   return (result);
}
