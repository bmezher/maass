\\ given coefficients and eigenvalue, evaluates the Maass form F
\\ at any z = [x1, x2, x3, y1, y2]

\\ modify the location stuff below as needed
\\ set up globals as in this header, then use maassF to evaluate a Maass 
\\ form

\\ what do we consider 0?
eps = 10^(-8);

\\ need this in order to evaluate things
read("hecke.gp");
read("coeff9x9.txt");

\\ don't bother computing the whittaker function if either y1 or y2 exceed 
\\ this value
whCutOff = 15; 

read("mellinInvMethod.gp");
allocatemem();
allocatemem();
allocatemem();

\\ for lifts
mellinSetup(makerfromlift(R));

\\ for generic forms
\\ mellinSetup(r);

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ functions to compute things -- qualified personnel only \\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

gl2sum(xx,yy,mm) = {
   local(myD,t,sml=0,result=0,y2cutOff);
   
   \\determinant -- will be constant throughout
   myD = (mm[1]*yy[1])^2*mm[2]*yy[2];
   computeWeights(myD);
   
   \\ whittaker function is too small if y2 is smaller than this
   y2cutOff = myD/(whCutOff^2);

   result += cos(2*Pi*mm[2]*xx[2])*cos(2*Pi*mm[1]*xx[1])*whittakerEvaluate(mm[2]*yy[2]);
  
   c=1;
   while(sml<2,
      t = dsum(xx,yy[2],mm,c);
      result += t;
      if(abs(t/result)<eps,sml++,sml=0);
      c++;
   );
   return(result);
}

dsum(xx,yy2,mm,c) = {
   local(extgcd,zz,z,d,t1,flag,term,result);

\\   print("c: ",c);

   zz = c*(xx[2]+I*yy2);
   result = 0;
   flag=1;
 
   \\ this shoudl result in smallest possible y[1]
   d = round(-c*xx[2]);
   while(flag,
      extgcd=bezout(c,d);
      if(extgcd[3]==1,
         z = zz + d;
         t1 = mm[2]*yy2/abs(z)^2;
         if(t1>=y2cutOff,
            term = whittakerEvaluate(t1);
\\            print(t1 ," -> ", term);
\\            input();
            term *= cos(2*Pi*mm[1]*(c*xx[3]+d*xx[1]));
            term *= cos(2*Pi*mm[2]*(extgcd[2]-real(1/z))/c);
            result += term;
         ,
            flag = 0;
         );
      );
      d++;
   );

   flag = 1;
   d = round(-c*xx[2])-1;
   while(flag,
      extgcd=bezout(c,d);
      if(extgcd[3]==1,
         z = zz + d;
         t1 = mm[2]*yy2/abs(z)^2;
         if(t1>=y2cutOff,
            term = whittakerEvaluate(t1);
            term *= cos(2*Pi*mm[1]*(c*xx[3]+d*xx[1]));
            term *= cos(2*Pi*mm[2]*(extgcd[2]-real(1/z))/c);
            result += term;
         ,
            flag = 0;
         );
      );
      d--;
   );

   result;
}

\\ evaluates F(z), a GL(3) maass form
\\ xx = [x[1],x[2],x[3]]; yy = [y[1],y[2]]
\\     (1  x[1]  x[3]) (y[1]y[2]  0   0)
\\ z = (0   1    x[1]) (   0    y[1]  0)
\\     (0   0     1  ) (   0      0   0)
\\ if given a filename, will save intermediate results there to
\\ be later evaluated by evaluateFromFile
maassF(xx,yy,writeToFile = 0) = {
   local(evals=[],m1=1,m2=1,res,sml);

   until(m2 == 2,
      m2 = 1;
      evals = concat(evals,[[]]);
      res = 10^6;
      while(abs(res)>eps^0.8,
         res = gl2sum(xx,yy,[m1,m2]);
         evals[m1] = concat(evals[m1],res);
         m2++;
      );
      m1++;
   );
   
   if (writeToFile,
      write(writeToFile,"xx=",xx,";");
      write(writeToFile,"yy=",yy,";");
      write(writeToFile,"evals=",evals,";");
   );

   \\actually evaluate things
   res = 0;
   m1 = 1;
   while(m1<=length(evals),
      m2 = 1;
      while(m2<=length(evals[m1]),
         res += getCoeff(m1,m2,A)*evals[m1][m2]/(m1*m2);
         m2++;
      );
      m1++;
   );

   return(res);   
}

\\ evaluate Maass form using the info stored in a file
evaluateFromFile(filename) = {
   local(m1,m2,res,evals,xx,yy);
   read(filename);

   res = 0;
   m1 = 1;
   while(m1<=length(evals),
      m2 = 1;
      while(m2<=length(evals[m1]),
         res += getCoeff(m1,m2,A)*evals[m1][m2]/(m1*m2);
         m2++;
      );
      m1++;
   );

   return(res);   
}
