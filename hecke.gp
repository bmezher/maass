\\ computes A(m1,m2) given a vector of A(m1,1)'s
\\ horrible code... but it's 2 am....

getCoeff(m1,m2,vecAs) = {
   local(d);
   if(m2==1,return(vecAs[m1]));
   if(m1==1,return(conj(vecAs[m2])));
   d = gcd(m1,m2);
   return(vecAs[m1]*conj(vecAs[m2])-if(d>1,sumdiv(d,i,if(i>1,getCoeff(m1/i,m2/i,vecAs),0)),0));
}


\\ return a vector of A(m,1)'s for a given vector of A(p,1)'s
getVecAs(primeAs) = {
   local(n,a_1,a,a1,a2,result,p,ppow);
   \\ how long should the vector be?
   n = prime(length(primeAs));
   until(isprime(n),n++);
    
   result = vector(n-1,i,1);

   for(i=1,length(primeAs),
      a_1=0; a=0; a2=primeAs[i]; a1 = 1; p = prime(i); ppow = p;
     while(ppow < n,
         forstep(j=ppow,n-1,ppow,result[j] *= a2/a1);
         ppow *= p;
         a_1 = a; a = a1; a1 = a2; a2 = a1*primeAs[i]-a*conj(primeAs[i])+a_1;
      ); 
   );

   return(result);
}
