\\ find a fundamental domain equivalent of Z,
\\ Z = [X,Y], where X=[x1,x2,x3],Y=[y1,y2] from
\\ the Iwasawa form

\\ translation by 1 in all x's
gl3T1(Z) = {Z[1][1]+=1;Z}
gl3T3(Z) = {Z[1][3]+=1;Z}
gl3T1i(Z) = {Z[1][1]-=1;Z}
gl3T3i(Z) = {Z[1][3]-=1;Z}
gl3T2i(Z) = {Z[1][2]-=1;Z[1][3]-=Z[1][1];Z}
gl3T2(Z) = {Z[1][2]+=1;Z[1][3]+=Z[1][1];Z}

\\ make x's between -0.5 and 0.5 
gl3T(Z) = {
   local(xrnd);
   xrnd = round(Z.getX[2]);
   Z[1][3] -= xrnd*Z[1][1];
   xrnd = round(Z.getX);
   Z[1] -= xrnd;
   Z
}

\\ flip in y2
gl3S2(Z) = {
   local(sq);
   sq = Z[1][2]^2+Z[2][2]^2;
   [[Z[1][3],-Z[1][2]/sq,-Z[1][1]],[Z[2][1]*sqrt(sq),Z[2][2]/sq]]
}

\\ flip in y1
gl3S1(Z) = {
   local(sq,myx2,myx1);
   sq = Z[1][1]^2+Z[2][1]^2;
   myx2 = Z[1][3]-Z[1][1]*Z[1][2];
   myx1 = - Z[1][1]/sq;
   [[myx1,myx2,myx1*myx2-Z[1][2]],[Z[2][1]/sq,Z[2][2]*sqrt(sq)]]
}

\\ create Z from list of x's and y's vectors
makeZ(x1,x2,x3,y1,y2) = [[x1,x2,x3],[y1,y2]];

Z.getX = Z[1];
Z.getY = Z[2];

\\ find an equivalent point with maximum possible y1 y2
flipBack(Z) = {
   local(cnt);

   if(debug,print(Z));

   \\ don't want to go on forever
   cnt = 0;

   while(cnt<10,
      if(Z.getY[1]<1,
         Z = gl3S1(gl3T(Z));
         cnt++;
         ,
         if(Z.getY[2]<1,
            Z = gl3S2(gl3T(Z));
            cnt++;
            ,
            return(gl3T(Z)) \\ both greater than 1
         );
      );
   );
   return(gl3T(Z));   
}
