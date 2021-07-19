function u0 = Get_IC(P,M)
u0 = zeros(4*M,1);
[Pe,alphaA,alphaB,betaA,betaB,gamma_s,psi,vOH,L,yAF,C] = Get_Params(P);
x_Astar0 = x_star(betaA,yAF,betaB);
x_Bstar0 = x_star(betaB,1-yAF,betaA);

for i = 1:4*M
   if i <= M
       u0(i) = 0.21; % yA = 0.21
   elseif i <= 2*M
        u0(i) = x_Astar0; 
   elseif i <= 3*M
        u0(i) = x_Bstar0;
   elseif i<=4*M
       u0(i) = 1;  % vbar = 1
   end
end

end