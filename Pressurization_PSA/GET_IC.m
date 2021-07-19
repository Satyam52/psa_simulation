function u0 = GET_IC(M,P,Ax)

[~,~,~,betaA,betaB,~,~,~,~,~,~] = Get_Params(P);


u0 = zeros(4*M+8,1);
u0(1:M+2) = 0.8;
u0(M+3:2*M+4) = x_star(betaA,0.8,betaB);
u0(2*M+5:3*M+6) = x_star(betaB,0.2,betaA);
u0(3*M+7:end)= 4.26/1.07;

end
