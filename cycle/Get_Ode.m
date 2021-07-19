function dUdt = Get_Ode(t,u,PL,dPdt,PH,M)

P=PH;

% Getting Parameters
[Pe,alphaA,alphaB,betaA,betaB,gamma_s,psi,~,~,yAF,~]=Get_Params(PH);
[Ai,Ax,Bx,~] = Get_matrix(M,Pe);
Pm = 1/Pe;
dPdt = 0;

% u(1)...u(M+2) is yA
% u(M+3)...u(2*(M+2)) is xA
% u(2M+5)..u(3*(M+2)) is xB
% u(3M+7)..u(4*(M+2)) is vbar

%% for yA at z=0, i=1 ,
dUdt = zeros(8*(M+2),1);
term1 = 0;
term2 = 0;
for i=2:M+1
    term1 = term1 + (Ai(3)*Ax(M+2,i) - Ax(1,i))*u(i);
    term2 = term2 + (Ax(M+2,i)*u(i) );  % eS
end
dUdt(1) = - Ai(5)*term1-Ai(4)*term2 + Ai(5)*Pe*u(3*(M+2)+1)*yAF - u(1); % algebric equation

% for yA at i=2...M+1
for j=2:M+1
term1 = 0;
term2 = 0;
term3 = 0;
for i=2:M+1
    term1 = term1 + (Pm*Bx(j,i)-u(3*M+6+j)*Ax(j,i))*u(i);
    term2 = term2 + (Ai(3)*Ax(M+2,i)-Ax(1,i))*u(i);
    term3 = term3 + Ax(M+2,i)*u(i);
end
dUdt(j) = term1-Ai(5)*(Pm*Bx(j,1)-u(3*M+6+j)*Ax(j,1))*term2...
         +Ai(1)*(Pm*Bx(j,M+2)-u(3*M+6+j)*Ax(j,M+2))*term2...
         -Ai(4)*(Pm*Bx(j,1)-u(3*M+6+j)*Ax(j,1))*term3...
         +Ai(5)*(Pm*Bx(j,1)-u(3*M+6+j)*Ax(j,1))*Pe*yAF*u(3*M+6+1)...
         -Ai(1)*(Pm*Bx(j,M+2)-u(3*M+6+j)*Ax(j,M+2))*Pe*yAF*u(3*M+6+1)...
         +psi*((u(j)-1)*alphaA*(betaA*u(j)/(1+betaA*u(j)+betaB*(1-u(j)))-u(M+2+j))...
         +gamma_s*u(j)*alphaB*(betaB*(1-u(j))/(1+betaA*u(j)+betaB*(1-u(j)))-u(2*M+4+j)));  
end

% boundry condition, i=M+2
term1 = 0;
for i=2:M+1
term1 = term1 + (Ai(3)*Ax(M+2,i)- Ax(1,i))*u(i);
end
dUdt(M+2) = Ai(1)*term1 - Ai(1)*Pe*yAF*u(3*M+6+1) - u(M+2);  % Ai(1)*Pe*yAF*u(3*M+6+1) maybe correction


%% for xA, xB
% Eq for j=1 and j=M+2 are not in book, but should work
for j=1:M+2
     dUdt(M+2+j)=alphaA*(betaA*u(j)/(1+betaA*u(j)+betaB*(1-u(j)))-u(M+2+j));
     dUdt(2*M+4+j)=alphaB*(betaB*(1-u(j))/(1+betaA*u(j)+betaB*(1-u(j)))-u(2*M+4+j));
end

%% for velocity
% velocity at 1 is vOH so vbar = 1

dUdt(3*M+6+1) = u(3*M+6+1)-1;

 % Eq. (B.12), PSA textbook, Kent

 for j =2:M+1
     term1 = 0;
     for i=2:M+1
         term1=term1+(Ax(j,i)-Ax(M+2,i)*Ax(j,M+2)/Ax(M+2,M+2))*u(3*M+6+i);
     end
     dUdt(3*M+6+j)= -term1-psi*(alphaA*(betaA*u(j)/(1+betaA*u(j)+betaB*(1-u(j)))-u(M+2+j))...
                   +gamma_s*alphaB*(betaB*(1-u(j))/(1+betaA*u(j)+betaB*(1-u(j)))-u(2*M+4+j)))...
                   -(Ax(j,1)-Ax(M+2,1)*Ax(j,M+2)/Ax(M+2,M+2))*u(3*M+6+1)...
                   -1/P*dPdt;
 end


% boundry condition at z=1, i=M+2
sum1 = 0; 
for i=2:M+1
    sum1=sum1+Ax(M+2,i)*u(3*M+6+i);
end
dUdt(4*M+8)=-Ax(M+2,1)/Ax(M+2,M+2)*u(3*M+6+1)-sum1/Ax(M+2,M+2)-u(4*M+8);






%% For Bed 1
[~,~,~,betaA,betaB,~,psi,~,~,~,~]=Get_Params(PL);

% u(4M+8+1)...u(5M+10) is yA
% u(5M+11)...u(6*(M+2)) is xA
% u(6M+13)..u(7*(M+2)) is xB
% u(7M+15)..u(8*(M+2)) is vbar

%% for yA at z=0, i=1 
term1 = 0;
term2 = 0;
for i=2:M+1
    term1 = term1 + (Ai(3)*Ax(M+2,i) - Ax(1,i))*u(4*M+8+i);
    term2 = term2 + (Ax(M+2,i)*u(4*M+8+i));
end
dUdt(4*M+8+1) = -Ai(5)*term1 -Ai(4)*term2 + Ai(5)*Pe*u(7*M+14+1)*u(M+2) - u(4*M+8+1); % algebric equation

% for yA at i=2...M+1
for j=2:M+1
term1 = 0;
term2 = 0;
term3 = 0;
for i=2:M+1
    term1 = term1 + (Pm*Bx(j,i)-u(7*M+14+j)*Ax(j,i))*u(4*M+8+i);
    term2 = term2 + (Ai(3)*Ax(M+2,i)-Ax(1,i))*u(4*M+8+i);
    term3 = term3 + Ax(M+2,i)*u(4*M+8+i);
end
dUdt(4*M+8+j) = term1-Ai(5)*(Pm*Bx(j,1)-u(7*M+14+j)*Ax(j,1))*term2...
         +Ai(1)*(Pm*Bx(j,M+2)-u(7*M+14+j)*Ax(j,M+2))*term2...
         -Ai(4)*(Pm*Bx(j,1)-u(7*M+14+j)*Ax(j,1))*term3...
         +Ai(5)*(Pm*Bx(j,1)-u(7*M+14+j)*Ax(j,1))*Pe*u(7*M+14+1)*u(M+2)...
         -Ai(1)*(Pm*Bx(j,M+2)-u(7*M+14+j)*Ax(j,M+2))*Pe*u(7*M+14+1)*u(M+2)...
         +psi*((u(4*M+8+j)-1)*alphaA*(betaA*u(4*M+8+j)/(1+betaA*u(4*M+8+j)+betaB*(1-u(4*M+8+j)))-u(5*M+10+j))...
         +gamma_s*u(4*M+8+j)*alphaB*(betaB*(1-u(4*M+8+j))/(1+betaA*u(4*M+8+j)+betaB*(1-u(4*M+8+j)))-u(6*M+12+j)));  
end

% boundry condition, i=M+2
term1 = 0;
for i=2:M+1
term1 = term1 + (Ai(3)*Ax(M+2,i)- Ax(1,i))*u(4*M+8+i);
end
dUdt(5*M+10) = Ai(1)*term1 - Ai(1)*Pe*u(7*M+14+1)*u(M+2) - u(5*M+10);  


%% for xA, xB
% Eq for j=1 and j=M+2 are not in book, but should work
for j=1:M+2
     dUdt(5*M+10+j)=alphaA*(betaA*u(4*M+8+j)/(1+betaA*u(4*M+8+j)+betaB*(1-u(4*M+8+j)))-u(5*M+10+j));
     dUdt(6*M+12+j)=alphaB*(betaB*(1-u(4*M+8+j))/(1+betaA*u(4*M+8+j)+betaB*(1-u(4*M+8+j)))-u(6*M+12+j));
end


%% for velocity
% velocity at 1 is vOH so vbar = vbar exit of bed 2

dUdt(7*M+14+1) = u(7*M+14+1) - u(4*M+8);

 % Eq. (B.12), PSA textbook, Kent

 for j =2:M+1
     term1 = 0;
     for i=2:M+1
         term1=term1+(Ax(j,i)-Ax(M+2,i)*Ax(j,M+2)/Ax(M+2,M+2))*u(7*M+14+i);
     end
     dUdt(7*M+14+j)= -term1 -psi*(alphaA*(betaA*u(4*M+8+j)/(1+betaA*u(4*M+8+j)+betaB*(1-u(4*M+8+j)))-u(5*M+10+j))...
                   +gamma_s*alphaB*(betaB*(1-u(4*M+8+j))/(1+betaA*u(4*M+8+j)+betaB*(1-u(4*M+8+j)))-u(6*M+12+j)))...
                   -(Ax(j,1)-Ax(M+2,1)*Ax(j,M+2)/Ax(M+2,M+2))*u(7*M+14+1)...
                   -1/P*dPdt;
 end


% boundry condition at z=1, i=M+2
sum1 = 0; 
for i=2:M+1
    sum1=sum1+Ax(M+2,i)*u(7*M+14+i);
end
dUdt(8*M+16)=-Ax(M+2,1)/Ax(M+2,M+2)*u(7*M+14+1)-sum1/Ax(M+2,M+2)-u(8*M+16);


end