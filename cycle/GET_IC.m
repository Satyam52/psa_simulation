function u0 = GET_IC(M,P,Ax,step_no)

[~,alphaA,alphaB,betaA,betaB,gamma_s,psi,~,~,yA,~]=Get_Params(P);

u0 = zeros(4*M+8,1);

if step_no == 4 % purge
    yAF = 0.21;
    u0(1:M+2) = yAF;
    u0(M+3:2*M+4) = x_star(betaA,yAF,betaB);
    u0(2*M+5:3*M+6) = x_star(betaB,1-yAF,betaA);
    u0(3*M+7:4*M+8) = 1;
    
elseif step_no==1 % pressurization
    
    yAF= yA;
    
    for i=1:M+2
    u0(i) = yAF;
    end

    % for xA
    for i=M+3:2*M+4
        u0(i) = x_star(betaA,yAF,betaB);
    end

    % for xB
    for i=2*M+5:3*M+6
        u0(i) = x_star(betaB,1-yAF,betaA);
    end


    % Eq. (B.12), PSA textbook, Kent
    % solve system of linear equation D*s=E, s=inv(D)*E
    % D is MxM conatning coefficients of vbar(i=2:M+1) from LHS of B.12
    % E is Mx1 matrix containing RHS of B.12
    % s is solution to  vbar(i=2:M+1), the M unknowns

    D = zeros(M);
    E = zeros(M,1);

     for j =2:M+1
         for i=1:M
         D(j-1,i)=(Ax(j,i+1)-Ax(M+2,i+1)*Ax(j,M+2)/Ax(M+2,M+2));
         end
         E(j-1)=-psi*(alphaA*(betaA*u0(j)/(1+betaA*u0(j)+betaB*(1-u0(j)))-u0(M+2+j))...
            +gamma_s*alphaB*(betaB*(1-u0(j))/(1+betaA*u0(j)+betaB*(1-u0(j)))-u0(2*M+4+j)))...
             -(Ax(j,1)-Ax(M+2,1)*Ax(j,M+2)/Ax(M+2,M+2));
     end

     s=D\E;

     u0(3*M+6+1) = 1;
     for j=2:M+1
         u0(3*M+6+j)=s(j-1);
     end

      % Eq. (B.15), PSA textbook, Kent ,Boundry value

    sum1 = 0; 
    for i=1:M+1
        sum1=sum1+Ax(M+2,i)*u0(3*M+6+i);
    end
    u0(4*M+8)=-sum1/Ax(M+2,M+2);


end


end
