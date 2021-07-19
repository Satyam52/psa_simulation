clear;clc;close all;

PH = 4.26; % high adsorption pressure in atm
PL=1.07; % purge pressure in atm
M = 5;
[Pe,alphaA,alphaB,betaA,betaB,gamma_s,psi,vOH,L,yAF,C]=Get_Params(PH);
[Ai,Ax,Bx,z] = Get_matrix(M,Pe);
u0 = GET_IC(M,PL,Ax);
% u0 = load('u.mat').importingU;



% construct mass matrix for ODE solver
MM = zeros(4*M+8);
tmin=0;
tmax=0.2*160;
taumin = tmin*vOH/L;
taumax = tmax*vOH/L;
steps=100;
tspan = taumin:(taumax-taumin)/(steps-1):taumax;

for j=2:M+1
    MM(j,j)=1;
end
for j=M+3:3*M+6
     MM(j,j)=1;
end

AbsTol = 1e-16*ones(1,4*M+8);
options = odeset('Mass',MM,'RelTol',1e-8,'AbsTol',AbsTol,"Stats",'on');
tic
[t,u] = ode15s(@(t,u) Get_Ode(t,u,PH,M),tspan,u0,options);
toc

% for i=1:6
%     [t,u] = ode15s(@(t,u) Get_Ode(t,u,1.45+i*.4,M),tspan,u0,options);
%     subplot(3,2,i);
%     plot(t*L/vOH,u(:,7));
%     title("Pressure",1.45+i*.4)
% end

for i=2:7
    subplot(3,2,i-1);
    plot(t*L/vOH,u(:,i));
    xlabel('time') 
    ylabel('y_O_2') 
    title('At point',i);
end
U = u(end,:);
save('ufinal.mat','U')

bcfor1 = u(:,[7 28]);
save('bcfor1.mat','bcfor1')


% plot([1,2,3,4,5,6,7],u(20,22:28));