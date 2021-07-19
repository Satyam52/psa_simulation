clear;clc;close all;

PH = 4.26; % high adsorption pressure in atm
PB = 1;    % atmospheric pressure in atm
M = 5;
[Pe,alphaA,alphaB,betaA,betaB,gamma_s,psi,vOH,L,yAF,C]=Get_Params(PH);
[Ai,Ax,Bx,z] = Get_matrix(M,Pe);

% Getting initial condition from previous step
u0 = load('ufinal.mat').U;

% construct mass matrix for ODE solver
MM = zeros(4*M+8);
tmin=0;
tmax=0.3*160;
taumin = tmin*vOH/L;
taumax = tmax*vOH/L;
tauP = taumax - taumin;
dPdtau = (PB - PH)/tauP;

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
[t,u] = ode15s(@(t,u) Get_Ode(t,u,PB,dPdtau,PH,M),tspan,u0,options);
toc

for i=2:7
    subplot(3,2,i-1);
    plot(t*L/vOH,u(:,14+i));
    xlabel('time') 
    ylabel('y_O_2') 
    title('At point',i);
end
U = u(end,:);
save('finalblowdown.mat','U');
