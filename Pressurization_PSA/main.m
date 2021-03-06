clear;clc;close all;

PH = 4.26*1.01325; % pressure in bar , high adsorption
PL = 1.07*1.01325; % purge pressure bar
M = 5;
[Pe,alphaA,alphaB,~,~,~,~,vOH,L,~,~] = Get_Params(PH);
[Ai,Ax,Bx,z] = Get_matrix(M,Pe);

% Intitial concentration
u0 = GET_IC(M,PL,Ax);

% construct mass matrix for ODE solver
MM = zeros(4*M+8);
for j=2:M+1
    MM(j,j)=1;
end
for j=M+3:3*M+6
     MM(j,j)=1;
end

% time params
tmin=0;
tmax=0.3*160;
taumin = tmin*vOH/L;
taumax = tmax*vOH/L;

tauP = taumax - taumin;
dPdt = (PH-PL)/tauP;

steps=100;
tspan = taumin:(taumax-taumin)/(steps-1):taumax;


AbsTol = 1e-16*ones(1,4*M+8);
options = odeset('Mass',MM,'RelTol',1e-8,'AbsTol',AbsTol,"Stats",'on');
tic
[t,u] = ode15s(@(t,u) Get_Ode(t,u,PL,dPdt,PH,M,alphaA,alphaB),tspan,u0,options);
toc
for i=1:7
    subplot(4,2,i);
    plot(t*L/vOH,u(:,21+i));
    xlabel('time') 
    ylabel('y_O_2') 
    title('At point',i);
end

