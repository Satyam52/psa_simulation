clear;clc;close all;

PH = 4.26; % high adsorption pressure in atm
PL = 1.07; % purge pressure in atm
M = 5;
[Pe,alphaA,alphaB,betaA,betaB,gamma_s,psi,vOH,L,yAF,C]=Get_Params(PH);
[Ai,Ax,Bx,z] = Get_matrix(M,Pe);

% Getting intial condition from previous step
u0 = load('finalblowdown.mat').U + 0.01;

% construct mass matrix for ODE solver
MM = zeros(4*M+8);
for j=2:M+1
    MM(j,j)=1;
end
for j=M+3:3*M+6
     MM(j,j)=1;
end

% time parmas
tmin=0;
tmax=0.2*160;
taumin = tmin*vOH/L;
taumax = tmax*vOH/L;
steps=100;
tspan = taumin:(taumax-taumin)/(steps-1):taumax;

bcfor1 = load('bcfor1.mat').bcfor1;

AbsTol = 1e-16*ones(1,4*M+8);
options = odeset('Mass',MM,'RelTol',1e-8,'AbsTol',AbsTol,"Stats",'on');


tic
[t,u] = ode15s(@(t,u) Get_Ode(t,u,PL,0,PH,M,alphaA,alphaB),tspan,u0,options);
toc


for i=2:7
    subplot(3,2,i-1);
    plot(t*L/vOH,u(:,i));
    xlabel('time') 
    ylabel('y_O_2') 
    title('At point',i);
end

