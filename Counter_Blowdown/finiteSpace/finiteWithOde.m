clear;clc;close all;

P = 1.48; % atm
PL = 1.07;
[Pe,alphaA,alphaB,betaA,betaB,gamma_s,psi,vOH,L,yAF,C] = Get_Params(P);
x_Astar0 = x_star(betaA,yAF,betaB);
x_Bstar0 = x_star(betaB,1-yAF,betaA);
M = 50;
dz =1/(M-1);


% intial conditions at t = 0
% yA = 0.21 , yB = 0.79 and xA = x_Astar, vbar = 1
% 1 to M --- > yA
% M+1 to 2M ---> xA
% 2M+1 to 3M ---> xB
% 3M+1 to 4M ---> vbar
% 
U0 = Get_IC(PL,M);

MM = zeros(4*M);
for i=2:M-1
    MM(i,i)=1;
end 

for i=M+1:3*M
  MM(i,i) = 1;
end

tmin = 0;
tmax = 1000; % in seconds
taumin = tmin*vOH/L;
taumax = tmax*vOH/L;
tspan = linspace(taumin,taumax,100);

AbsTol = 1e-16*ones(1,4*M);
options = odeset('Mass',MM,'Stats','on','RelTol',1e-8,'AbsTol',AbsTol);

% solve DAE
tic
[t,U] = ode15s(@(t,U) InitializeOde(t,U,P,M),tspan,U0,options);
toc

plot(t*L/vOH,U(:,M))


% 1 to M --- > yA
% M+1 to 2M ---> xA
% 2M+1 to 3M ---> xB
% 3M+1 to 4M ---> vbar
function dUdt = InitializeOde(t,U,P,M)
[Pe,alphaA,alphaB,betaA,betaB,gamma_s,psi,vOH,L,yAF,C] = Get_Params(P);
dz =1/(M-1);
Pm = 1/Pe;
dUdt = zeros(1,4*M);

% usually dU(i)/dt, 
% 1.For yA
dUdt(1) = U(1) - 0.21 + dz*Pe*U(3*M+1)*(0.21-U(1)) ;  %% boundry condition @ t = 0
for i=2:M-1
    dUdt(i) =  Pm*(U(i-1)-2*U(i)+U(i+1))/(dz*dz) - U(3*M+i)*(U(i) - U(i-1))/(dz) + ...
               psi*(alphaA*(U(i)-1)*(betaA*U(i)/(1+betaA*U(i)+betaB*(1-U(i)))-U(M+i)) + ...
               gamma_s*alphaB*U(i)*(betaB*(1-U(i))/(1+betaA*U(i)+betaB*(1-U(i))) - U(2*M+i)));
end
dUdt(M) = U(M)-U(M-1);


% 2.For xA
for i=M+1:2*M
   dUdt(i) = alphaA*(betaA*U(i-M)/(1+betaA*U(i-M)+betaB*(1-U(i-M)))-U(i));
end

% 3.For xB
for i=2*M+1:3*M
   dUdt(i) = alphaB*(betaB*(1-U(i-2*M))/(1+betaA*U(i-2*M)+betaB*(1-U(i-2*M)))-U(i));
end

% 4.For vbar, only algebric
dUdt(3*M+1) = U(3*M+1) - 1;
for i=3*M+2:4*M-1
    dUdt(i) = U(i-1)-psi*dz*(alphaA*(1-U(i-3*M))*(betaA*U(i-3*M)/(1+betaA*U(i-3*M)+betaB*(1-U(i-3*M)))-U(i-2*M)) ...
              +gamma_s*alphaB*U(i-3*M)*(betaB*(1-U(i-3*M))/(1+betaA*U(i-3*M)+betaB*(1-U(i-3*M)))-U(i-M)))...
              - U(i);
end
dUdt(4*M) = U(4*M) - U(4*M-1);

dUdt = dUdt';
end
















