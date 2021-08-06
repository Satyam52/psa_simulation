clear;clc;close all;

PH = 4.26*1.01325;  % high adsorption pressure in atm
PL = 1.07*1.01325;  % purge pressure in atm
PB = 1;     
M = 5;
[Pe,alphaA,alphaB,betaA,betaB,gamma_s,psi,vOH,L,yAF,C]=Get_Params(PH);
[Ai,Ax,Bx,z] = Get_matrix(M,Pe);

% construct mass matrix for ODE solver
MM = zeros(8*M+16);
for j=2:M+1
    MM(j,j)=1;
end
for j=M+3:3*M+6
     MM(j,j)=1;
end

for j=4*M+10:5*M+9
    MM(j,j)=1;
end

for j=5*M+11:7*M+14
     MM(j,j)=1;
end


% set option for ode15s
AbsTol = 1e-16*ones(1,8*M+16);
options = odeset('Mass',MM,'RelTol',1e-8,'AbsTol',AbsTol,"Stats",'on');

% Writing for a cycle
% Get_Ode(t,u,PL,dPdt,PH,M,step_no,bed_no)
% 1. Bed2---> Pressurization and Bed1--->Blowdown


% 1. Bed2---> High pressure adsorption and Bed1--->Purge

    % bed2 - high P adsorption
    u0=[GET_IC(M,PL,Ax,1);GET_IC(M,PL,Ax,4)];
    tmin=0;
    tmax=0.3*160;
    taumin = tmin*vOH/L;
    taumax = tmax*vOH/L;
    steps=100;
    tspan = taumin:(taumax-taumin)/(steps-1):taumax;
    tic
    [t,u] = ode15s(@(t,u) Get_Ode(t,u,PL,0,PH,M),tspan,u0,options);
    toc
%     uexit = u(end,:);
%     save('ads_exit.mat','uexit');
    
    % for bed1 - purge
%     tic
%     [t,u] = ode15s(@(t,u) Get_Ode(t,u,PL,0,PH,M,4,1),tspan,u0,options);
%     toc

        


% for i=2:7
%     subplot(3,2,i-1);
%     plot(t*L/vOH,u(:,i)*100);
%     xlabel('time') 
%     ylabel('y_O_2') 
%     title('At point',i);
% end
% U = u(end,:);
% save('ufinal.mat','U')
