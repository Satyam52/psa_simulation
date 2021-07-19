function [Pe,alphaA,alphaB,betaA,betaB,gamma_s,psi,vOH,L,yAF,C]=Get_Params(P)


% Parameter values refer Farooq, 1989, Chem Eng. Sci, pp 2809
Pe = 500;
L = 35;         % bed length in cm
r = 1.75;       % column radius in cm
A = pi()*r^2;   % column area in cm2
eps=0.4;        % bed porosity  
P = P*1.01325;
T = 25+273.15;     % column temperature in K
Q = 66.7;            % feed flow rate in cm3/s at 1 atm and 25 C
QP = Q*1.01325/P;  % feed flow rate in cm3/s at column T and P
vOH = QP/A;        % interstial velocity at column inlet and P
R = 83.14;         % bar cc / mol/K
qAs = 5.26e-3;     % saturation capacity for component A,O2, mol/cc
qBs = 5.26e-3;     % saturation capacity for component B,N2, mol/cc
kA = 62;           % 1/s LDF rate constant for A, O2
kB = 19.7;         % 1/s LDF rate constant for B, N2
KA = 4.7;          % Langmuir constant for A, O2
KB = 14.8;         % Langmuir constant for B, N2
yAF = 0.21;        % yA in feed, O2


alphaA = L*kA/vOH; % dimensionless rate constant in Eq B.4, PSA textbook, Kent
alphaB = L*kB/vOH; % dimensionless rate constant in Eq B.4, PSA textbook, Kent
bA = KA/qAs;       % Langmuir parameter in cc/mol for A, O2
bB = KB/qBs;       % Langmuir parameter in cc/mol for B, N2
C = P/R/T;         % gas phase concentration in mol/cc at PH and T
betaA = bA*C;      % dimensionless constant in Eq B.4, PSA textbook, Kent
betaB = bB*C;      % dimensionless constant in Eq B.4, PSA textbook, Kent
gamma_s= qBs/qAs;  % dimensionless ratio of saturation capacities
psi = (1-eps)/eps*R*T/P*qAs;      % dimensionless constant in Eq B.5, PSA textbook, Kent

end