%% Assign initial variables
Rbar = 8.314;
M = 28.97;
R = Rbar/M;
alpha = 3.653;
beta = -1.337*10^-3;
gamma = 3.294*10^-6;
rho = -1.913*10^-9;
epsilon = .2763*10^-12;
etaComp = .85;
etaTurb = .85;
etaReg = .8;
Rc = 9.9;
T1 = 288.15;
P1 = 101.325;
%% Standard equations to use below
%dsFCN = @(T) ((R((alpha*log(T/t))+(beta*(T-t))+(gamma*(T^2-t^2))+(rho*(T^3-t^3))+(epsilon*(T^4-t^4))))-(R*log(Rc)));
%dhFCN = @(T) (R*((alpha*(T-t))+(beta*(T^2-t^1))+(gamma*(T^3-t^3))+(rho*(T^4-t^4))+(epsilon*(T^5-t^5)))) - dh;
%% Solve for state 1
% solve from 0 to T1
h1 = R*((alpha*(T1))+(beta/2*(T1^2))+(gamma/3*(T1^3))+(rho/4*(T1^4))+(epsilon/5*(T1^5)));
%% Solve for state 2
% solve for ideal T2
dsFCN1i = @(T) ((R*((alpha*log(T/T1))+(beta*(T-T1))+(gamma/2*((T)^2-(T1^2)))+(rho/3*((T)^3-(T1^3)))+(epsilon/4*((T)^4-(T1^4)))))-(R*log(Rc)));
T2i = fzero(dsFCN1i,500);
% solve change in h and use efficiency to find actual change in h
dh12i = R*((alpha*(T2i-T1))+(beta/2*(T2i^2-T1^2))+(gamma/3*(T2i^3-T1^3))+(rho/4*(T2i^4-T1^4))+(epsilon/5*(T2i^5-T1^5)));
h2i = h1 + dh12i;
iWin = dh12i;
aWin = dh12i/etaComp;
dh12a = aWin;
h2a = (dh12i/etaComp) + h1;
% solve for actual T2 using the actual change in h
dhFCN1a = @(T) (R*((alpha*(T-T1))+(beta/2*(T^2-T1^2))+(gamma/3*(T^3-T1^3))+(rho/4*(T^4-T1^4))+(epsilon/5*(T^5-T1^5)))) - dh12a;
T2a = fzero(dhFCN1a,500);
ds12a = ((R*((alpha*log(T2a/T1))+(beta*(T2a-T1))+(gamma*((T2a)^2-(T1^2)))+(rho*((T2a)^3-(T1^3)))+(epsilon*((T2a)^4-(T1^4)))))-(R*log(Rc)));
%% Solving state 3 and 4 implicitly
m = 17.9; %kg/s
Power = 0;
T3 = 1340; %initial guess
while Power <= 4600
    % solve from 0 to T3
    h3 = R*((alpha*(T3))+(beta/2*(T3^2))+(gamma/3*(T3^3))+(rho/4*(T3^4))+(epsilon/5*(T3^5)));
    % solve for ideal T4
    dsFCN2i = @(T) ((R*((alpha*log(T/T3))+(beta*(T-T3))+(gamma/2*((T)^2-(T3^2)))+(rho/3*((T)^3-(T3^3)))+(epsilon/4*((T)^4-(T3^4)))))-(R*log(1/Rc)));
    T4i = fsolve(dsFCN2i,500);
    % solve change in h and use efficiency to find actual change in h
    dh34i = R*((alpha*(T3-T4i))+(beta/2*(T3^2-T4i^2))+(gamma/3*(T3^3-T4i^3))+(rho/4*(T3^4-T4i^4))+(epsilon/5*(T3^5-T4i^5)));
    h4i = h3 - dh34i;
    iWout = dh34i;
    aWout = dh34i*etaTurb;
    dh34a = aWout;
    h4a = h3 - (dh34a);
    % solve for actual T4 using the actual change in h
    dhFCN2a = @(T) (R*((alpha*(T3-T))+(beta/2*(T3^2-T^2))+(gamma/3*(T3^3-T^3))+(rho/4*(T3^4-T^4))+(epsilon/5*(T3^5-T^5)))) - dh34a;
    T4a = fzero(dhFCN2a,500);
    ds34a = ((R*((alpha*log(T4a/T3))+(beta*(T4a-T3))+(gamma*((T4a)^2-(T3^2)))+(rho*((T4a)^3-(T3^3)))+(epsilon*((T4a)^4-(T3^4)))))-(R*log(1/Rc)));
    Wnet = aWout - aWin;
    Power = m*Wnet;
    disp(T3)
    error = (Power - 4600)/4600;
    disp(error)
    disp(Power)
    T3 = T3 + .00025;
end
%% Solve for state 5 and 6

h6s = h2a;
qRegens = h4a - h6s;
h6a = h4a - (etaReg*(qRegens));
dh61a = h6a - h1;
qRegenA = h4a - h6a;
qInOld = h3 - h2a;
qInNew = qInOld - qRegenA;
qOut = dh61a;
tE = ((aWout - aWin)/qInNew) * 100;
QIn = qInNew * m;   %kJ/s
QIn_hr = QIn * 3600;    %kJ/hr





