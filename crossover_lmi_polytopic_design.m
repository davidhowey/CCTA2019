% Crossover system, Extended State Observer
% Solution LMI Polytopic Problem
% Mosek and Yalmip needed


clear all
close all

ops = sdpsettings('verbose',1,'savesolveroutput',1,'showprogress',0,...
    'solver','mosek-sdp',...
     'mosek.MSK_DPAR_INTPNT_TOL_PFEAS',1e-11,...
     'mosek.MSK_DPAR_INTPNT_TOL_DFEAS',1e-11,...
     'mosek.MSK_DPAR_INTPNT_TOL_REL_GAP',1e-11);

%Parameters
R = 8.314472;       %[J K^-1 mol^-1]
T = 22 + 273;       %[K]
Far = 96485;        %[C/mol]
 
V_res=17.6e-3;      %[L]
c_0=0.1;            %[mol/L]
dot_V= 9e-3/60;     %[L/s]
V_cell = 1.6408*1.3408*(.125*2.54)/(10^3); % volume of one half of reactor chamber in L
epsil=0.87;         %[-]
k_mt=3.3685e-6;     %[L/s] (slope = -k_mt)
E0_cell=2.2;        %[V]  (equilibrium voltage)

%Functions
I=@(t) 0*t;
dot_Nx=@(z) k_mt*c_0*(z);                  %[L/s] (slope = -k_mt)


%Space State Matrices 
Acell=[ 0 0; 1/(epsil*V_cell), -1/(epsil*V_cell)];
B=-(1/c_0)*[1/(Far*V_res); 1/(epsil*Far*V_cell)];
E=-(1/c_0)*[1/V_res; 1/(epsil*V_cell)];
C=[0 1];

% dot_V_min=10e-3/60; %[L/s]
% dot_V_max=20e-3/60; %[L/s]

dot_V_min=9e-3/60;
dot_V_max=9e-3/60; %[L/s]

Qdom=[0.25*dot_V_min, 2*dot_V_max];

%Adaptive Observer Design
[n0,n0]=size(Acell);
[n0,m]=size(E);
[r,n0]=size(C);

%Extended Model

% l=4;
% lambda=[0.25 0.025 0.0025];

%case 1
l=3;
lambda=[0.5 0.05];

% %case 2
% l=2;
% lambda=[5];

M=diag(lambda);
varrho=1e-4;

A=@(q) [q*Acell varrho*E zeros(n0,l-1); zeros(l-1,n0) zeros(l-1,m) M; zeros(1,n0) 0 zeros(1,l-1)]; 
E_e=[E; zeros(l,r)];
C_e=[C,zeros(r,l)];
B_e=[B; zeros(l,1)];
Iv= [eye(n0,n0); zeros(l,n0)];
Ibar=Iv*Iv';

%order
n=l+n0;

P = sdpvar(n,n,'symmetric');
W = sdpvar(n,n,'symmetric');
Z= sdpvar(n,r);
balpha=sdpvar(1,1);
gamma_z=sdpvar(1,1);
gamma_f=sdpvar(1,1);
In=eye(n);
beta=1e-4;

CT=[...
    P>=0,...
    W>=blkdiag([1e-3, 0; 0 1e-6], 1e-3*eye(l,l)),...
    [gamma_z*In  Z; Z' gamma_z]>=0,...
    gamma_z>=0,...
    balpha>=1e-3,...
    ];

for i=1:length(Qdom)
    
       CT=[CT,...
           [-A(Qdom(i))'*P-P*A(Qdom(i))+C_e'*Z'+Z*C_e-(beta*Ibar+W), P; ...
           P, balpha*In]>=0,...
           ];
           
end

kappa_z=0.01;

sol=optimize(CT,[balpha+kappa_z*gamma_z],ops);
P_v=value(P);
W_v=value(W);
Z_v=value(Z);
balpha_v=value(balpha)

eig(P_v)
eig(W_v)
An=A(Qdom(2));
eig(An'*P_v+P_v*An-C_e'*Z_v'-Z_v*C_e+ (1/balpha_v)*P_v*Ibar*P_v+(beta*Ibar+W_v))
L_v=inv(P_v)*Z_v



save('obs_gains_poly', 'L_v', 'M', 'varrho');




 
