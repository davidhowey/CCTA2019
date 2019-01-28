% Crossover system
% Data and Extended State Observer 
% LMI (polytopic) design

clear all
close all

%Parameters
R = 8.314472;       %[J K^-1 mol^-1]
T = 22 + 273;       %[K]
Far = 96485;        %[C/mol]
 
V_res=17.6e-3;      %[L]
c_0=0.1;            %[mol/L]
V_cell = 1.6408*1.3408*(.125*2.54)/(10^3); % volume of one half of reactor chamber in L
epsil=0.87;         %[-]
k_mt=3.3685e-6;     %[L/s] (slope = -k_mt)
E0_cell=2.2;        %[V]  (equilibrium voltage)

%Functions
I=@(t) 0*t;
dot_Nx=@(z) k_mt*c_0*(z);                  %[L/s] (slope = -k_mt)


%Data
%'t_data', 'v_data'
load crossover_data.mat

lini=100;
lend=13000;
intv=lini:20:lend;

t_data=t_data(intv);
t_data=t_data-t_data(1);
v_data=v_data(intv);

%flowrate
Q=@(t) (9e-3/60) + 0.*t;


%Space State Matrices 
A=@(t) [ 0 0; Q(t)/(epsil*V_cell), -Q(t)/(epsil*V_cell)];
B=-(1/c_0)*[1/(Far*V_res); 1/(epsil*Far*V_cell)];
E=-(1/c_0)*[1/V_res; 1/(epsil*V_cell)];
C=[0 1];



%Simulation
%---------------------------------------------------------------
%polytopic observer design: 'L_v', 'M', 'varrho'
load obs_gains_poly.mat

%---------------------------------------------------------------
%Cross over system:
cross_sys=@(t,x) [A(t)*x(1:2,:)+ E*dot_Nx(x(2,:)) + B*I(t)];

%initial condition
x0=[1;1];  %[SOC, SOC_cell]
tspan=t_data;

%simulation
[tout,xsol] = ode45(@(t,x) cross_sys(t,x),tspan,x0);

%outputs
SOC=xsol(:,1);
SOC_cell=xsol(:,2);
Vout=E0_cell+(R*T*2/Far)*log(SOC_cell./(1-SOC_cell));

%---------------------------------------------------------------
 

%---------------------------------------------------------------
%observer
ltime=length(t_data);

%Extended Model
l=1+length(M);
Psi=@(y) 0.5*(1+y);
Tinv=@(y) blkdiag(diag([1,1]),varrho*(1/Psi(y))*eye(l,l));

%orders
n0=2; %# of battery states
m=1; %# of regresors
r=1; %# of outputs

A_e=@(t,y) [A(t) E*Psi(y) zeros(n0,l-1); zeros(l-1,n0) zeros(l-1,m) M; zeros(1,n0) 0 zeros(1,l-1)]; 
E_e=[E; zeros(l,r)];
C_e=[C,zeros(r,l)];
B_e=[B; zeros(l,1)];

%observer initial condition
xe0=[0.87; 0.85];
xini=[xe0;zeros(l,1)];


%-----------------------------------------
%System: Process and Adpative Observer
obs_cross_sys=@(t,x,y) [A_e(t,x(2))*x + B_e*I(t)+ Tinv(x(2))*L_v*(y-x(2))];

                  
%Step simulation
SOCe(1)=xini(1);
SOCe_cell(1)=xini(2);
w_cell(:,1)=xini(3:end);
the_cell(1,1)=w_cell(1,1);


%Data
fact=exp((Far/(2*R*T))*(v_data-E0_cell));
%inverse_measurement
y=fact./(1+fact);

dot_Nxe(1)=Psi(y(1))*the_cell(1,1);


for k=1:ltime-1,

k
tspank=linspace(t_data(k),t_data(k+1),3);

[toute,xsole] = ode45(@(t,x) obs_cross_sys(t,x,y(k+1)),tspank,xini);

SOCe(k+1)=xsole(end,1);
SOCe_cell(k+1)=xsole(end,2);
w_cell(:,k+1)=xsole(end,3:end)';
the_cell(1,k+1)=w_cell(1,k+1);
dot_Nxe(k+1)=Psi(y(k+1))*the_cell(1,k+1);

xini(1,1)=SOCe(k+1);
xini(2,1)=SOCe_cell(k+1);
xini(3:end,1)=w_cell(:,k+1);

Vout_e=E0_cell+(R*T*2/Far)*log(SOCe_cell(1:k+1)./(1-SOCe_cell(1:k+1)));


%Figures
figure(1)

subplot(231);
plot(tout(1:k+1)/3600,SOC(1:k+1),'g--',t_data(1:k+1)/3600,SOCe,'r','LineWidth',2); 
title('SOC');
xlabel('Time[hrs]');drawnow

subplot(232);
plot(t_data(1:k+1)/3600,w_cell(:,1:k+1));
title('Parameters');
xlabel('Time[hrs]');drawnow

subplot(233);
plot(t_data(1:k+1)/3600,y(1:k+1),'g--',t_data(1:k+1)/3600,SOCe_cell,'r','LineWidth',2);  
title('SOC_{cell}');
xlabel('Time[hrs]');drawnow

subplot(223);
plot(t_data(1:k+1)/3600, dot_Nx(SOC_cell(1:k+1)),'g--',t_data(1:k+1)/3600,dot_Nxe(1:k+1),'r','LineWidth',2);
title('Cross-over Rate');
xlabel('Time[hrs]'); drawnow

subplot(224);
plot(t_data/3600, Vout,'g--',t_data/3600,v_data,'r', t_data(1:k+1)/3600,Vout_e(1:k+1),'b:','LineWidth',2);
title('Voltage: V_{out}');
xlabel('Time[hrs]'); drawnow

end


%save crossover_obs_data.mat
save crossover_obs_data_poly.mat

%------------
figure(3)
plot(tout/3600,SOC,'g-.',t_data/3600,SOCe,'r','LineWidth',2); 
title('(a)');
xlabel('Time[hrs]');
ylabel('SOC');
 
figure(4)
plot(tout(1:60)/3600,SOC(1:60),'g-.',t_data(1:60)/3600,SOCe(1:60),'r','LineWidth',2); 

figure(5)
plot(t_data/3600,y,'g-',t_data/3600,SOCe_cell,'r.','LineWidth',2);  
title('(b)');
xlabel('Time[hrs]');
ylabel('SOC_{cell}');

figure(7)
plot(t_data(1:20)/3600,y(1:20),'g-',t_data(1:20)/3600,SOCe_cell(1:20),'r.','LineWidth',2);  

figure(8)
plot(t_data/3600,the_cell);
title('(c)');
ylabel('Parameters');
xlabel('Time[hrs]');

figure(9)
plot(t_data/3600, dot_Nx(SOC_cell),'g-',t_data/3600,dot_Nxe,'r','LineWidth',2);
title('(d)');
ylabel('Cross-over Rate');
xlabel('SOC_{cell}');

figure(10)
plot(tout(1:60)/3600, dot_Nx(SOC_cell(1:60)),'g-',tout(1:60)/3600,dot_Nxe(1:60),'r','LineWidth',2);


figure(11)
plot(t_data/3600, Vout,'g',t_data/3600,v_data,'b--', t_data/3600,Vout_e,'b','LineWidth',2);
ylabel('Voltage: V_{out}');
xlabel('SOC_{cell}');
