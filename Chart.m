%% Stability and safety charts of connected cruise control
clear; clc;

%% Problem setup
xi=0.5;     % first-order lag [s]
Dst=5;      % standstill distance [m]
kappa=0.6;	% range policy gradient [1/s]
Dsf=1;      % safe distance [m]
T=1/0.6;	% safe time headway [s]
C=xi/T;     % acceleration gain [1]
vbar=15;	% maximum speed [m/s]
abar=10;    % maximum acceleration [m/s^2]

% Range of frequencies for plant and string stability boundaries
Om=0:2*pi/200:2*pi;
om=0:2*pi/200:2*pi;

% Plot settings
IA=[-0.2,1.2];	LA='$A$ (1/s)';
IB=[-0.2,1.2];	LB='$B$ (1/s)';
IC=[0,1];       LC='$C$ (1)';
Ixi=[0,2];      Lxi='$\xi$ (s)';
blue=[0,0,1];
darkgreen=[0,170,0]/256;
darkred=[230,0,0]/256;

%% Stability boundaries
% s=0 plant stability boundary
AP1=@(B)0*B;
% s=jOm plant stability boundary
AP2=@(Om)Om.^2/kappa;
BP2=@(Om)Om.^2*xi-AP2(Om);
% om=0 string stability boundaries
AS1=@(B)0*B;
AS2=@(B)2*((1-C)*kappa-B);
% om>0 string stability boundaries
AS3=@(om)om.^2*xi+(1-C)*((1+C)/2/xi-kappa)+sqrt((om.^2*xi+(1-C)*((1+C)/2/xi-kappa)).^2-om.^4*xi^2);
BS3=@(om)om.^2*xi+(1-C^2)/2/xi-AS3(om);
AS4=@(om)om.^2*xi+(1-C)*((1+C)/2/xi-kappa)-sqrt((om.^2*xi+(1-C)*((1+C)/2/xi-kappa)).^2-om.^4*xi^2);
BS4=@(om)om.^2*xi+(1-C^2)/2/xi-AS4(om);

%% Safety boundaries
A1=@(B)(abs(1/T-xi/T^2-B)*vbar+abs(xi/T-C)*abar)/kappa/(Dst-Dsf);
A2=@(B)(1-xi/T)^2/4/xi+0*B;
C0=@(xi)xi/T;
C1=@(xi)xi/T - kappa*(Dst-Dsf)*(1-xi/T).^2/4./xi/abar;
C2=@(xi)xi/T + kappa*(Dst-Dsf)*(1-xi/T).^2/4./xi/abar;
xi1=@(C)T+0*C;

%% Stability and safety chart
% plot chart
figure(2); clf; hold on; box on;
BB=linspace(IB(1),IB(2),1001);
% plot s=0 plant stability boundary
plot(IB,AP1(IB),'--','Color',darkred,'LineWidth',2,'HandleVisibility','off');
% plot s=jOm plant stability boundary
plot(BP2(Om),AP2(Om),'Color',darkred,'LineWidth',2,'HandleVisibility','off');
% plot om=0 string stability boundaries
plot(IB,AS1(IB),'--','Color',blue,'LineWidth',2,'HandleVisibility','off');
plot(IB,AS2(IB),'--','Color',blue,'LineWidth',2,'HandleVisibility','off');
% plot om>0 string stability boundaries
plot(BS3(om),AS3(om),'Color',blue,'LineWidth',2,'HandleVisibility','off');
plot(BS4(om),AS4(om),'Color',blue,'LineWidth',2,'HandleVisibility','off');
% plot safety boundaries
plot(BB,A1(BB),'Color',darkgreen,'LineWidth',2,'HandleVisibility','off');
plot(IB,A2(IB),'Color',darkgreen,'LineWidth',2,'HandleVisibility','off');
PlotFinalize({LB,LA},[IB,IA]);

% plot critical lag
figure(3); clf; hold on; box on;
CC=linspace(IC(1),IC(2),1001);
XI=linspace(Ixi(1),Ixi(2),1001);
plot(C0(Ixi),Ixi,'--','Color',darkgreen,'LineWidth',1,'HandleVisibility','off');
plot(IC,xi1(IC),'--','Color',darkgreen,'LineWidth',1,'HandleVisibility','off');
plot(C1(XI),XI,'Color',darkgreen,'LineWidth',2,'HandleVisibility','off');
plot(C2(XI),XI,'Color',darkgreen,'LineWidth',2,'HandleVisibility','off');
PlotFinalize({LC,Lxi},[IC,Ixi]);

%% Save results
filename=['chart','_xi',num2str(xi,'%3.2f'),'_C',num2str(C,'%3.2f')];
filename=strrep(filename,'.','');
% figure(2);
% saveas(gcf,[filename,'.fig']);
% saveas(gcf,[filename,'.svg']);