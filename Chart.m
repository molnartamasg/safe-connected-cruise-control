%% Stability and safety charts of connected cruise control
clear; clc;

%% Problem setup
Dst=5;      % standstill distance [m]
kappa=0.6;	% range policy gradient [1/s]
C=0;        % acceleration gain [1]
CBF='TH';	% choice of CBF, 'TH' (with C=0) or 'TT'
Dsf=1;      % safe distance [m]
T=1/0.6;	% safe time headway / time to conflict [s]
vbar=15;	% maximum speed [m/s]
gamma=@(vL)sqrt(20*vL); % acceleration bound for TTC, aL>=-gamma(vL)

% Range of frequencies for plant and string stability boundaries
Om=0:2*pi/200:2*pi;

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
BP2=@(Om)-AP2(Om);
% om=0 string stability boundaries
AS1=@(B)0*B;
AS2=@(B)2*((1-C)*kappa-B);

%% Safety boundaries
switch CBF
    case 'TH'
        A=@(B)abs(1/T-B)*vbar/kappa/(Dst-Dsf);
        B=linspace(IB(1),IB(2),201);
    case 'TTC'
        AA=linspace(IA(1),IA(2),201);
        BB=linspace(IB(1),IB(2),201);
        VVL=linspace(0,vbar,201);
        [A,B,vL]=meshgrid(AA,BB,VVL);
        S=min((1/T-B+A).*vL-(1-C)*gamma(vL) + min(0,B-1/T)*vbar + A*kappa*(Dst-Dsf),[],3);
end

%% Stability and safety chart
% plot chart
figure(2); clf; hold on; box on;
% plot s=0 plant stability boundary
plot(IB,AP1(IB),'Color',darkred,'LineWidth',2,'HandleVisibility','off');
% plot s=jOm plant stability boundary
plot(BP2(Om),AP2(Om),'Color',darkred,'LineWidth',2,'HandleVisibility','off');
% plot om=0 string stability boundaries
plot(IB,AS1(IB),'Color',blue,'LineWidth',2,'HandleVisibility','off');
plot(IB,AS2(IB),'Color',blue,'LineWidth',2,'HandleVisibility','off');
% plot safety boundaries
switch CBF
    case 'TH'
        plot(B,A(B),'Color',darkgreen,'LineWidth',2,'HandleVisibility','off');
    case 'TTC'
        contour(B(:,:,1),A(:,:,1),S,[0,0],...
                    'Color',darkgreen,'LineWidth',2,'HandleVisibility','off');
end
PlotFinalize({LB,LA},[IB,IA]);

%% Save results
filename=['chart_',CBF,'_vbar',num2str(vbar),'_C',num2str(C,'%3.2f')];
filename=strrep(filename,'.','');
% figure(2);
% saveas(gcf,[filename,'.fig']);
% saveas(gcf,[filename,'.svg']);