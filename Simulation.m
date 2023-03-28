%% Simulation and safety evaluation of connected cruise control
clear; clc;

%% Problem setup
% Setup and parameters
par.p=@(v)0*v.^2;           % resistance term [m/s^2]	% 0.1+3e-4*v.^2
% desired controller
par.Dst=5;                  % standstill distance [m]
par.kappa=0.6;              % range policy gradient [1/s]
par.vmax=15;                % speed limit [m/s]
par.A=0.4;                  % distance related gain [1/s]
par.B=0.3;                  % velocity difference gain [1/s]
par.C=0;                    % acceleration gain [1]
% safety
par.filter='off';           % safety filter 'on', 'off'
par.Dsf=1;                  % safe distance [m]
par.T=1/0.6;                % safe time headway [s]
par.alpha=@(h)h;            % class-K function (xi=0)

% Simulation settings
t0=0;                       % start time [s]
tend=20;                    % end time [s]
dt=0.01;                    % time step [s]
t=t0:dt:tend;               % time [s]

% Initial conditions
v0=15;                      % CAV's velocity
vL0=v0;                     % leader's velocity
D0=v0/par.kappa+par.Dst;    % distance
x0=[D0;v0;vL0];             % state

% Plot settings
It=[0,10];  Lt='time, $t$ (s)';
ID=[0,35];  LD='distance, $D$ (m)';
Iv=[0,15];  Lv='speed, $v$ (m/s)';
Ia=[-10,2]; La='acceleration, $a$ (m/s$^2$)';
Ih=[-2,4];	Lh='CBF, $h$ (m/s)';
Iu=[-10,2]; Lu='input, $u$ (m/s$^2$)';
purple=[170,0,170]/256;
orange=[255,170,0]/256;
black=[0,0,0];
blue=[0,0,1];
darkgreen=[0,170,0]/256;
darkred=[230,0,0]/256;

%% Simulation
% Simulate system
sol=ode45(@(t,x)rhs(t,x,par),[t0,tend],x0,odeset('RelTol',1e-6));
x=deval(sol,t);

% Evaluate solution
[D,v,vL]=states(x);
a=gradient(v)./gradient(t);
aL=accL(t);
[u,ud,h]=k(t,x,par);

%% Plot results
figure(1); clf;

% Distance
subplot(2,3,1); hold on; box on;
plot(t,D,'Color',blue,'LineWidth',2,'HandleVisibility','off');
PlotFinalize({Lt,LD},[It,ID]);

% Velocity
subplot(2,3,2); hold on; box on;
plot(t,vL,'Color',darkred,'LineWidth',2,'DisplayName','$v_{\rm L}$');
plot(t,v,'Color',blue,'LineWidth',2,'DisplayName','$v$');
PlotFinalize({Lt,Lv},[It,Iv]);

% Acceleration
subplot(2,3,3); hold on; box on;
plot(t,aL,'Color',darkred,'LineWidth',2,'DisplayName','$a_{\rm L}$');
plot(t,a,'Color',blue,'LineWidth',2,'DisplayName','$a$');
PlotFinalize({Lt,La},[It,Ia]);

% Distance-velocity plot
subplot(2,3,4); hold on; box on;
DD=linspace(ID(1),ID(2),101);
Vdes=min(par.kappa*(DD-par.Dst),par.vmax);	% range policy
Vsafe=(DD-par.Dsf)/par.T;                   % safe set boundary
plot(DD,Vdes,'--','Color',black,'LineWidth',1,'DisplayName','$V$')
plot(DD,Vsafe,'Color',darkgreen,'LineWidth',2,'DisplayName','$\partial S$')
plot(D,v,'Color',blue,'LineWidth',2,'DisplayName','$v$')
PlotFinalize({LD,Lv},[ID,Iv]);

% CBF
subplot(2,3,5); hold on; box on;
plot(It,[0,0],'k','LineWidth',1,'HandleVisibility','off');
plot(t,h,'Color',blue,'LineWidth',2,'HandleVisibility','off');
PlotFinalize({Lt,Lh},[It,Ih]);

% Control input
subplot(2,3,6); hold on; box on;
plot(t,ud,'Color',orange,'LineWidth',2,'DisplayName','$k_{\rm d}(x)$');
plot(t,u,'Color',blue,'LineWidth',2,'DisplayName','$k(x)$');
PlotFinalize({Lt,Lu},[It,Iu]);

%% Save results
filename=['sim_',par.filter,'_A',num2str(par.A,'%3.2f'),...
          '_B',num2str(par.B,'%3.2f'),'_C',num2str(par.C,'%3.2f')];
filename=strrep(filename,'.','');
% saveas(gcf,[filename,'.fig']);
% saveas(gcf,[filename,'.svg']);

%% Functions for dynamics
% States
function [D,v,vL] = states(x)
    D = x(1,:);
    v = x(2,:);
    vL = x(end,:);
end

% Lead vehicle's acceleration
function aL = accL(t)
    aL = -10*(t-3).*(3<=t & t<=4)...
         -10*(4<t & t<=4.5)...
         +(10*(t-4.5)-10).*(4.5<t & t<=5.5);
end

% System model
function [f,g] = sys(t,x,par)
    [~,v,vL] = states(x);
    aL = accL(t);
    f = [vL-v; -par.p(v); aL];
    g = repmat([0; 1; 0], 1,size(x,2));
end

% Right-hand side
function dxdt = rhs(t,x,par)
    [f,g] = sys(t,x,par);
    u = k(t,x,par);
    dxdt = f+g*u;
end

%% Functions for control
% Desired controller
function ud = kd(t,x,par)
    [D,v,vL] = states(x);
    aL = accL(t);
    ud = par.A*(min(par.kappa*(D-par.Dst),par.vmax)-v) + ...
         par.B*(min(vL,par.vmax)-v) + ...
         par.C*aL;
end

% CBF evaluation
function [h,gradh] = CBF(x,par)
    [D,v] = states(x);
    h = (D-par.Dsf)/par.T - v;
    gradh = repmat([1/par.T; -1; 0], 1,size(x,2));
end

% Safe controller
function [us,h] = ks(t,x,par)
    [h,gradh] = CBF(x,par);
    [f,g] = sys(t,x,par);
    Lfh = dot(gradh,f,1);
    Lgh = dot(gradh,g,1);
    us = -(Lfh+par.alpha(h))./Lgh;
end

% Controller
function [u,ud,h] = k(t,x,par)
    % desired controller
    ud = kd(t,x,par);
    % safety filter
    switch par.filter
        case 'off'
            h = CBF(x,par);
            u = ud;
        case 'on'
            [us,h] = ks(t,x,par);
            u = min(ud,us);
    end
end