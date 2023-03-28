%% Safe sets of connected cruise control
clear; clc;

%% Problem setup
Dsf=1;      % safe distance [m]
T=1/0.6;	% time headway / time to conflict [s]

% Plot settings
ID=[0,35];  LD='distance, $D$ (m)';
Iv=[0,15];  Lv='speed, $v$ (m/s)';
IvL=[0,15];	LvL='speed, $v_{\rm L}$ (m/s)';
vLs=[5,10];
blue=[0,0,1];
darkgreen=[0,170,0]/256;
darkred=[230,0,0]/256;
lightblue=[212,218,249]/256;
lightgreen=[185,217,190]/256;
lightred=[248,150,150]/256;

%% Safe set boundaries
DD=@(v)Dsf+0*v;
TH=@(v)Dsf+T*v;
TTC=@(v,vL)Dsf+T*(v-vL);

%% Safe set illustrations
% 2D illustration
figure(3); clf; hold on; box on;
plot(DD(Iv),Iv,'Color',blue,'LineWidth',2,'HandleVisibility','off');
plot(TH(Iv),Iv,'Color',darkgreen,'LineWidth',2,'HandleVisibility','off');
plot(TTC(Iv.',vLs),Iv,'Color',darkred,'LineWidth',2,'HandleVisibility','off');
PlotFinalize({LD,Lv},[ID,Iv]);

% 3D illustration
figure(4); clf; hold on; box on;
[V,VL]=meshgrid(Iv,IvL);
surf(DD(V),V,VL,'FaceColor',lightblue,'FaceAlpha',0.75,...
                'EdgeColor',blue,'LineWidth',1,'HandleVisibility','off');
surf(TH(V),V,VL,'FaceColor',lightgreen,'FaceAlpha',0.75,...
                'EdgeColor',darkgreen,'LineWidth',1,'HandleVisibility','off');
surf(TTC(V,VL),V,VL,'FaceColor',lightred,'FaceAlpha',0.75,...
                    'EdgeColor',darkred,'LineWidth',1,'HandleVisibility','off');
PlotFinalize({LD,Lv,LvL},[ID,Iv,IvL]);
pbaspect([1,0.5,0.5]);
view([-15,45]);