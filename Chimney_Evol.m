%% CHIMNEY PROGRAM

clear all
close all
clc

%% DATI

L2 = 3;         % [m] Lunghezza parte esterna canna fumaria
L1 = 6;         % [m] Lunghezza parte interna canna fumaria
L_t = L1+L2;    % [m ]Lunghezza totale canna fumaria

d_i = 0.2;      % [m] Diametro canna fumaria
r_i = d_i/2;    % [m] Raggio canna fumaria
s = 0.05;       % [m] Spessore canna fumaria
d_e = d_i+2*s;  % [m] Diametro esterno canna fumaria
r_e = d_e/2;    % [m] Raggio esterno canna fumaria

rho_fumi = 0.69;    % [kg/m^3] Densità fumi
rho_acciaio = 7900; % [kg/m^3] Densità condotto (materiale AISI 1010)

c_condotto = 900;   % [J/(kg K)] Calore specifico condotto
cp_fumi=1029;       % [J/(kg K)] Calore specifico a pressione costante fumi

lambda_acciaio = 49.8;  % [W/(m K)] Conducibilità condotto (materiale AISI 1010) 
lambda_fumi=0.0395;     % [W/(m K)] Conducibilità fumi

w_media = [0 1.0262]; % [m/s] Velocità media fumi nello scambiatore di calore
w_parabolica = @(x,y) [0.*y , (-3*w_media(2)/r_i^2*x.^2+3*w_media(2))];

T_out = -8+273;     % [K] Temperatura esterna -8°C
T_amb = 18+273;     % [K] Temperatura interna 18°C
T_fumi = 316+273;   % [K] Temperatura ingresso fumi 316.15°C

ha_in = 13.1;       % [W/(m^2 K)] Coefficiente scambio termico parte interna diviso per rho_condotto e c_condotto
ha_out = 10;        % [W/(m^2 K)] Coefficiente scambio termico parte esterna diviso per rho_condotto e c_condotto

mu_acciaio = lambda_acciaio/(rho_acciaio*c_condotto);
mu_fumi = lambda_fumi/(rho_fumi*cp_fumi);

% %STUDIO SENSITIVITA'
% w_media = [0  0.2];  % [m/s] velocità fluido con profilo costante
% w_parabolica = @(x,y) [0.*y , (-3*w_media(2)/r_i^2*x.^2+3*w_media(2))]; %[m/s] velocità parabolica del fluido
%     %nel condotto, garantendo la stessa portata del profilo lineare

%% DEFINIZIONE GEOMETRIA

dom_fumi = regions.rectN([0 0],[r_i L_t], 'beta', w_parabolica,'rho',1, 'mu',mu_fumi,'name','fumi');
dom_acciaio1 = regions.rectN ([r_i 0],[r_e,L1],'rho',1, 'mu',mu_acciaio,'name','condotto interno');
dom_acciaio2 = regions.rectN ([r_i L1],[r_e, L_t],'rho',1, 'mu',mu_acciaio,'name','condotto esterno');

dom = [dom_fumi,dom_acciaio1,dom_acciaio2];

% figure
% dom.draw
% axis equal
%% CONDIZIONI AL CONTORNO

dom(1).Borders(1).Bc(4) = boundaries.none;
dom(2).Borders(1).Bc(2) = boundaries.none;
dom(2).Borders(1).Bc(3) = boundaries.none;
dom(3).Borders(1).Bc(1) = boundaries.none;
dom(3).Borders(1).Bc(2) = boundaries.none;

dom(1).Borders(1).Bc(1) = boundaries.dirichlet(1);

dom(1).Borders(1).Bc(2) = boundaries.neumann(0);
dom(1).Borders(1).Bc(3) = boundaries.neumann(0);
dom(3).Borders(1).Bc(3) = boundaries.neumann(0);
dom(2).Borders(1).Bc(1) = boundaries.neumann(0);

dom(2).Borders(1).Bc(4) = boundaries.robin(ha_in, ha_in*T_amb);
dom(3).Borders(1).Bc(4) = boundaries.robin(ha_out, ha_out*T_out);

%% MESH

Me = mesh2D(dom,1e-5);

DirichletNodes_Start = Me.find(@(x,y)((y==0)&(x<=r_i)),'d'); %Trova i nodi che sono nel bordo di Dirichlet
Dof = Me.Nodes.Dof>0;

[D, bconst, bvar] = BuildStiffEvol_cil (Me);

%% SOLUZIONE STAZIONARIA

T_Staz0 = D\(bconst+bvar*T_amb); %Soluzione stazionaria inziale

%Plot della soluzione iniziale
TT=zeros(size(Me.Nodes.X));
TT(Dof>0)=T_Staz0;              %Nei gdl metto come temperatura iniziale la T_stazionaria 
TT(DirichletNodes_Start)=T_amb; %Pongo inizialmente nel bordo in ingresso la T_amb
figure
Me.draw(TT,'hidemesh')
shading interp;
colormap jet
ylabel(colorbar(),'Temperatura [K]');
title('SOLUZIONE STAZIONARIA')

%% EVOLUTZIONE NEL TEMPO - MATRICE DI MASSA

t_transitorio = 3600;   % Tempo che il sistema impiega per arrivare a regime
t_end = 4000;           %Tempo di fine valutazione

%Legge di variazione LINEARE che porta il sistema nella stato stazionario
%in un tempo t_transitorio:
% fDirichlet = @(t)T_amb + (t<=t_transitorio)*(T_fumi - T_amb)/t_transitorio*t +... 
%     (t>t_transitorio)*(T_fumi - T_amb);

%Legge di variazione PARABOLICA che porta il sistema nella stato stazionario
    %in un tempo t_transitorio:
fDirichlet = @(t) T_amb + (t<=t_transitorio).*sqrt((T_fumi-T_amb)^2/t_transitorio.*t) + ...
    (t>t_transitorio)*(T_fumi-T_amb);

%% MATRICE DI MASSA

[M, mvar] = BuildMassEvol_cil (Me , r_i);

%% EULERO IMPLICITO

figure

disp('Soluzione svolta con Eulero implicito');


T  = T_Staz0;
dt = 100;       %intervallo di tempo dello studio

A = (M + D*dt); 

tic
T_iniziale = fDirichlet(0);

for k=1:t_end/dt
    
    t=k*dt;
    T_new = fDirichlet(t);
    T = A\(M*T+dt*(bconst+bvar*T_new)-mvar*(T_new-T_iniziale));
    TT(Dof)=T;
    TT(DirichletNodes_Start) = T_new;
    T_iniziale = T_new;
    
    hold off; 
    
        Me.draw(TT,'hidemesh');
        shading interp;
        colormap jet
        zlim([T_out-10 T_fumi+10]);
        caxis([T_out-10 T_fumi+10]);
        title(['t = ' num2str(t) 's']);    
        xlabel('r dir [m]');
        ylabel('z dir [m]');
        ylabel(colorbar(),'Temperature [K]');
        
        drawnow();

end
toc

%% CALCOLO MEDIA DELLA TEMPERATURA

%SEZIONE INGRESSO
z_ingresso = 0;

xx = linspace(0, r_e, 1000)';
yy = zeros(size(xx))+z_ingresso;

Tinterp = Me.interpolate(TT,[xx,yy]);

T_media_in = sum(Tinterp)/size(Tinterp,1); %Temperatura media nella sezione in ingresso

%SEZIONE INTERMEDIA
z_intemedia = L1;

xx = linspace(0, r_e, 1000)';
yy = zeros(size(xx))+z_intemedia;

Tinterp = Me.interpolate(TT,[xx,yy]);

T_media_inter = sum(Tinterp)/size(Tinterp,1); %Temperatura media nella sezione intermedia

%SEZIONE USCITA
z_out = L_t;

xx = linspace(0, r_e, 1000)';
yy = zeros(size(xx))+z_out;

Tinterp = Me.interpolate(TT,[xx,yy]);

T_media_out = sum(Tinterp)/size(Tinterp,1); %Temperatura media nella sezione in ingresso

fprintf('\n Le Temperature medie calcolate in tre sezioni differenti (ingresso, sezione intemedia, uscita) sono, dopo il transitorio: \n %f K in z=%f m \n %f K in z=%f m \n %f K in z=%f m \n',T_media_in,z_ingresso,T_media_inter,z_intemedia,T_media_out,z_out)
