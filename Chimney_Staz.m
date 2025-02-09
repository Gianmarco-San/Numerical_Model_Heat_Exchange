%% CHIMNEY PROGRAM

clear all
close all
clc


%% DATI
%tutti le grandezze sono espresse in unità di misura internazionali (SI)

%GEOMETRICI
L1 = 6;      % [m] lunghezza tratto interno edificio
L2 = 3;      % [m] lunghezza tratto esterno edificio
L_t = L1+L2; % [m] lunghezza totale
r_e = 0.15;  % [m] raggio
s = 0.05;    % [m] spessore lamiera canna fumaria 
r_i = r_e-s; % [m] raggio interno

%TEMPERATURE
T_out = -8+273;   % [K] Temperatura dell'ambiente esterno
T_amb = 18+273;   % [K] Temperatura ambiente interno all'abitazione
T_fumi = 316+273; % [K] Temperatura dei fumi in ingressi nella canna fumaria

ha_in = 13.1; % [W/(m^2 K)] Coeff convettivo aria interno alla casa
ha_out = 10;  % [W/(m^2 K)] Coeff convettivo aria esterno all'abitazione

%FUMI
cp_fumi = 1029;         % [J/(kg K)] Calore specififico a pressione costante fumi
rho_fumi = 0.69;        % [kg/m^3] Densità fumi
lambda_fumi = 0.0395;   % [W/(m K)] Conducibilità termica fumi

%ACCIAIO
c_acciaio = 900;        % [J/(kg K)] Colore specico acciaio
rho_acciaio = 7900;     % [kg/m^3] Densità acciaio
lamdba_acciaio = 49.8;  % [W/(m K)] Conducibilità termica acciaio

%VELOCITA'
w_media = [0  1.0262];  % [m/s] velocità fluido con profilo costante
w_parabolica = @(x,y) [0.*y , (-3*w_media(2)/r_i^2*x.^2+3*w_media(2))]; %[m/s] velocità parabolica del fluido
    %nel condotto, garantendo la stessa portata del profilo lineare
    
% %STUDIO SENSITIVITA'
% w_media = [0  0.2];  % [m/s] velocità fluido con profilo costante
% w_parabolica = @(x,y) [0.*y , (-3*w_media(2)/r_i^2*x.^2+3*w_media(2))]; %[m/s] velocità parabolica del fluido
%     %nel condotto, garantendo la stessa portata del profilo lineare

%% DEFINIAMO LA GEOMETRIA

A_fumi = regions.rectN([0 0] , [r_e-s L_t]);    %Regione passaggio fumi
A_in = regions.rectN([r_i 0], [r_e L1]);        %Regione tubo interno alla casa
A_out = regions.rectN([r_i L1], [r_e L_t]);     %Regione tubo esterna alla casa

dom = [A_in A_out A_fumi];

dom(1).Name = 'Lamiera in';
dom(2).Name = 'Lamiera out';
dom(3).Name = 'Fumi';

%% CONDIZIONI AL CONTORNO

% axis equal
% dom.draw('e')

dom(3).Borders(1).Bc(2) = boundaries.neumann(0);
dom(3).Borders(1).Bc(3) = boundaries.neumann(0);
dom(1).Borders(1).Bc(1) = boundaries.neumann(0);
dom(2).Borders(1).Bc(3) = boundaries.neumann(0);

dom(3).Borders(1).Bc(1) = boundaries.dirichlet(T_fumi);

dom(1).Borders(1).Bc(4) = boundaries.robin(ha_in, ha_in*T_amb);
dom(2).Borders(1).Bc(4) = boundaries.robin(ha_out, ha_out*T_out);

dom(1).Borders(1).Bc(3) = boundaries.none;
dom(2).Borders(1).Bc(1) = boundaries.none;
dom(3).Borders(1).Bc(4) = boundaries.none;
dom(1).Borders(1).Bc(2) = boundaries.none;
dom(2).Borders(1).Bc(2) = boundaries.none;

% figure
% axis equal
% dom.draw('bc')

%% Propietà delle regioni

%FUMI

dom(3).rho = 1;     
dom(3).mu = lambda_fumi / (rho_fumi*cp_fumi); 
dom(3).beta = w_parabolica;

%LAMIERA:
    %interna

dom(1).rho = 1;
dom(1).mu = lamdba_acciaio/(rho_acciaio*c_acciaio);

    %esterna
dom(2).rho = dom(1).rho;
dom(2).mu = dom(1).mu;

%% MESH

Me = mesh2D(dom, 10^(-5.5));

% figure
% axis equal
% Me.draw

%% COSTRUIRE MATRICE DI RIGIDEZZA

[D,b] = BuildStiff_Cylindrical_Staz(Me);

% %SOLUZIONE CON GRADIENTE CONIUGATO
% opts.type = 'ilutp';
% opts.droptol = 1e-3;
% [L,U] = ilu(D,opts);
% T = bicg(D,b,1e-10,1000,L,U);

% SOLUZIONE CON METODO DI ELIMINAZIONE DI GAUSS
T = D\b;

TT=Me.copyToAllNodes(T);

%PLOT DISTRIBUZIONE TEMPERATURA
figure;
Me.draw(TT,'hidemesh'); 
%Render grafico
shading interp; 
colormap jet
xlabel('r dir [m]');
ylabel('z dir [m]');
zlabel('Temperature [K]');

%PLOT PROFILO TEMPERATURA (sez. ingresso, intermedia, uscita)
figure
hold on

    y_in = 0;
    xx = linspace(0,r_e,1000)'; %x equispaziati, mille punti
    yy = zeros(size(xx))+y_in; 
    Tinterp = Me.interpolate(TT,[xx,yy]); 
    plot(xx,Tinterp);

 y_intermedio = L1; 
    xx = linspace(0,r_e,1000)'; %x equispaziati, mille punti
    yy = zeros(size(xx))+y_intermedio; 
    Tinterp = Me.interpolate(TT,[xx,yy]); 
    plot(xx,Tinterp);

 y_out = L_t; 
    xx = linspace(0,r_e,1000)'; %x equispaziati, mille punti
    yy = zeros(size(xx))+y_out; 
    Tinterp = Me.interpolate(TT,[xx,yy]); 
    plot(xx,Tinterp);


xlabel('x dir [m]');
ylabel('Temperature [K]');
grid;
legend('z=0 m','z=6 m','y=9 m');

%FLUSSO
[Tx,Ty]=Me.gradient(TT); 
    figure; 
quiver(Me.Nodes.X,Me.Nodes.Y,Tx,Ty) 
title('Flux');
xlabel('x dir [m]');
ylabel('y dir [m]');
   
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

fprintf('\n Le Temperature medie calcolate in tre sezioni differenti (ingresso, sezione intemedia, uscita) sono, nel caso stazionario: \n %f K in z=%f m \n %f K in z=%f m \n %f K in z=%f m \n',T_media_in,z_ingresso,T_media_inter,z_intemedia,T_media_out,z_out)




