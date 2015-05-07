%author JK, TW, PR
clear all; clc;
%functions:
deg2rad = @(deg) pi * deg/180;
%source:

%%%%%%%%%%%%%%%%
% READ ME/TODO %
%%%%%%%%%%%%%%%%
% (0) Pole Ez powinno byæ rozbite na sk³adowe Ezx i Ezy - potrzebne do
% PMLi - do zrobienia przy implementacji warunków brzegowych
% --------------------
% (1) Na koniec powinno to wygl¹daæ w ten sposób (brak skali :P):
%  __________________
% |   ___PML______   |
% |  |  __SF___   |  |
% |  | | TF-SF |  |  |
% |  | |_______|  |  |
% |  |____________|  |
% |__________________|
% warstwa PML ma gruboœæ d, obszar obliczeniowy N, gdzie N = n + 2*n_scat
% ostatecznie pola powinny byæ zdefiniowane na macierzach o wymiarach
% (N+2*d)*(N+2*d) [pole Ezy i Ezx], (N+2*d)*(N+2*d+1) [Hx], (N+2*d+1)*(N+2*d) [Hy]
% --------------------
% (2) Wszystkie wspó³czynniki typu: (1 - sigma*dt/eps_0)/(1 + sigma*dt/eps_0)
% powinny byæ zdefiniowane jako macierze o wymiarach N*N - wspó³czynniki 
% aktualizacyjne. Jest o istotne,je¿eli chcemy by czêœæ obszaru 
% obliczeniowego odwzorowywa³a jakiœ inny materia³.
% --------------------
% (3) Rysowanie wykresów przerobione na cienowane powierzchnie - wydaje siê
% dzia³aæ nieco lepiej, w ka¿dym razie jest kolorowo.


sigma = 0;%1e-5;                    % przewodnosc dla wody
mu_0 = 1.2566370614359171*1e-6;       % przenikalnoœæ magnetyczna pró¿ni
eps_0 = 8.8541878176203892*1e-12;     % przenikalnoœæ elektryczna pró¿ni
mu_r = 1.0;                         % wzglêdna przenikalnoœæ magnetyczna
eps_r = 1.0;                        % wzglêdna przenikalnoœæ elektryczna
cc = sqrt (1/mu_0/eps_0/mu_r/eps_r);    % prêdkoœæ œwiat³a
                                    
L = 1;                              % rozmiar pola [metr? milimetr?]
N = 100;                            % rozmiar siatki - liczba wêz³ów
dx = L/N;                           % krok przestrzenny
dy = L/N;                               
dt = dx /cc/sqrt(2);                % krok czasowy
% x = linspace(0, L, N);
% y = linspace(0, L, N);
T_end = 1000 ;                      % koniec symulacji
t = 0;

% inicjalizacja pól
Ez = zeros(N,N);                    % notacja nie macierzowa
Hx= zeros(N,N-1);                   % to znaczy uzywamy (wymiar_poziomy_x, wymiar pionowy_y)
Hy= zeros(N-1,N);                   % [tw] N+1 zamiast N-1

% PMC - ? --> bêd¹ zast¹pione przez 
HxN = zeros(N,1);                   % wektor poziomy
HxS = zeros(N,1);
HyW = zeros(1,N);                   % wektor pionowy
HyE = zeros(1,N);

E0 = 1;                             % moc Ÿród³a     
m0 = [1,1];                         % po³o¿enie Ÿród³a - naro¿nik SW
m0_inc = 4;                         % po³o¿enia Ÿród³a w bazie padaj¹cej fali
inc_angle = 45 ;                    % k¹t pod którym pada fala [stopnie]
angle = deg2rad(inc_angle);         % przeliczenie k¹ta na radiany

E_inc = zeros(ceil(N*sqrt(2)) + m0_inc,1);
H_inc = zeros(ceil(N*sqrt(2)) + m0_inc,1);
k_inc = [cos(angle), sin(angle)];

% nei rozumiem co to jest za wspolczynnik predkosci fazowych z tego ksera
V = 1;
%V = ?

%indeksy:
%E(i,j) = E(i-0.5, j+0,5)
%Hx(i,j) = Hx(i-0.5, j+1)
%Hy(i,j) = Hy(i, j+0.5)
%H_inc(i) = H_inc(i-0.5)
%E_inc(i) = E_inc(i)

%inicjalizacja wykresu i zmiennych
% figure('position', [0,0,600,1000]);
etemp2 = 0;
etemp1 = 0;
htemp1 = 0;
htemp2 = 0;
%zrodlo:
g = @(t) exp(-.5*(( t/dt - 40 ) / 12 ) ^2);

%kopie pol
Eze = Ez;
Hye = [HyW; Hy; HyE];
Hxe = [HxS, Hx, HxN]; 

%% PML
eta = sqrt(mu_0*mu_r/eps_r/eps_0);  % impedancja
r = 1e-5;                           % maksymalne dopuszczalne "odbicie"
d = 20;                             % szerokoœæ warstwy PML

sigmaE_max = -1.5*log(r)/d/eta;
sigmaH_max = sigmaE_max*eta^2;

sigmaE_l = zeros(d, N);
sigmaH_l = zeros(d, N);

sigmaE_r = zeros(d, N);
sigmaH_r = zeros(d, N);

for i = 1:d
    sigmaE_l(i,:) = sigmaE_max*((d-i)/d)^2;
    sigmaE_r(i,:) = sigmaE_max*((i-1)/d)^2;
    
    sigmaH_l(i,:) = sigmaH_max*((d-i)/d)^2;
    sigmaH_r(i,:) = sigmaH_max*((i-1)/d)^2;
end

cEE_l = (2*eps_0*eps_r - dt.*sigmaE_l)./(2*eps_0*eps_r + dt.*sigmaE_l);
cHH_l = (2*mu_0*mu_r - dt.*sigmaH_l)./(2*mu_0*mu_r + dt.*sigmaH_l);
cEH_l = 2*dt./(2*eps_0*eps_r + dt.*sigmaE_l)/dx;
cHE_l = 2*dt./(2*mu_0*mu_r + dt.*sigmaH_l)/dx;

cEE_r = (2*eps_0*eps_r - dt.*sigmaE_r)./(2*eps_0*eps_r + dt.*sigmaE_r);
cHH_r = (2*mu_0*mu_r - dt.*sigmaH_r)./(2*mu_0*mu_r + dt.*sigmaH_r);
cEH_r = 2*dt./(2*eps_0*eps_r + dt.*sigmaE_r)/dx;
cHE_r = 2*dt./(2*mu_0*mu_r + dt.*sigmaH_r)/dx;

%% WYKRESY
figure;
set(gcf,'outerposition',[100   100   1200   600]);
figDim = get(gcf,'Position');

lblDim = [140 15];

time = uicontrol('style','text');
set(time,'BackgroundColor',[0.8 0.8 0.8],...
    'Position',[0.5*figDim(3)-0.5*lblDim(1) figDim(4)-lblDim(2) lblDim(1) lblDim(2)],...
    'String', 't = 0.0e0 [s]',...
    'FontSize',9);

x = linspace(1,N,N);
y = linspace(1,N,N);
[Y, X] = meshgrid(x, y);

subplot(1,3,1);
plotE = surf(X, Y, Eze(:,:));
title('E_z [V/m]');
axis equal; 
axis([1 N 1 N]);
xlabel('x'); ylabel('y');
shading interp;

subplot(1,3,2);
plotHx = surf(X, Y, Hxe(:,2:end));
title('H_x [A/m]');
axis equal; 
axis([1 N 1 N]);
xlabel('x'); ylabel('y');
shading interp;

subplot(1,3,3);
plotHy = surf(X, Y, Hye(2:end,:));
title('H_y [A/m]');
axis equal; 
axis([1 N 1 N]);
xlabel('x'); ylabel('y');
shading interp;

%% PÊTLA
while (t < T_end)
    %source:
    drawnow 
    %incident wave - policzona jednowymiarowa fala padajaca
    E_inc(2:end) = E_inc(2:end) + dt/(eps_0*V*dx) * (H_inc(2:end) - H_inc(1:end-1));
    E_inc(2) = E0 * g(t);   %zrodlo jest w drugim wyrazie i jest to m0-2 wiec m0 = 2 + 2 =4
    E_inc(1) = etemp2;
    etemp2 = etemp1;
    etemp1 = E_inc(2);    
    H_inc(1:end-1) = H_inc(1:end-1) + dt/(mu_0*V*dx) * (E_inc(2:end) - E_inc(1:end-1));
    H_inc(end) = htemp2;
    htemp2= htemp1;
    htemp1 = H_inc(end-1);
    
    %laczenie warunkow brzegowych (niepotrzebne)
    Hye = [HyW; Hye(2:end-1,:); HyE];
    Hxe = [HxS, Hxe(:,2:end-1), HxN];
    t = t + dt ;  %update czasu
    
    %TE
    Eze = (1 - sigma*dt/eps_0)/(1 + sigma*dt/eps_0) * Eze + ...
        dt/(eps_0*dx)/(1+sigma*dt/(2*eps_0))*(Hye(2:end,:)- Hye(1:end-1,:))+ ...
        dt/(eps_0*dy)/(1+sigma*dt/(2*eps_0))*(Hxe(:,1:end-1) - Hxe(:,2:end));
    Hxe(:,2:end-1) = (1 - sigma*dt/2/mu_0)/(1 + sigma*dt/2/mu_0) * Hxe(:,2:end-1)+...
        (dt/(mu_0*dy))/(1+sigma*dt/(2*mu_0))*(Eze(:,1:end-1) - Eze(:,2:end));
    Hye(2:end-1,:) = (1 - sigma*dt/2/mu_0)/(1 + sigma*dt/2/mu_0) * Hye(2:end-1,:)+...
        (dt/(mu_0*dx))/(1+sigma*dt/(2*mu_0))*(Eze(2:end,:) - Eze(1:end-1,:));
    
    %Magia z oddzialywaniem m0_inc = 2fali padajacej na nasza probke, cos tu nie gra
    %bo tracimy symetrie...
    
    n = [1:N];
    dN = k_inc * rc(n, N , m0); %wektory d dla kazdego boku
    dS = k_inc * rc(n, 1 , m0); %poludnie
    dW = k_inc * rc(1, n, m0); %
    dE = k_inc * rc(N, n, m0);

    indN = floor(dN);
    indS = floor(dS);
    indW = floor(dW);
    indE = floor(dE);
    
    
    dN = dN - indN;
    dS = dS - indS;
    dW = dW - indW;
    dE = dE - indE;

    HxN = ((1-dN).*H_inc(m0_inc+indN)' + dN.*H_inc(m0_inc+indN+1)') * sin(angle);
    HxS = ((1-dS).*H_inc(m0_inc+indS)' + dS.*H_inc(m0_inc+indS+1)') * sin(angle);
    HyW = ((1-dW).*H_inc(m0_inc+indW)' + dW.*H_inc(m0_inc+indW+1)') * cos(angle);
    HyE = ((1-dE).*H_inc(m0_inc+indE)' + dE.*H_inc(m0_inc+indE+1)') * cos(angle);
    HxN = HxN';
    HxS = HxS';
    
    %NorthFace
    Eze(:,end) = Eze(:,end) + dt/(dx*eps_0).* HxN;
    
    %southFace
    Eze(:,1) = Eze(:,1) + dt/(dx*eps_0).* HxS;
    
    %westFace;
    Eze(1,:) = Eze(1,:) + dt/(dx*eps_0).* HyW;
    
    %eastFace:
    Eze(end,:) = Eze(end,:) + dt/(dx*eps_0).* HyE;
%     subplot(2,1,1)
%     %plot natezenia pola E
%     contourf(Eze)
%     axis equal
%     %contourf(Eze, 'LevelList', [-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5])
%     subplot(2,1,2)
%     %plot tej fali padajacej
%     plot(E_inc)
%     axis auto
    set(plotE,'Zdata',Eze(:,:));
    set(plotHx,'Zdata',Hxe(:,2:end));
    set(plotHy,'Zdata',Hye(2:end,:));
    set(time,'String',sprintf('t = %0.2e [s]',t));
    refresh;
    pause(0.01);
end


