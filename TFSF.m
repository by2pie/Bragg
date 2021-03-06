%author JK, TW, PR
clear all; clc;
%functions:
deg2rad = @(deg) pi * deg/180;
%source:

%%%%%%%%%%%%%%%%
% READ ME/TODO %
%%%%%%%%%%%%%%%%
% (0) Pole Ez powinno by� rozbite na sk�adowe Ezx i Ezy - potrzebne do
% PMLi - do zrobienia przy implementacji warunk�w brzegowych
% --------------------
% (1) Na koniec powinno to wygl�da� w ten spos�b (brak skali :P):
%  __________________
% |   ___PML______   |
% |  |  __SF___   |  |
% |  | | TF-SF |  |  |
% |  | |_______|  |  |
% |  |____________|  |
% |__________________|
% warstwa PML ma grubo�� d, obszar obliczeniowy N, gdzie N = n + 2*n_scat
% ostatecznie pola powinny by� zdefiniowane na macierzach o wymiarach
% (N+2*d)*(N+2*d) [pole Ezy i Ezx], (N+2*d)*(N+2*d+1) [Hx], (N+2*d+1)*(N+2*d) [Hy]
% --------------------
% (2) Wszystkie wsp�czynniki typu: (1 - sigma*dt/eps_0)/(1 + sigma*dt/eps_0)
% powinny by� zdefiniowane jako macierze o wymiarach N*N - wsp�czynniki 
% aktualizacyjne. Jest o istotne,je�eli chcemy by cz�� obszaru 
% obliczeniowego odwzorowywa�a jaki� inny materia�.
% --------------------
% (3) Rysowanie wykres�w przerobione na cienowane powierzchnie - wydaje si�
% dzia�a� nieco lepiej, w ka�dym razie jest kolorowo.

%(0) Coś przesuwa się symetrycznie, być może brakuje PML a moze strzelone sa indeksy  


sigma = 0;%1e-5;                    % przewodnosc dla wody
mu_0 = 1.2566370614359171*1e-6;       % przenikalno�� magnetyczna pr�ni
eps_0 = 8.8541878176203892*1e-12;     % przenikalno�� elektryczna pr�ni
mu_r = 1.0;                         % wzgl�dna przenikalno�� magnetyczna
eps_r = 1.0;                        % wzgl�dna przenikalno�� elektryczna
cc = sqrt (1/mu_0/eps_0/mu_r/eps_r);    % pr�dko�� �wiat�a
     



L = 1;                              % rozmiar pola [metr? milimetr?]
N_TF_SF = 100;                            % rozmiar TF_SF
N_SF = 20;                         % rozmiar pola SF
N = N_TF_SF + 2*N_SF
interior = false(N,N)
interior(N_SF+1:N-N_SF,N_SF+1:N-N_SF) = true
dx = L/N;                           % krok przestrzenny
dy = L/N;                               
dt = dx /cc/sqrt(2);                % krok czasowy
T_end = 1000 ;                      % koniec symulacji
t = 0;


eps(interior) = eps_0 * eps_r
mu(interior) = mu_0 * mu_r


Ez = zeros(N,N);                    %notacja nie macierzowa
Hx= zeros(N,N+1);                   %to znaczy uzywamy (wymiar_poziomy_x, wymiar pionowy_y)
Hy= zeros(N+1,N);


E0 = 1;                             % moc �r�d�a     
m0 = [1,1];                         % po�o�enie �r�d�a - naro�nik SW
m0_inc = 4;                         % po�o�enia �r�d�a w bazie padaj�cej fali
inc_angle = 45 ;                    % k�t pod kt�rym pada fala [stopnie]
angle = deg2rad(inc_angle);         % przeliczenie k�ta na radiany

E_inc = zeros(ceil(N*sqrt(2)) + m0_inc,1);
H_inc = zeros(ceil(N*sqrt(2)) + m0_inc,1);
k_inc = [cos(angle), sin(angle)];


V = 1;

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
Hye = Hy
Hxe = Hx

%% PML
eta = sqrt(mu_0*mu_r/eps_r/eps_0);  % impedancja
r = 1e-5;                           % maksymalne dopuszczalne "odbicie"
d = 20;                             % szeroko�� warstwy PML

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
plotE = surf(X, Y, Eze(:,:)*0);
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

%% P�TLA
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
    
    t = t + dt ;  %update czasu

     %TE
    Eze = (1 - sigma*dt/eps_0)/(1 + sigma*dt/eps_0) * Eze + ...
        dt/(eps_0*dx)/(1+sigma*dt/(2*eps_0))*(Hye(2:end,:)- Hye(1:end-1,:))+ ...
        dt/(eps_0*dy)/(1+sigma*dt/(2*eps_0))*(Hxe(:,1:end-1) - Hxe(:,2:end));
    
    Hxe(:,2:end-1) = (1 - sigma*dt/2/mu_0)/(1 + sigma*dt/2/mu_0) * Hxe(:,2:end-1)+...
        (dt/(mu_0*dy))/(1+sigma*dt/(2*mu_0))*(Eze(:,1:end-1) - Eze(:,2:end));
    
    Hye(2:end-1,:) = (1 - sigma*dt/2/mu_0)/(1 + sigma*dt/2/mu_0) * Hye(2:end-1,:)+...
        (dt/(mu_0*dx))/(1+sigma*dt/(2*mu_0))*(Eze(2:end,:) - Eze(1:end-1,:));
    
    n = [1:N_TF_SF];
    
    dNe = k_inc * rc(n, N_TF_SF , m0); %wektory d dla kazdego boku
    dSe = k_inc * rc(n, 1 , m0); %poludnie
    dWe = k_inc * rc(1, n, m0); %
    dEe = k_inc * rc(N_TF_SF, n, m0);

    dNh = k_inc * rc(n, N_TF_SF+0.5 , m0)+0.5; %wektory d dla kazdego boku
    dSh = k_inc * rc(n, 1 , m0)+0.5; %poludnie
    dWh = k_inc * rc(1, n, m0)+.5; %
    dEh = k_inc * rc(N_TF_SF+0.5, n, m0)+.5;

    indNe = floor(dNe);
    indSe = floor(dSe);
    indWe = floor(dWe);
    indEe = floor(dEe);
    
    indNh = floor(dNh);
    indSh = floor(dSh);
    indWh = floor(dWh);
    indEh = floor(dEh);
    
    dNe = dNe - indNe;
    dSe = dSe - indSe;
    dWe = dWe - indWe;
    dEe = dEe - indEe;

    dNh = dNh - indNh;
    dSh = dSh - indSh;
    dWh = dWh - indWh;
    dEh = dEh - indEh;

    HxN = ((1-dNh).*H_inc(m0_inc+indNh)' + dNh.*H_inc(m0_inc+indNh+1)') * sin(angle);
    HxS = ((1-dSh).*H_inc(m0_inc+indSh)' + dSh.*H_inc(m0_inc+indSh+1)') * sin(angle);
    HyW = -((1-dWh).*H_inc(m0_inc+indWh)' + dWh.*H_inc(m0_inc+indWh+1)') * cos(angle);
    HyE = -((1-dEh).*H_inc(m0_inc+indEh)' + dEh.*H_inc(m0_inc+indEh+1)') * cos(angle);
    HxN = HxN';
    HxS = HxS';
    
    EzN = ((1-dNe).*H_inc(m0_inc+indNe)' + dNe.*H_inc(m0_inc+indNe+1)');
    EzS = ((1-dSe).*H_inc(m0_inc+indSe)' + dSe.*H_inc(m0_inc+indSe+1)');
    EzW = ((1-dWe).*H_inc(m0_inc+indWe)' + dWe.*H_inc(m0_inc+indWe+1)');
    EzE = ((1-dEe).*H_inc(m0_inc+indEe)' + dEe.*H_inc(m0_inc+indEe+1)');
    
    EzN = EzN';
    EzS = EzS';
    %NorthFace
    Eze(N_SF+1:end-N_SF,end-N_SF) = Eze(N_SF+1:end-N_SF,end-N_SF) + dt/(dx*eps_0).* HxN;
    Hx(N_SF+1:end-N_SF,end-N_SF-1) =  Hx(N_SF+1:end-N_SF,end-N_SF-1) + dt/(dx*mu_0).* EzN;
    %southFace
    Eze(N_SF+1:end-N_SF,N_SF+1) = Eze(N_SF+1:end-N_SF,N_SF+1) - dt/(dx*eps_0).* HxS;
    Hx(N_SF+1:end-N_SF,N_SF+1) = Hx(N_SF+1:end-N_SF,N_SF+1) - dt/(dx*mu_0).* EzS;
    %westFace
    Eze(N_SF+1,N_SF+1:end-N_SF) = Eze(N_SF+1,N_SF+1:end-N_SF) - dt/(dx*eps_0).* HyW;
    Hy(N_SF+1,N_SF+1:end-N_SF) = Hy(N_SF+1,N_SF+1:end-N_SF) - dt/(dx*mu_0).* EzW;
    %eastFace:
    Eze(end-N_SF,N_SF+1:end-N_SF) = Eze(end-N_SF,N_SF+1:end-N_SF) + dt/(dx*eps_0).* HyE;
    Hy(end-N_SF-1,N_SF+1:end-N_SF) = Hy(end-N_SF-1,N_SF+1:end-N_SF) + dt/(dx*mu_0).* EzE;

    set(plotE,'Zdata',Eze(:,:)^2);
    set(plotHx,'Zdata',Hxe(:,2:end));
    set(plotHy,'Zdata',Hye(2:end,:));
    set(time,'String',sprintf('t = %0.2e [s]',t));
    refresh;
    pause(0.01);
end


