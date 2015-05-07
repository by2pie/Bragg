%author JK, TW, PR

%functions:
deg2rad = @(deg) pi * deg/180;
%source:

sigma = 1e-5;                       %przewodnosc dla wody
mu = 1.2566370614359171*1e-6;       %mu
eps = 8.8541878176203892*1e-12;     %epsilon
cc = sqrt (1/ mu / eps);            %predkosc swiatla
N = 200;                            %rozmiar siatki
N_scatter = 20;
Ez = zeros(N+N_scatter,N+N_scatter);                    %notacja nie macierzowa
Hx= zeros(N+N_scatter,N+N_scatter+1);                   %to znaczy uzywamy (wymiar_poziomy_x, wymiar pionowy_y)
Hy= zeros(N+N_scatter+1,N+N_scatter);
                                    
L = 1;                              %rozmiar naszego pola
dx = L/N;
dy = L/N;                               
dt = dx /cc/sqrt(2);                   %krok czasowy
x = linspace(0, L, N);
y = linspace(0, L, N);
T_end = 1000 ;                      %koniec symulacji
t = 0;
HxN = zeros(N,1);                   %vector poziomy
HxS = zeros(N,1);
HyW = zeros(1,N);                   %wektor pionowy
HyE = zeros(1,N);


E0 = 1; %moc zrodla     
m0 = [1,1]; %polozenie zrodla naroznik SW
m0_inc = 4 %polozenia zrodla w bazie padajacej fali
inc_angle = 45 ;    %fala nadciaga pod kate 45 stopni
angle = deg2rad(inc_angle);  

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
figure('position', [0,0,600,1000]);
etemp2 = 0;
etemp1 = 0;
htemp1 = 0;
htemp2 = 0;
%zrodlo:
g = @(t) exp(-.5*(( t/dt - 40 ) / 12 ) ^2);

%kopie pol
Eze = Ez;
Hye = Hy;
Hxe = Hx; 

%glowna petla
while (t < T_end)
    %source:
    drawnow 
    %incident wave - policzona jednowymiarowa fala padajaca
    E_inc(2:end) = E_inc(2:end) + dt/(eps*V*dx) * (H_inc(2:end) - H_inc(1:end-1));
    E_inc(2) = E0 * g(t);   %zrodlo jest w drugim wyrazie i jest to m0-2 wiec m0 = 2 + 2 =4
    E_inc(1) = etemp2;
    etemp2 = etemp1;
    etemp1 = E_inc(2);    
    H_inc(1:end-1) = H_inc(1:end-1) + dt/(mu*V*dx) * (E_inc(2:end) - E_inc(1:end-1));
    H_inc(end) = htemp2;
    htemp2= htemp1;
    htemp1 = H_inc(end-1);
    
    %laczenie warunkow brzegowych (niepotrzebne)

    t = t + dt ;  %update czasu
    size(Eze)
    size(Hxe)
    size(Hye)
    
    
    
    %TE
    'Eze'
    size(Eze(2:end,2:end))
    'Hx'
    size(Hxe(:,1:end-1) - Hxe(:,2:end))
    'Hy'
    size(Hye(2:end,:)- Hye(1:end-1,:))
    pause
    Eze = (1 - sigma*dt/eps)/(1 + sigma*dt/eps) * Eze + ...
        dt/(eps*dx)/(1+sigma*dt/(2*eps))*(Hye(2:end,:)- Hye(1:end-1,:))+ ...
        dt/(eps*dy)/(1+sigma*dt/(2*eps))*(Hxe(:,1:end-1) - Hxe(:,2:end));
    Hxe(:, 2:end) = (1 - sigma*dt/2/mu)/(1 + sigma*dt/2/mu) * Hxe +...
        (dt/(mu*dy))/(1+sigma*dt/(2*mu))*(Eze(:,1:end-1) - Eze(:,2:end));
    Hye(2:end, :) = (1 - sigma*dt/2/mu)/(1 + sigma*dt/2/mu) * Hye +...
        (dt/(mu*dx))/(1+sigma*dt/(2*mu))*(Eze(2:end,:) - Eze(1:end-1,:));
    
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
    Eze(:,end) = Eze(:,end) + dt/(dx*eps).* HxN;
    
    %southFace
    Eze(:,1) = Eze(:,1) + dt/(dx*eps).* HxS;
    
    %westFace;
    Eze(1,:) = Eze(1,:) + dt/(dx*eps).* HyW;
    
    %eastFace:
    Eze(end,:) = Eze(end,:) + dt/(dx*eps).* HyE;
    subplot(2,1,1)
    %plot natezenia pola E
    contourf(Eze)
    axis equal
    %contourf(Eze, 'LevelList', [-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5])
    subplot(2,1,2)
    %plot tej fali padajacej
    plot(E_inc)
    axis auto
end


