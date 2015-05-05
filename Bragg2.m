%author JK, TW, PR

%functions:
deg2rad = @(deg) pi * deg/180;
%source:

sigma = 1e-6;                       %przewodnosc dla wody
mu = 1.2566370614359171*1e-6;       %mu
eps = 8.8541878176203892*1e-12;     %epsilon
cc = sqrt (1/ mu / eps);            %predkosc swiatla
N = 200;                            %rozmiar siatki
Nout = 20;
Ez = zeros(N+2*Nout,N+2*Nout);                    %notacja nie macierzowa
Hx= zeros(N+2*Nout,N+2*Nout+1);                   %to znaczy uzywamy (wymiar_poziomy_x, wymiar pionowy_y)
Hy= zeros(N+2*Nout+1,N+2*Nout);
                                    
L = 1;                              %rozmiar naszego pola
dx = L/N;
dy = L/N;                               
dt = dx /cc/sqrt(2);                   %krok czasowy
x = linspace(0, L, N);
y = linspace(0, L, N);
T_end = 1000 ;                      %koniec symulacji
t = 0;
HxN = zeros(N+Nout,1);                   %vector poziomy
HxS = zeros(N+Nout,1);
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
%Hye = [HyW; Hy; HyE];
%Hxe = [HxS, Hx, HxN]; 
Hye = Hy
Hxe = Hx

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
    %Hye = [HyW; Hye(2:end-1,:); HyE];
    %Hxe = [HxS, Hxe(:,2:end-1), HxN]; 
    t = t + dt ;  %update czasu
    
    %TE
    Eze = (1 - sigma*dt/eps)/(1 + sigma*dt/eps) * Eze + ...
        dt/(eps*dx)/(1+sigma*dt/(2*eps))*(Hye(2:end,:)- Hye(1:end-1,:))+ ...
        dt/(eps*dy)/(1+sigma*dt/(2*eps))*(Hxe(:,1:end-1) - Hxe(:,2:end));
    Hxe(:,2:end-1) = (1 - sigma*dt/2/mu)/(1 + sigma*dt/2/mu) * Hxe(:,2:end-1)+...
        (dt/(mu*dy))/(1+sigma*dt/(2*mu))*(Eze(:,1:end-1) - Eze(:,2:end));
    Hye(2:end-1,:) = (1 - sigma*dt/2/mu)/(1 + sigma*dt/2/mu) * Hye(2:end-1,:)+...
        (dt/(mu*dx))/(1+sigma*dt/(2*mu))*(Eze(2:end,:) - Eze(1:end-1,:));
    
    %Magia z oddzialywaniem m0_inc = 2fali padajacej na nasza probke, cos tu nie gra
    %bo tracimy symetrie...
    
    n = [1:N];
    dNe = k_inc * rc(n, N , m0); %wektory d dla kazdego boku
    dSe = k_inc * rc(n, 1 , m0); %poludnie
    dWe = k_inc * rc(1, n, m0); %
    dEe = k_inc * rc(N, n, m0);

    dNh = k_inc * rc(n, N+0.5 , m0)+0.5; %wektory d dla kazdego boku
    dSh = k_inc * rc(n, 1 , m0)+0.5; %poludnie
    dWh = k_inc * rc(1, n, m0)+.5; %
    dEh = k_inc * rc(N+0.5, n, m0)+.5;

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
    Eze(Nout+1:end-Nout,end-Nout) = Eze(Nout+1:end-Nout,end-Nout) + dt/(dx*eps).* HxN;
    Hx(Nout+1:end-Nout,end-Nout-1) =  Hx(Nout+1:end-Nout,end-Nout-1) + dt/(dx*mu).* EzN;
    %southFace
    Eze(Nout+1:end-Nout,Nout+1) = Eze(Nout+1:end-Nout,Nout+1) - dt/(dx*eps).* HxS;
    Hx(Nout+1:end-Nout,Nout+1) = Hx(Nout+1:end-Nout,Nout+1) - dt/(dx*mu).* EzS;
    %westFace
    Eze(Nout+1,Nout+1:end-Nout) = Eze(Nout+1,Nout+1:end-Nout) - dt/(dx*eps).* HyW;
    Hy(Nout+1,Nout+1:end-Nout) = Hy(Nout+1,Nout+1:end-Nout) - dt/(dx*mu).* EzW;
    %eastFace:
    Eze(end-Nout,Nout+1:end-Nout) = Eze(end-Nout,Nout+1:end-Nout) + dt/(dx*eps).* HyE;
    Hy(end-Nout-1,Nout+1:end-Nout) = Hy(end-Nout-1,Nout+1:end-Nout) + dt/(dx*mu).* EzE;
    
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

