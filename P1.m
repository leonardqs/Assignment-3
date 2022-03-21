clear 
close all
 

%global C

% C.q_0 = 1.60217653e-19;             % electron charge
% C.hb = 1.054571596e-34;             % Dirac constant
% C.h = C.hb * 2 * pi;                % Planck constant
% C.m_0 = 9.10938215e-31;             % electron mass
% C.kb = 1.3806504e-23;               % Boltzmann constant
% C.eps_0 = 8.854187817e-12;          % vacuum permittivity
% C.mu_0 = 1.2566370614e-6;           % vacuum permeability
% C.c = 299792458;                    % speed of light
% C.g = 9.80665;                      % metres (32.1740 ft) per s²

q_0 = 1.60217653e-19;             % electron charge
hb = 1.054571596e-34;             % Dirac constant
h = hb * 2 * pi;                % Planck constant
m_0 = 9.10938215e-31;             % electron mass
kb = 1.3806504e-23;               % Boltzmann constant
eps_0 = 8.854187817e-12;          % vacuum permittivity
mu_0 = 1.2566370614e-6;           % vacuum permeability
c = 299792458;                    % speed of light
g = 9.80665;                      % metres (32.1740 ft) per s²


Length = 200e-9;                % Boundary length
Width = 100e-9;                % Boundary width
timeStep = 7e-15;                     % Time Step
totalParticle = 3e3;
displayParticle = 10;              % Display Particle in the box
SimulationTime = 200;                       % Simulation Time
T = 300;                            % Default Temperature
t_mn = 0.2e-12;                     % Mean time between collision 

VoL = 0.1; % Left side of the Area has Boundary Voltage = VoL
E = VoL/Length
F = E*q_0
Acceleration = F/m_0

mass_eff =  0.26*m_0;
initial_Vth = sqrt((2* kb *T)/mass_eff);

% Generate the initial electron position
x_pos = (Length).*rand(1,totalParticle);    % Initial x position
y_pos = (Width).*rand(1,totalParticle);    % Initlal Y position     




angle = 2*pi*rand(1,totalParticle);
Vx=(initial_Vth).*cos(angle); % velocity * random direction
Vy=(initial_Vth).*sin(angle);

VTot = sqrt(Vx.^2 + Vy.^2);
time = 0;
Ix = 0;

particle_Time = zeros(1,totalParticle);


% exponential scattering probability
Pscat = 1- exp(-(timeStep/t_mn));


VTot=sqrt(Vx.^2 + Vy.^2);


for n = 1: SimulationTime
    100*n/SimulationTime
    
    totalTime = timeStep*SimulationTime;

    % Update the electron position and velocity each time step
    prev_Ix = Ix;
    prev_time = time;
    time = prev_time + timeStep;
    
    %prev_tmn = tmn;

    LastX = x_pos;
    LastY = y_pos;

    
    Sx = Vx.*timeStep+ (1/2)*Acceleration*timeStep*timeStep;
    Vx = Vx + Acceleration*timeStep;
    Sy = Vy.*timeStep;
    
    x_pos = LastX + Sx;
    y_pos = LastY + Sy;
    

    
    LastX(x_pos<0) = Length;
    x_pos(x_pos<0)= x_pos(x_pos<0)+Length;

    
    LastX(x_pos>Length) = 0;
    x_pos(x_pos>Length)=x_pos(x_pos>Length)-Length;

    Vy(y_pos<=0) = -Vy(y_pos<=0);
    y_pos(y_pos<=0)=0;
    
    Vy(y_pos>=Width) = -Vy(y_pos>=Width);
    y_pos(y_pos>=Width)=Width;

    shouldScat =  Pscat>rand(n,1);
    particle_Time(Pscat>rand(n,1))=0;
    particle_Time = particle_Time + timeStep;
    Vx(shouldScat)=initial_Vth/sqrt(2)*randn(sum(shouldScat),1);%(1,sum(shouldScat));
    Vy(shouldScat)=initial_Vth/sqrt(2)*randn(sum(shouldScat),1);%(1,sum(shouldScat));
    VTot=sqrt(Vx.^2 + Vy.^2);


    figure(1)
    col = hsv(10);

    xlim([0 200e-9]);
    ylim([0 100e-9]);
    xlabel('Length (m)');
    ylabel('Width (m)');
    for k = 1:displayParticle
        plot([LastX(k) x_pos(k)],[LastY(k) y_pos(k)],'color',col(k,:));
        hold on
    end
   pause(0.01)
   figure(2)
   Ix = q_0*totalParticle*mean(Vx)*Width;
   plot([prev_time time],[prev_Ix Ix],'r');
   hold on
   title('Current vs Time');
   xlabel('Time (s)');
   ylabel('Current (A)');
end

figure(3)
electronMapX = linspace(0,Length, 101);
electronMapY = linspace(0,Width, 51);
MapX = linspace(0,Length, 100);
MapY = linspace(0,Width, 50);
ElectronDensity = histcounts2(y_pos,x_pos,electronMapY,electronMapX);

surf(MapX,MapY,ElectronDensity),title('Electron density map');
xlabel('Length (m)');
ylabel('Width (m)');
hold on

Temp = ((VTot(:).^2).* mass_eff)./(2*kb);
figure(4)
MapX = linspace(0,Length,totalParticle);
MapY = linspace(0,Width,totalParticle);
[X,Y] = meshgrid(MapX,MapY);
Tem = griddata(x_pos,y_pos,Temp,X,Y);
surf(MapX,MapY,Tem),title('Temperature map');


%surf(MapX,MapY,TemperarureMap),title('Temperature');
xlabel('Length (m)');
ylabel('Width (m)');
hold on


E = VoL/Length
F = E*q_0
Acceleration = F/(m_0)
