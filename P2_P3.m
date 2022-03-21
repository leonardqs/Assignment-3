clear
close all


VoL = 0.8; % Left side of the Area has Boundary Voltage = VoL
VoR = 0; % Right side of the Area has Boundary Voltage = VoR
VoT = 0; % Top side of the Area has Boundary Voltage = VoT
VoB = 0; % Bottom side of the Area has Boundary Voltage = VoB
Pixel = 5e8; % Number of mesh per unit length/width

Lb = 0.4e-7; % Box Length
Wb = 0.4e-7; % Box Width
xBox = Lb*Pixel; % Box length Pixel Number
yBox = Wb*Pixel; % Box width Pixel Number
Length = 2e-7; % Area Length
Width = 1e-7; % Area Width
nx = Length*Pixel; % Area length Pixel Number
ny = Width*Pixel; % Area wdith Pixel Number
Sigma = 1; % Outside Box Area Conductivity
BoxSigma = 0.01; % Inside Box Area Conductivity
G = sparse(nx*ny); % G matrix has size(nx*ny,nx*ny)
B = zeros(nx*ny,1); % B is the product of G matrix * V
Conductivity = Sigma*ones(nx,ny); % Conductivity of the entire area

q_0 = 1.60217653e-19;             % electron charge
hb = 1.054571596e-34;             % Dirac constant
h = hb * 2 * pi;                % Planck constant
m_0 = 9.10938215e-31;             % electron mass
kb = 1.3806504e-23;               % Boltzmann constant
eps_0 = 8.854187817e-12;          % vacuum permittivity
mu_0 = 1.2566370614e-6;           % vacuum permeability
c = 299792458;                    % speed of light
g = 9.80665;                      % metres (32.1740 ft) per s²



timeStep = 6e-15;                     % Time Step
totalParticle = 1e4;
displayParticle = 10;              % Display Particle in the box
SimulationTime = 1000;                       % Simulation Time
T = 300;                            % Default Temperature
t_mn = 0.2e-12;                     % Mean time between collision 
loop_index = 10000;
time = 0;




% exponential scattering probability
Pscat = 1- exp(-(timeStep/t_mn));



for iRow = 1:nx
    for jColumn = 1:ny

        if iRow>=(nx-xBox)/2 && iRow<=((nx+xBox)/2) && jColumn<=yBox
            Conductivity(iRow,jColumn) = BoxSigma;
        elseif iRow>=(nx-xBox)/2 && iRow<=((nx+xBox)/2) && jColumn<=ny && jColumn>ny-yBox
            Conductivity(iRow,jColumn) = BoxSigma;
        end
    end
end



% for jColumn = 1:ny
%     for iRow = 1:nx
for iRow = 1:nx
    for jColumn = 1:ny
        n = jColumn+(iRow-1)*ny;
        % Left side Boundary Condition
       
        if iRow == 1      
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = VoL;
        % Right side Boundary Condition
           
        elseif iRow == nx
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = VoR;

        % Bottom side Boundary Condition
        elseif jColumn == 1    
            nxm = jColumn+((iRow-1)-1)*ny;
            nxp = jColumn+((iRow+1)-1)*ny;
            nyp = (jColumn+1)+(iRow-1)*ny;

            rxm = (Conductivity(iRow,jColumn) + Conductivity(iRow-1,jColumn))/2;
            rxp = (Conductivity(iRow,jColumn) + Conductivity(iRow+1,jColumn))/2;
            ryp = (Conductivity(iRow,jColumn) + Conductivity(iRow,jColumn+1))/2;
            
            G(n,n) = -(rxm+rxp+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nyp) = ryp;
            
        % Top side Boundary Condition
        elseif jColumn == ny 
            nxm = jColumn+((iRow-1)-1)*ny;
            nxp = jColumn+((iRow+1)-1)*ny;
            nym = (jColumn-1)+(iRow-1)*ny;

            rxm = (Conductivity(iRow,jColumn) + Conductivity(iRow-1,jColumn))/2;
            rxp = (Conductivity(iRow,jColumn) + Conductivity(iRow+1,jColumn))/2;
            rym = (Conductivity(iRow,jColumn) + Conductivity(iRow,jColumn-1))/2;
            
            G(n,n) = -(rxm+rxp+rym);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;

        else 
            nxm = jColumn+((iRow-1)-1)*ny;
            nxp = jColumn+((iRow+1)-1)*ny;
            nym = (jColumn-1)+(iRow-1)*ny;
            nyp = (jColumn+1)+(iRow-1)*ny;

            rxm = (Conductivity(iRow,jColumn) + Conductivity(iRow-1,jColumn))/2;
            rxp = (Conductivity(iRow,jColumn) + Conductivity(iRow+1,jColumn))/2;
            rym = (Conductivity(iRow,jColumn) + Conductivity(iRow,jColumn-1))/2;
            ryp = (Conductivity(iRow,jColumn) + Conductivity(iRow,jColumn+1))/2;
            
            G(n,n) = -(rxm+rxp+rym+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp;
        end

    end

end


% 3D Plot of Voltage V(x,y)
figure('name','Voltage V(x,y)')

Vn = G\B; % Find (ny*ny:1) size, V = G\B
% Mapping the V to the matrix size of nx*ny
for iRow = 1: nx
    for jColumn = 1:ny
         n = jColumn+(iRow-1)*ny;
         V(iRow,jColumn) = Vn(n);    
    end
end

x = linspace(0,Length,nx);
y = linspace(0,Width,ny);
[X,Y] = meshgrid(x,y);

surf(x,y,V')
view(31,23)
xlim([0 Length])
ylim([0 Width])
title("3D Plot of Voltage V(x,y)")
xlabel("Length")
ylabel("Width")
zlabel("Voltage (V)")

% 3D Plot of Electric Field E(x,y)
[Ex,Ey] = gradient(V');
figure('name','Electric Field E(x,y)')
quiver(X,Y,-Ex,-Ey)
xlim([0 Length])
ylim([0 Width])
title("3D Plot of Electric Field E(x,y)")
xlabel("Length")
ylabel("Width")
zlabel("Electric Field (V/m)")

figure('name','Electric Field E(x)')
surf(X,Y,-Ex)
xlim([0 Length])
ylim([0 Width])
title("3D Plot of Electric Field E(x)")
xlabel("Length")
ylabel("Width")
zlabel("Electric Field (V/m)")
view(2)

figure('name','Electric Field E(y)')
surf(X,Y,-Ey)
xlim([0 Length])
ylim([0 Width])
title("3D Plot of Electric Field E(y)")
xlabel("Length")
ylabel("Width")
zlabel("Electric Field (V/m)")
view(2)

%E = sqrt((Ex.^2)+ (Ey.^2));

% 3D Plot of Conductivity Cond(x,y)
figure('name','Conductivity Cond(x,y)')
cond = Conductivity';
surf(X,Y,cond)
xlim([0 Length])
ylim([0 Width])
title("3D Plot of Conductivity Cond(x,y)")
xlabel("Length")
ylabel("Width")
zlabel("Conductivity (ohm*m)")

% 3D Plot of Current Density J(x,y)
figure('name','Current density J(x,y)')
Jx = cond.*(-Ex);
Jy = cond.*(-Ey);


quiver(X,Y,Jx,Jy)
xlim([0 Length])
ylim([0 Width])
title("3D Plot of Current Density J(x,y)")
xlabel("Length")
ylabel("Width")
zlabel("Current Density (A/Area)")



mass_eff =  0.26*m_0;
initial_Vth = sqrt((2* kb *T)/mass_eff);
particle_Time = zeros(1,totalParticle);


x_pos = zeros(1,totalParticle);
y_pos = zeros(1,totalParticle);



index = 1;
for k = 1:loop_index
    rand_x = 2*rand;
    rand_y = rand;
    if(rand_x < 0.8 || rand_x > 1.2) || ((rand_y > 0.4 && rand_y < 0.6) && (rand_x >= 0.8 || rand_x <= 1.2))
        x_pos(1,index) = 1e-7*rand_x;
        y_pos(1,index) = Width*rand_y;
        index = index + 1;
    end
    if index == totalParticle + 1
        break
    end
end    
% Generate the initial electron position
% x_pos = (Length).*rand(1,totalParticle);    % Initial x position
% y_pos = (Width).*rand(1,totalParticle);    % Initlal Y position     



 Vx = initial_Vth.*randn(1,totalParticle);
 Vy = initial_Vth.*randn(1,totalParticle);
% angle = 2*pi*rand(1,totalParticle);
% Vx=(initial_Vth).*cos(angle); % velocity * random direction
% Vy=(initial_Vth).*sin(angle);

% VoL = 0.8; % Left side of the Area has Boundary Voltage = VoL
% E = VoL/Length
% F = E*q_0
% Acceleration = F/m_0




xm = linspace(0,Length,100);
ym = linspace(0,Width,50);
[XM,YM] = meshgrid(xm,ym);

figure(7)
rectangle('position',[0.8e-7 0e-7 0.4e-7 0.4e-7]);
rectangle('position',[0.8e-7 0.6e-7 0.4e-7 0.4e-7]);
hold on
col = hsv(10);
xlim([0 200e-9]);
ylim([0 100e-9]);
xlabel('Length (m)');
ylabel('Width (m)');
for k = 1: SimulationTime
    100*k/SimulationTime
    
    totalTime = timeStep*SimulationTime;
    % Update the electron position and velocity each time step
    prev_time = time;
    time = prev_time + timeStep;
    %prev_tmn = tmn;
    LastX = x_pos;
    LastY = y_pos;
%    
    Ex_p = interp2(XM,YM,-Ex,x_pos,y_pos);
    Ey_p = interp2(XM,YM,-Ey,x_pos,y_pos);
% E = VoL/Length
% F = E*q_0
% Acceleration = F/(m_0)
    XAcceleration = (Ex_p*q_0)/m_0;
    YAcceleration = (Ey_p*q_0)/m_0;

    Sx = Vx.*timeStep+ (1/2)*timeStep*timeStep*XAcceleration;%*8.564e+17;
    Vx = Vx + timeStep*XAcceleration;%*8.564e+17;
    Sy = Vy.*timeStep + (1/2)*timeStep*timeStep*YAcceleration;%*8.564e+7
    Vy = Vy + timeStep*YAcceleration;%*8.564e+7
    
    x_pos = LastX + Sx;
    y_pos = LastY + Sy;
    
    

    VTot=sqrt(Vx.^2 + Vy.^2);
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
    for n = 1:totalParticle

      
        if (x_pos(1,n) >= 0.8e-7 && x_pos(1,n) <= 1e-7 ) && (y_pos(1,n) <= 0.4e-7 || y_pos(1,n) >= 0.6e-7) && LastX(1,n) < 0.8e-7
            x_pos(1,n) = 0.8e-7;
            LastX(1,n) = 0.8e-7;
            Vx(1,n) = -(Vx(1,n));
    
        elseif (x_pos(1,n) >= 1e-7 && x_pos(1,n) <= 1.2e-7) && (y_pos(1,n) <= 0.4e-7 || y_pos(1,n) >= 0.6e-7) && LastX(1,n) > 1.2e-7  
            x_pos(1,n) = 1.2e-7;
            LastX(1,n) = 1.2e-7;
            Vx(1,n) = -(Vx(1,n));
    
        elseif (y_pos(1,n) <= 0.4e-7) && (x_pos(1,n) >= 0.8e-7 && x_pos(1,n) <= 1.2e-7)
            y_pos(1,n) = 0.4e-7;
            LastY(1,n) = 0.4e-7;
            Vy(1,n) = -(Vy(1,n));
        elseif (y_pos(1,n) >= 0.6e-7) && (x_pos(1,n) >= 0.8e-7 && x_pos(1,n) <= 1.2e-7)
            y_pos(1,n) = 0.6e-7;
            Vy(1,n) = -(Vy(1,n));
        end
    
    end
    for p = 1:displayParticle
        plot([LastX(p) x_pos(p)],[LastY(p) y_pos(p)],'color',col(p,:));
        hold on
    end
    
    pause(0.01)
end

figure(8)
electronMapX = linspace(0,Length, 101);
electronMapY = linspace(0,Width, 51);
MapX = linspace(0,Length, 100);
MapY = linspace(0,Width, 50);
ElectronDensity = histcounts2(y_pos,x_pos,electronMapY,electronMapX);

surf(MapX,MapY,ElectronDensity),title('Electron density map');
xlabel('Length (m)');
ylabel('Width (m)');
hold on
