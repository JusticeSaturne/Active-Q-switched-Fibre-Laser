%==========================================================================
% Values for the variables in simulation
%==========================================================================
tic
clear all
close all
clc
lambdas = 1550e-9;         % signal wavelength
sigmaas = 3.15e-25;        % Absorption cross section at signal wavelength from MORASSE
sigmaes = 4.44e-25;        % Emission cross section at signal wavelength from MORASSE
NA = 0.2;                  % from MORASSE
dcore = 2.3e-6;            % fiber core diameter
L = 3;                     % Fiber length
alfap = 0.005;             % attenuation coefficient at pump wavelength
alfas = 0.005;             % Attenuation coefficient at signal wavelength
n0 = 4e24;                 % Doping concentration of the fiber from MORASSE
c0 = 3e8;                  % Speed of light in vacuum
gamas = 0.7;
h = 6.626e-34;             % Plank constant
nfiber = 1.45;             % Refractive index of the fiber
vg = c0/nfiber;            % Group velocity in doped fiber
step_number = 50;
deltax = L/step_number;         % longitudianl section in doped fiber
deltat = L/(step_number*vg);     % sections in the time domain
deltalambda = 3e-6;
Aco = pi*(dcore/2)^2;      % transverse section of the doped fiber
lambdap = 980e-9;           % pump wavelength
gamap = 0.83;
sigmaep = 1.8e-25;        % Absorption cross section at pump wavelength
sigmaap = 2.44e-24;         % Emission cross section at pump wavelength
tao = 10e-3;              % Erbium metastable  state lifetime
Rmax = 0.9;
total_loss = 0.2032;      % 3 dB filter, 3 dB circulator, 0.8 dB WDM, 0.1 dB Coupler
leng = L;
coupler = Rmax;
Nt = n0;
sections = step_number+1;
round_trip_time = L/vg;
% compute steady state values distribution
pumppower = 80e-3;         % Pump power in watt
Pump = pumppower;
[PP,Ps,N2x,N1x,N,tt] = LaserOutput(Pump,leng,coupler,Nt,sections);
% Initial conditions for doped fiber
%--------------------------------------------------------------------------
s = 1:step_number+1;
n2(s,1) = N2x;
n1(s,1) = n0-N2x;
ppf(s,1) = PP;
psf(s,1) = Ps;
gain(s,1) = 0.0;
xc(s,1) = 0.0;
%==========================================================================
% Finite difference scheme implementation
%==========================================================================

ppf(step_number,1) = pumppower; 
temps = deltat*1e6;
Time = deltat:deltat:temps;
t = 0;   % very important to be used in the boundary conditions function
tr = 40e-9;
to = 2e-6;
ti = 50e-6;
rep_rate = 35e3;
per = 1/rep_rate;
for k = 1:length(Time)
    ppf(1,2) = pumppower;
    for s = 2:step_number+1
        n2(s,2) = n2(s,1)+deltat*((gamap*lambdap/(h*c0*Aco))*(sigmaap*(n0-n2(s,1))-sigmaep*n2(s,1))*ppf(s,1)+(gamas*lambdas/(h*c0*Aco))*(sigmaas*(n0-n2(s,1))-sigmaes*n2(s,1))*psf(s,1)-n2(s,1)/tao);
        ppf(s,2) = ppf(s-1,1)+deltax*(gamap*(sigmaep*n2(s-1,1)-sigmaap*(n0-n2(s-1,1)))*ppf(s-1,1)-alfap*ppf(s-1,1));
        psf(s,2) = psf(s-1,1)+deltax*(gamas*(sigmaes*n2(s-1,1)-sigmaas*(n0-n2(s-1,1)))*psf(s-1,1)+2*gamas*sigmaes*n2(s-1,1)*h*c0*c0*deltalambda/(lambdas)^3-alfas*psf(s-1,1));
        n1(s,2) = n0-n2(s,1);
        gain(s,2) = gamas*(sigmaes*n2(s-1,1)-sigmaas*(n0-n2(s-1,1)));
    end
    n2(1,2) = n2(1,1)+deltat*((gamap*lambdap/(h*c0*Aco))*(sigmaap*(n0-n2(1,1))-sigmaep*n2(1,1))*ppf(1,1)+(gamas*lambdas/(h*c0*Aco))*(sigmaas* (n0-n2(1,1))-sigmaes*n2(1,1))*psf(1,1)-n2(1,1)/tao);
    n1(1,2) = n0-n2(1,1); % on complete le vecteur n2 dans la colonne 2 avec la derniere valeur SecNo+1
    % Boundary conditions)
    %====================
    t = t + deltat;
    %R = feval('boundary',Rmax,tr,to,t);
    if (t < ti)
        R = 0;
    elseif (t > ti && t < ti+tr)
        R = (Rmax/tr)*(t-ti);
    elseif (t > (ti+tr) && t < (ti+tr+to)) || (t > per+ti+tr && t < per+ti+tr+to) || (t > 2*per+ti+tr && t < 2*per+ti+tr+to) || (t > 3*per+ti+tr && t < 3*per+ti+tr+to)
        R = Rmax;
    elseif (t > ti+per && t < ti+per+tr)
        R = (Rmax/tr)*(t-(ti+per));
    elseif (t > ti+2*per && t < ti+2*per+tr)
        R = (Rmax/tr)*(t-(ti+2*per));
    elseif (t > ti+3*per && t < ti+3*per+tr)
        R = (Rmax/tr)*(t-(ti+3*per));
    else
        R = 0;
    end
    psf(1,2) = psf(step_number+1,1)*R*total_loss;
    for s = 1:step_number+1;
        n2(s,1) = n2(s,2);
        n1(s,1) = n1(s,2);
        ppf(s,1) = ppf(s,2);
        psf(s,1) = psf(s,2);
        gain(s,1) = gain(s,2);
        xc(s) = (s-1)*deltax;
    end
    N2(k) = n2(step_number+1)/10e23;
    N1(k) = n1(step_number+1)/10e23;
    Psf(k) = psf(step_number+1)*(1-R);
    Ppf(k) = ppf(step_number+1);
    Gain(k) = gain(step_number+1)*5;
    z(k) = k;
    time(k) = t*1e4;
    reflectivity(k) = R;
end
%plot(qq,N2,'g');
%  hold on
% figure(1)
plot(time,Psf);
hold on
plot(time,Gain,'r')
% figure(2)
%plot(time,reflectivity,'r');
% hold on
% plot(xc,n1,'m');
% plot(TimeSec:TimeSec:tt*TimeSec,N1,'k-');
%hold on
%plot(z,N2,'r-');
% plot(TimeSec:TimeSec:tt*TimeSec,Psf);
% hold on
% plot(TimeSec:TimeSec:tt*TimeSec,Gain,'r-');
% axis([0.00145 0.00155 0 6])
pic = max(Psf);
% plot(xc,psp)
% hold on
% plot(xc,psm)
% plot(xc,ppm)
%plot(xc,psp)
%hold on

toc