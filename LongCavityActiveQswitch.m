function [N1,N2,Psf,Ppf,Gain,time] =...
    LongCavityActiveQswitch(Pump_Power,UndopedFL,Rise_Time,FiberLength,Core_D,Ner)
%==========================================================================
%Part I: Provide initial values for the variables in simulation
%==========================================================================
% Function LongCavityActiveQswitch implements a finite difference method to compute
% output characteristics of an Erbium doped active Q-switched fibre laser
% using an acousto-optic modulator as switching element. This function
% include a length of undoped fibre in the cavity
% Inputs:
%  - Pump_Power: Input pump power
%  - UndopedFL: Length of undoped fibre
%  - Rise_Time: Rise time of the acousto-optic modulator
%  - FiberLength:Length of the Erbium doped fibre
%  - Core_D: Core diameter of the doped fibre
%  - Ner: Erbium doped fibre concentration
%  
%
% Outputs:
%   - N1: Ground level population density variation with time
%   - N2: excited state population density variation with time
%   - Gain variation as a function of time
%   - Gain: Gain variation as a function of time
%   - Ppf: Residual pump
%   - Psf: Q-switched fibre laser output power
%   - time: computation time
% Comments:
%   - Erbium doped gain medium is considered as a two level system and the
%     resulting rate equations are solved with simple analytical method using
%     substitution
%   -Initial conditions of the finite difference scheme are obtained by
%     solving the rate equations first in steady state with a shooting
%     algorithm
%
% References:
%
% Written by Justice Sompo, University of Johannesburg, South Africa
%
%
%
close all
clc
% Check inputs
if nargin < 6
    Ner = 4e24;
    if nargin < 5
        Core_D = 2.3e-6;
        if nargin < 4
            FiberLength = 3;
            if nargin < 3
                Rise_Time = 40e-9;
                if nargin < 2
                    UndopedFL = 12;
                    if nargin < 1
                        error(message('MATLAB:LaserOutput:NotEnoughInputs'));
                    end
                end
            end
        end
    end
end
% Values for the variables in simulation
%==========================================================================
lambdas = 1550e-9;         % Signal wavelength
sigmaas = 3.15e-25;        % Absorption cross section at signal wavelength
sigmaes = 4.44e-25;        % Emission cross section at signal wavelength
NA = 0.2;                  
dcore = Core_D;            % Fiber core diameter
L = FiberLength;           % Fiber length
L1 = UndopedFL;            % Length of undoped fiber
alfap = 0.005;             % attenuation coefficient at pump wavelength
alfas = 0.005;             % Attenuation coefficient at signal wavelength
n0 = Ner;                  % Doping concentration of the fiber
c0 = 3e8;                  % Speed of light in vacuum
gamas = 0.7;
h = 6.626e-34;             % Plank constant
nfiber = 1.45;             % Refractive index of the fiber
vg = c0/nfiber;            % Group velocity in doped fiber
step_number = 50;
deltax = L/step_number;          % Longitudianl section in doped fiber
deltat = L/(step_number*vg);     % Sections in the time domain
step_number_1 = 100;
deltax1 = L1/step_number_1;
deltat1 = L1/(step_number*vg);
deltalamda = 3e-6;
Aco = pi*(dcore/2)^2;       % Area of the doped fiber core
lambdap = 980e-9;           % Pump wavelength
gamap = 0.83;
sigmaep = 1.8e-25;          % Absorption cross section at pump wavelength
sigmaap = 2.44e-24;         % Emission cross section at pump wavelength
tao = 10e-3;                % Erbium excited metastable state lifetime
Rmax = 0.9;
% Adding this to generate initial condition
leng = L;
coupler = Rmax;
Nt = n0;
sections = step_number+1;
% Compute steady state values distribution
pumppower = Pump_Power;         % Pump power in watts
Pump = pumppower;
[PP,Ps,N2x] = LaserOutput(Pump,leng,coupler,Nt,sections);
% Initial conditions for doped fiber
%--------------------------------------------------------------------------
s = 1:step_number+1;
n2(s,1) = N2x;
n1(s,1) = n0-N2x;
ppf(s,1) = PP;
psf(s,1) = Ps;
gain(s,1) = 0.0;
xc(s,1) = 0.0;
%Initial conditions for undoped fiber
%--------------------------------------------------------------------------
s1 = 1:step_number_1+1;
ppfun(s1,1) = 0;
psfun(s1,1) = 0;
%==========================================================================
%Part II: Iterative part
%==========================================================================

ppf(1,1) = pumppower;
temps = deltat*400000;
Time = deltat:deltat:temps;
t = 0;   % very important to be used in the boundary conditions function
tr = Rise_Time;
to = 5e-6;
longueur = 0;
for k = 1:length(Time)
    ppf(1,2) = pumppower;
    for s = 2:step_number+1
        n2(s,2) = n2(s,1)+deltat*((gamap*lambdap/(h*c0*Aco))*(sigmaap*(n0-n2(s,1))...
            -sigmaep*n2(s,1))*ppf(s,1)+(gamas*lambdas/(h*c0*Aco))*(sigmaas*(n0-n2(s,1))...
            -sigmaes*n2(s,1))*psf(s,1)-n2(s,1)/tao);
        ppf(s,2) = ppf(s-1,1)+deltax*(gamap*(sigmaep*n2(s-1,1)...
            -sigmaap*(n0-n2(s-1,1)))*ppf(s-1,1)-alfap*ppf(s-1,1));
        psf(s,2) = psf(s-1,1)+deltax*(gamas*(sigmaes*n2(s-1,1)...
            -sigmaas*(n0-n2(s-1,1)))*psf(s-1,1)...
            +2*gamas*sigmaes*n2(s-1,1)*h*c0*c0* deltalamda/(lambdas)^3)-alfas*psf(s-1,1);
        n1(s,2) = n0-n2(s,1);
        gain(s,2) = gamas*(sigmaes*n2(s-1,1)-sigmaas*(n0-n2(s-1,1)));
    end
    n2(1,2) = n2(1,1)+deltat*((gamap*lambdap/(h*c0*Aco))*(sigmaap*(n0-n2(1,1))...
        -sigmaep*n2(1,1))*ppf(1,1)+(gamas*lambdas/(h*c0*Aco))*(sigmaas* (n0-n2(1,1))...
        -sigmaes*n2(1,1))*psf(1,1)-n2(1,1)/tao);
    n1(1,2) = n0-n2(1,1);
    
    ppfun(1,1) = ppf(step_number+1,2);
    psfun(1,1) = psf(step_number+1,2);
    for s1 = 2:step_number_1+1
        longueur = longueur+deltax1;
        psfun(s1,2) = psfun(s1-1,1)+deltax1*(-alfas*psfun(s1-1,1));
    end
    % Boundary conditions)
    %====================
    t = t + deltat;
    R = feval('boundary',Rmax,tr,to,t);
    psf(1,2) = psfun(step_number_1+1,2)*R;
    % Updating initial conditions
    %============================
    for s1 = 1:step_number_1+1
        psfun(s1,1) = psfun(s1,2);
    end
    for s = 1:step_number+1;
        n2(s,1) = n2(s,2);
        n1(s,1) = n1(s,2);
        ppf(s,1) = ppf(s,2);
        psf(s,1) = psf(s,2);
        gain(s,1) = gain(s,2);
        xc(s) = (s-1)*deltax;
    end
    N2(k) = n2(step_number+1);
    N1(k) = n1(step_number+1);
    Psf(k) = psf(step_number+1)*(1-R);
    Ppf(k) = ppf(step_number+1);
    Gain(k) = gain(step_number+1);
    z(k) = k;
    time(k) = t*1e4;
    reflectivity(k) = R;
end
figure
subplot(2,2,1)
plot(time,Psf,'Linewidth',2);
xlabel('Time (Microseconds)')
ylabel('Pulse Power (W)')
subplot(2,2,2)
plot(time,Gain,'r','Linewidth',2)
xlabel('Time (Microseconds)')
ylabel('Gain (m^-1)')
subplot(2,2,3)
plot(time,N2,'m','Linewidth',2)
xlabel('Time (Microseconds)')
ylabel('N2 (ions/m^3)')
subplot(2,2,4)
plot(time,N1,'k','Linewidth',2)
xlabel('Time (Microseconds)')
ylabel('N1 (ions/m^3)')
toc