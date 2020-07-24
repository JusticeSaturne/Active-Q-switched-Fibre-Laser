function [N1,N2,Psf,Ppf,Gain,time] =...
    MultiPulseActiveQswitch(Pump_Power,RepRate,Rise_Time,FiberLength,Core_D,Ner)

% Function MultipleActiveQswitch implemente a finite difference method to compute
% output characteristics of an Erbium doped active Q-switched fibre laser
% using an acousto-optic modulator as switching element. A pulse train
% which number depends on the repetition frequency is obtained
%
% Inputs:
%  - Pump_Power: Input pump power
%  - L:Length of the Erbium doped fibre
%  - Core_D: Core diameter of the doped fibre
%  - Rise_Time: Rise time of the acousto-optic modulator
%
% Outputs:
%   - n1: Ground level population density distribution along the cavity
%   - n2: Excited state population density distribution along the cavity
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
%   solving the rate equations first in steady state with a shooting
%   algorithm
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
                    RepRate = 20e3;
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
lambdas = 1550e-9;         % signal wavelength
sigmaas = 3.15e-25;        % Absorption cross section at signal wavelength
sigmaes = 4.44e-25;        % Emission cross section at signal wavelength
NA = 0.2;                  
dcore = Core_D;            % fiber core diameter
L = FiberLength;           % Fiber length
alfap = 0.005;             % Attenuation coefficient at pump wavelength
alfas = 0.005;             % Attenuation coefficient at signal wavelength
n0 = Ner;                  % Doping concentration of the fiber from MORASSE
c0 = 3e8;                  % Speed of light in vacuum
gamas = 0.7;
h = 6.626e-34;             % Plank constant
nfiber = 1.45;             % Refractive index of the fiber
vg = c0/nfiber;            % Group velocity in doped fiber
step_number = 50;          % Number of sections
deltax = L/step_number;          % Longitudianl section in doped fiber
deltat = L/(step_number*vg);     % Sections in the time domain
deltalambda = 3e-6;
Aco = pi*(dcore/2)^2;       % Transverse section of the doped fiber
lambdap = 980e-9;           % Pump wavelength
gamap = 0.83;
sigmaep = 1.8e-25;          % Absorption cross section at pump wavelength
sigmaap = 2.44e-24;         % Emission cross section at pump wavelength
tao = 10e-3;                % Erbium metastable  state lifetime
Rmax = 0.9;
total_loss = 0.2032;        % 3 dB filter, 3 dB circulator, 0.8 dB WDM, 0.1 dB Coupler
leng = L;
coupler = Rmax;
Nt = n0;
sections = step_number+1;
round_trip_time = L/vg;
% compute steady state values distribution
pumppower = Pump_Power;         % Pump power in watt
Pump = pumppower;
[PP,Ps,N2x] = LaserOutput(Pump,leng,coupler,Nt,sections); % The function
% LaserOuput is called to solve the rate equations in steady state and
% provide initial conditions for finite difference scheme
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
tr = Rise_Time;
to = 2e-6;
ti = 50e-6;
rep_rate = RepRate;
per = 1/rep_rate;
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
            +2*gamas*sigmaes*n2(s-1,1)*h*c0*c0*deltalambda/(lambdas)^3-alfas*psf(s-1,1));
        n1(s,2) = n0-n2(s,1);
        gain(s,2) = gamas*(sigmaes*n2(s-1,1)-sigmaas*(n0-n2(s-1,1)));
    end
    n2(1,2) = n2(1,1)+deltat*((gamap*lambdap/(h*c0*Aco))*(sigmaap*(n0-n2(1,1))...
        -sigmaep*n2(1,1))*ppf(1,1)+(gamas*lambdas/(h*c0*Aco))*(sigmaas* (n0-n2(1,1))...
        -sigmaes*n2(1,1))*psf(1,1)-n2(1,1)/tao);
    n1(1,2) = n0-n2(1,1);
    % Boundary conditions)
    %====================
    t = t + deltat;
    if (t < ti)
        R = 0;
    elseif (t > ti && t < ti+tr)
        R = (Rmax/tr)*(t-ti);
    elseif (t > (ti+tr) && t < (ti+tr+to)) || (t > per+ti+tr && t < per+ti+tr+to)...
            || (t > 2*per+ti+tr && t < 2*per+ti+tr+to) ||...
            (t > 3*per+ti+tr && t < 3*per+ti+tr+to)
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
    N2(k) = n2(step_number+1);
    N1(k) = n1(step_number+1);
    Psf(k) = psf(step_number+1)*(1-R);
    Ppf(k) = ppf(step_number+1);
    Gain(k) = gain(step_number+1);
    z(k) = k;
    time(k) = t*1e4;
    reflectivity(k) = R;
end
% figure
% subplot(2,2,1)
% plot(time,Psf,'Linewidth',2);
% xlabel('Time (Microseconds)')
% ylabel('Pulse Power (W)')
% subplot(2,2,2)
% plot(time,Gain,'r','Linewidth',2)
% xlabel('Time (Microseconds)')
% ylabel('Gain (m^-1)')
% subplot(2,2,3)
% plot(time,N2,'m','Linewidth',2)
% hold on
% plot(time,N1,'k','Linewidth',2)
% xlabel('Time (Microseconds)')
% ylabel('N1,N2 (ions/m^3)')
% subplot(2,2,4)
% plot(time,Ppf,'g','Linewidth',2)
% xlabel('Time (Microseconds)')
% ylabel('Pump Power (W)')
