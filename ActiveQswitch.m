function [N1z,N2z,Ppz,Psz,Gainz,xc,N1,N2,Psf,Ppf,Gain,time,reflectivity ] =...
    ActiveQswitch(Pump_Power,L,Core_D,Rise_Time)
% Function ActiveQswitch implemente a finite difference method to compute
% output characteristics of an Erbium doped active Q-switched fibre laser
% using an acousto-optic modulator as switching element
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
%   - gain: gain distribution along the cavity
%   - Gain variation as a function of time
%   - ppf: Pump power distribution along the cavity
%   - psf: laser power distribution along the cavity
%   - Gain: Gain variation as a function of time
%   - Ppf: Residual pump
%   - Psf: Q-switched fibre laser output power
%   - xc: sections number along the cavity
%   - time: computation time
% Comments:
%   - Erbium doped gain medium is considered as a two level system and the
% resulting rate equations are solved with simple analytical method using
% substitution
%   -Initial conditions of the finite difference scheme are obtained by
%   solving the rate equations first in steady state
%
% References:
%
% Written by Justice Sompo, University of Johannesburg, South Africa
%
%
%

%==========================================================================
% Values for the variables in simulation
%==========================================================================
close all
clc
% Check inputs
if nargin <4
    Rise_Time = 40e-9;
    if nargin < 3
        Core_D = 2.3e-6;
        if nargin < 2
            L = 5;
            if nargin < 1
                error(message('MATLAB:LaserOutput:NotEnoughInputs'));
            end
        end
    end
end
lamddas = 1550e-9;         % Signal wavelength in meters
sigmaas = 3.15e-25;        % Absorption cross section at signal wavelength m^2
sigmaes = 4.44e-25;        % Emission cross section at signal wavelength m^2
NA = 0.2;
dcore = Core_D;            % Fiber core diameter
alfap = 0.005;             % Attenuation coefficient at pump wavelength
alfas = 0.005;             % Attenuation coefficient at signal wavelength
n0 = 4e24;                 % Doping concentration of the fiber from ions
c0 = 3e8;                  % Speed of light in vacuum
gamas = 0.7;
h = 6.626e-34;             % Plank constant
nfiber = 1.45;             % Refractive index of the fiber
vg = c0/nfiber;            % Group velocity in doped fiber
step_number = 50;
deltax = L/step_number;          % Longitudianl section in doped fiber
deltat = L/(step_number*vg);     % Sections in the time domain
deltalambda = 3e-6;
Aco = pi*(dcore/2)^2;       % Transverse section of the doped fiber
lambdap = 980e-9;           % Pump wavelength
gamap = 0.83;
sigmaep = 1.8e-25;          % Absorption cross section at pump wavelength
sigmaap = 2.44e-24;         % Emission cross section at pump wavelength
tao = 10e-3;                % Erbium excited state lifetime
Rmax = 0.9;
total_loss = 0.2032;        % 3 dB filter, 3 dB circulator, 0.8 dB WDM, 0.1 dB Coupler
leng = L;
coupler = Rmax;
Nt = n0;
sections = step_number+1;
round_trip_time = L/vg;
% compute steady state values distribution
pumppower = Pump_Power; % Pump power in watts
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
%==========================================================================
% Finite difference scheme implementation
%==========================================================================
ppf(step_number,1) = pumppower;
temps = deltat*1e6;
Time = deltat:deltat:temps;
t = 0;                    % very important to be used in the boundary conditions function
tr = Rise_Time;
to = 5e-6;
for k = 1:length(Time)
    ppf(1,2) = pumppower;
    for s = 2:step_number+1
        n2(s,2) = n2(s,1)+deltat*((gamap*lambdap/(h*c0*Aco))*(sigmaap*(n0-n2(s,1))...
            -sigmaep*n2(s,1))*ppf(s,1)+(gamas*lamddas/(h*c0*Aco))*(sigmaas*(n0-n2(s,1))...
            -sigmaes*n2(s,1))*psf(s,1)-n2(s,1)/tao);
        ppf(s,2) = ppf(s-1,1)+deltax*(gamap*(sigmaep*n2(s-1,1)-sigmaap*(n0-n2(s-1,1)))*ppf(s-1,1)...
            -alfap*ppf(s-1,1));
        psf(s,2) = psf(s-1,1)+deltax*(gamas*(sigmaes*n2(s-1,1)-sigmaas*(n0-n2(s-1,1)))*psf(s-1,1)...
            +2*gamas*sigmaes*n2(s-1,1)*h*c0*c0*deltalambda/(lamddas)^3-alfas*psf(s-1,1));
        n1(s,2) = n0-n2(s,1);
        gain(s,2) = gamas*(sigmaes*n2(s-1,1)-sigmaas*(n0-n2(s-1,1)));
    end
    n2(1,2) = n2(1,1)+deltat*((gamap*lambdap/(h*c0*Aco))*(sigmaap*(n0-n2(1,1))...
        -sigmaep*n2(1,1))*ppf(1,1)+(gamas*lamddas/(h*c0*Aco))*(sigmaas* (n0-n2(1,1))...
        -sigmaes*n2(1,1))*psf(1,1)-n2(1,1)/tao);
    n1(1,2) = n0-n2(1,1);
    % Boundary conditions)
    %====================
    t = t + deltat;
    R = feval('boundary',Rmax,tr,to,t);
    psf(1,2) = psf(step_number+1,1)*R*total_loss;
    for s = 1:step_number+1;
        n2(s,1) = n2(s,2);
        n1(s,1) = n1(s,2);
        ppf(s,1) = ppf(s,2);
        psf(s,1) = psf(s,2);
        gain(s,1) = gain(s,2);
        xc(s) = (s-1)*deltax;
        N2z(s) = n2(s,2);
        N1z(s) = n1(s,2);
        Ppz(s) = ppf(s,2);
        Psz(s) = psf(s,2);
        Gainz(s) = gain(s,2);
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
end

