function [PP,Ps,N2,N1,N,t] = LaserOutput(Pump,leng,coupler,Nt,sections)
% This function solves the rate equations and propagation
% equations for a ring cavity fiber laser using the MATLAB ODE 45 routine
% and a shooting algorithm
%
%
% Inputs:
%  -Pump: Pump power
%  -leng: Length of the doped fibre
%  -coupler: Coupling ratio of the output coupler
%
% Outputs:
%  -PP : Pump power distribition
%  -Ps: Laser power
%  -N2: Population density at the metastable level
%  -N1: Population density at ground level
%  -N : Total population density of the doped fibre
%
% Comments:
%  - Computation is done in the ideal case of a two level energy system
%  equations analytically.

% References:
%
%
% Written by Justice Sompo, University of Johannesburg, South Africa
close all
clc
format long
% Check inputs
if nargin < 5
    sections = 200;
    if nargin <4
        %     Nt = 1.26e25;
        Nt = 8e24;
        if nargin < 3
            coupler = 0.1;
            if nargin < 2
                leng = 5;
                if nargin < 1
                    error(message('MATLAB:LaserOutput:NotEnoughInputs'));
                end
            end
        end
    end
end
%==========================================================================
%Rate equations parameters
%==========================================================================
sigma12S = 3.15e-25;
sigma21S = 4.44e-25;
sigma12P = 1.8e-25;
% Lifetimes
t21 = 10e-3;
A21 = 1/t21;
% Concentrations
Ner = Nt;
% Overlap factors and background losses
Gammap = 0.83;
Gammas = 0.7;
lambdap = 980e-9; % Pumb wavelength
lambdas = 1550e-9; %signal wavelength
% Other physical constants
r = 2.3e-6;
h = 6.626e-34;
cel = 3e8;
% Calculated parameters
A = pi*r^2;
F_s = cel/lambdas;
den1 = A*h*F_s;
F_p = cel/lambdap;
den2 = A*h*F_p;
%==========================================================================
% Implementation of the shooting secant algorithm to solve propagation
% equations in steady state
%==========================================================================
epsilon = 1e-6; % tolerance for computation
count = 0; % iniltialise iterations before the while loop
RL = coupler;  % reflectivity of mirror at L.
Pp = Pump;
% Default values of initial guess
% x1 = 10e-3; % first signal power guess
% x2 = 100e-3; % second signal power guess
x1 = 10e-3; % first signal power guess
x2 = 100e-3; % second signal power guess
L = leng;
%tspan = [0 6];
tspan = linspace(0,L,sections); % to obtain solutions at each position for 100 positions one can write linspace(0,L,100)
[t1,y1]=ode45(@propa,tspan,[Pp x1]); %y1 will have the values of Pp in its first column and the values of Ps in its second column.
[t2,y2]=ode45(@propa,tspan,[Pp x2]);
i = 1;
Psl1=y1(end,2); %  Ps1 power at z = L
Psl2=y2(end,2); %  Ps2 power at z = L
P01 = Psl1*RL;
P02 = Psl2*RL;
m1 = P01-x1;
m2 = P02-x2;
while (abs(x2-x1) > epsilon)
    tmp=x2;
    count = count+1;
    x2 = x1-(x1-x2)/(m1-m2)*m1;
    x1 = tmp;
    [t1,y1]=ode45(@propa,tspan,[Pp x1]);
    [t2,y2]=ode45(@propa,tspan,[Pp x2]);
    Psl1=y1(end,2); %  Ps power at z = L
    Psl2=y2(end,2); %  Ps power at z = L
    P01 = Psl1*RL;
    P02 = Psl2*RL;
    m1 = (P01-x1);
    m2 = (P02-x2);
    i = i+1;
    t = t2;
    PP = y2(:,1); % Pump power distribution
    Ps = y2(:,2); % Signal power distribution
    %Population distribution
    W12 = ((sigma12S*Gammas)*Ps)/den1;
    W21 = ((sigma21S*Gammas)*Ps)/den1;
    R12 = ((sigma12P*Gammap)*PP)/den2;
    N1 = Ner*(W21+A21)./(W12+R12+W21+A21);
    N2 = Ner*(W12+R12)./(W21+A21+W12+R12);
    N = N2./Ner;
end
%==========================================================================
% Plotting
%==========================================================================
% figure(1)
% plot(t,N2,'m','linewidth',2)
% hold on
% plot(t,N1,'g','linewidth',2)
% figure(2)
% plot(t,Ps,'r-','linewidth',2)
% hold on
% plot(t,PP,'linewidth',2)

