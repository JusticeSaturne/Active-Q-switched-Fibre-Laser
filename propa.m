function dpdz = propa(~,p)
% PROPA function represents the system of pump and laser propagation equations
% 
% Inputs:
% -The function takes the pump power and initial guesses of output powers
%  
% Outputs: 
%   dpdz(1): Represents the pump propagating power 
%   dpdz(2): Represents the laser propagating power
%  
% Comments:
%   - Computation is done in the ideal case of a two level energy system
%
% References: 
%  
% Written by Justice Sompo, University of Johannesburg, South Africa      
close all
clc
%
format longe
% Rate equation Parameters
sigma12S = 3.15e-25;   % Laser absorption cross section
sigma21S = 4.44e-25;   % Laser emission cross section
sigma12P = 1.8e-25;      % Pump absorption cross section
sigma21P = 0;          % Pump emission cross section
% Lifetimes
t21 = 10e-3;           % Lifetime of the metastable energy level of Erbium
A21 = 1/t21;
% Concentrations
%Ner = 1.2e25;          % Erbium ions concentration 
Ner = 8e24;
% Overlap factors and background losses
Gammap = 0.83;         % Overlap factor at pump wavelength
Gammas = 0.7;         % Overlap factor at signal wavelength
alphap = 0.005;        % Background loss at pump wavelength
alphas = 0.005;        % Background loss at laser wavelength
lambdap = 980e-9; % Pumb wavelenght
lambdas = 1550e-9; %signal wavelength 
% Other physical constants
r = 2.3e-6;           % Doped fibre core
h = 6.626e-34;        % Plank's constant
cel = 3e8;            % speed of light in free space
% Calculated parameters
A = pi*r^2;  
F_s = cel/lambdas;den1 = A*h*F_s;
F_p = cel/lambdap;  den2 = A*h*F_p;
% Fiber Bragg Grating Parameters
W12 = ((sigma12S*Gammas)*(p(2)))/den1;
W21 = ((sigma21S*Gammas)*(p(2)))/den1;
R12 = ((sigma12P*Gammap)*(p(1)))/den2;
%Rate equations solving
N1 = Ner*(W21+A21)./(W12+R12+W21+A21); % Population density of the ground energy level
N2 = Ner*(W12+R12)./(W21+A21+W12+R12); % Population density of themetastable energy level 
dpdz = zeros(2,1);
dpdz(1)= Gammap*(sigma21P*N2-sigma12P*N1)*p(1)-alphap*p(1);
dpdz(2)= Gammas*(sigma21S*N2-sigma12S*N1)*p(2)-alphas*p(2);

 