% This script plot differents values resulting from the main Active Q
% switched fibre laser function
tic
clear all
close all
clc
Pump_Power0 = 40e-3;
[n10,n20,ppf0,psf0,gain0,xc0,N10,N20,Psf0,Ppf0,Gain0,time0,reflectivity0 ]...
    = ActiveQswitch(Pump_Power0);

Pump_Power1 = 80e-3;
[n11,n21,ppf1,psf1,gain1,xc1,N11,N21,Psf1,Ppf1,Gain1,time1,reflectivity1 ]...
    = ActiveQswitch(Pump_Power1);

Pump_Power2 = 120e-3;
[n12,n22,ppf2,psf2,gain2,xc2,N12,N22,Psf2,Ppf2,Gain2,time2,reflectivity2 ]...
    = ActiveQswitch(Pump_Power2);

Pump_Power3 = 160e-3;
[n13,n23,ppf3,psf3,gain3,xc3,N13,N23,Psf3,Ppf3,Gain3,time3,reflectivity3 ]...
    = ActiveQswitch(Pump_Power3);
figure(1)
subplot(2,2,1)
plot(xc0,ppf0,'b','Linewidth',2)
xlabel('Position (m)')
ylabel('Pump Power (W)')
subplot(2,2,2)
plot(xc0,psf0,'Linewidth',2)
xlabel('Position (m)')
ylabel('Laser Power (W)')
subplot(2,2,3)
plot(xc0,n10,'m','Linewidth',2)
hold on
plot(xc0,n20,'r','Linewidth',2)
xlabel('Position (m)')
ylabel('Population Density (ions/m^3)')
subplot(2,2,4)
plot(xc0,gain0,'k','Linewidth',2)
xlabel('Position (m)')
ylabel('Gain (m^-1)')

figure(2)
subplot(2,2,1)
plot(time0,N20,'k','Linewidth',2)
xlabel('Time (Microseconds)')
ylabel('Gain (m^-1)')
subplot(2,2,2)
plot(time0,Psf0,'Linewidth',2)
xlabel('Time (Microseconds)')
ylabel('Laser Power (W)')
subplot(2,2,3)
plot(time0,N10,'m','Linewidth',2)
xlabel('Time (Microseconds)')
ylabel('N1 (ions/m^3)')
subplot(2,2,4)
plot(time0,Gain0,'r','Linewidth',2)
xlabel('Time (Microseconds)')
ylabel('N2 (ions/m^3)')

figure(3)
subplot(2,2,1)
plot(xc2,ppf2,'b','Linewidth',2)
xlabel('Position (m)')
ylabel('Pump Power (W)')
subplot(2,2,2)
plot(xc2,psf2,'Linewidth',2)
xlabel('Position (m)')
ylabel('Laser Power (W)')
subplot(2,2,3)
plot(xc2,n12,'m','Linewidth',2)
hold on
plot(xc2,n22,'r','Linewidth',2)
xlabel('Position (m)')
ylabel('Population Density (ions/m^3)')
subplot(2,2,4)
plot(xc2,gain2,'k','Linewidth',2)
xlabel('Position (m)')
ylabel('gain (m^-1)')


figure(4)
subplot(2,2,1)
plot(time2,N22,'b','Linewidth',2)
xlabel('Time (Microseconds)')
ylabel('N2 (ions/m^3)')
subplot(2,2,2)
plot(time2,Psf2,'Linewidth',2)
xlabel('Time (Microseconds)')
ylabel('Laser Power (W)')
subplot(2,2,3)
plot(time2,N12,'m','Linewidth',2)
xlabel('Time (Microseconds)')
ylabel('N1 (ions/m^3)')
subplot(2,2,4)
plot(time2,Gain2,'r','Linewidth',2)
xlabel('Time (Microseconds)')
ylabel('Gain (m^-1)')
toc
toc