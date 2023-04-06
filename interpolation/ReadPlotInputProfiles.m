close all;
clear all;
clc;

% Physical constants:
% 

%%

r=ncread('profiles.nc','gridR');
z=ncread('profiles.nc','gridZ');

ne=ncread('profiles.nc','ne');
figure; imagesc(r,z,ne);
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
xlabel('$r$ [m]','interpreter','Latex','fontSize',18);
ylabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input Density')
colorbar;



te=ncread('profiles.nc','te');
figure; imagesc(r,z,te);
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
xlabel('$r$ [m]','interpreter','Latex','fontSize',18);
ylabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input Te')
colorbar;

% vx=ncread('../TeX1/input/profilesProtoMPEX.nc','vx');
% vx=ncread('../TeX1/input/profilesProtoMPEX.nc','vy');
vz=ncread('profiles.nc','ve');
figure;imagesc(r,z,vz);
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
xlabel('$r$ [m]','interpreter','Latex','fontSize',18);
ylabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input Vz')
colorbar;
% figure; plot(vz(1,:))
ne1D=mean(ne,2);

figure, plot(r,ne1D)



