close all
clear all


file = strcat('output/positions.nc');
hitWall = ncread(file,'hitWall');
nHit = length(find(hitWall));
hasHit = find(hitWall);
notHit = find(hitWall==0);
x0 = ncread(file,'x');
y0 = ncread(file,'y');
z0 = ncread(file,'z');
vx0 = ncread(file,'vx');
vy0 = ncread(file,'vy');
vz0 = ncread(file,'vz');
time0=ncread(file,'time');
distTraveled = ncread(file,'distTraveled');
charge0 = ncread(file,'charge');
weight0 = ncread(file,'weight');
vtot = sqrt(vx0.^2 +vy0.^2 + vz0.^2);
E = 0.5*27*1.66e-27*vtot.^2/1.602e-19;
%E = 0.5*16*1.66e-27*vtot.^2/1.602e-19;
% figure(11)
% histogram(E)
% 
% figure()
% scatter3(x0(hasHit),y0(hasHit),z0(hasHit))
% hold on
% scatter3(x0(notHit),y0(notHit),z0(notHit))

r0 = sqrt(x0.^2 + y0.^2);

% figure
% scatter(r0,z0)

specFile = strcat('output/spec.nc');
Chargedens = ncread(specFile,'n');
gridR = ncread(specFile,'gridR');
gridZ = ncread(specFile,'gridZ');

% figure
slice1 = Chargedens(:,:,end)
% 
figure
h1 = pcolor(gridR,gridZ,log10(slice1'))
h1.EdgeColor = 'none';

% profilesFile = strcat('input/profiles.nc');
% br = ncread(profilesFile,'br');
% bt = ncread(profilesFile,'bt');
% bz = ncread(profilesFile,'bz');
% gridR = ncread(profilesFile,'gridR');
% gridZ = ncread(profilesFile,'gridZ');
% % br = reshape(br,1,[]);
% % br = reshape(br,length(gridZ),length(gridR));
% figure
% h1 = pcolor(gridR,gridZ,br')
% h1.EdgeColor = 'none';
% colorbar
% 
% figure
% h1 = pcolor(gridR,gridZ,bt')
% h1.EdgeColor = 'none';
% colorbar
% 
% figure
% h1 = pcolor(gridR,gridZ,bz')
% h1.EdgeColor = 'none';
% colorbar