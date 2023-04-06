close all;
clear all;
clc;

%% Read Plasma Profiles

fileName='psep1p5mw.h5';
% Grid info
% ---------
r=h5read(fileName,'/solps_like/r');
z=h5read(fileName,'/solps_like/z');

% B-field
% ---------------
br=h5read(fileName,'/bfield/b_r')';
bphi=h5read(fileName,'/bfield/b_phi')';
bz=h5read(fileName,'/bfield/b_z')';

% Electron profiles
% ----------------
ne=h5read(fileName,'/n_e/dens')';
vez=h5read(fileName,'/n_e/parr_flow')';
te=h5read(fileName,'/n_e/temp')';

% Dueterium profiles
% ------------------
ni=h5read(fileName,'/n_i/dens')';
viz=h5read(fileName,'/n_i/parr_flow')';
ti=h5read(fileName,'/n_i/temp')';

% Oxygen profiles
% ---------------

numOxygenSpecies=8;

for ss=1:numOxygenSpecies
% O+1 profiles

speciesDensity = ['/o_',int2str(ss),'/dens'];
speciesParrFlow = ['/o_',int2str(ss),'/parr_flow'];
speciesTemp = ['/o_',int2str(ss),'/temp'];


n_o{ss} = h5read(fileName,speciesDensity)';
t_o{ss} = h5read(fileName,speciesParrFlow)';
vz_o{ss} = h5read(fileName,speciesTemp)';

end

r=reshape(r,4,165*146);
z=reshape(z,4,165*146);




zW=h5read(fileName,'/z_wall_points');
rW=h5read(fileName,'/r_wall_points');

figure(1); patch(r,z,n_o{8}(:),'EdgeColor','k');
hold on;
 
plot(rW,zW,'r-', 'LineWidth', 2);
set(gca, 'ColorScale', 'log', 'FontSize',18)
xlabel('R[m]', 'Interpreter','latex')
ylabel('z[m]', 'Interpreter','latex')
hold off;




xlim([1.5 3.5])
ylim([-1 1])
axis equal
colormap('jet')
colorbar
% 
% 
%% SOLEDGE2D to GITR grid

%% Plasma (electron + D ions) profile

% Querry points
% -------------
rq = 2.2;
zq = -0.5;
rgrid = linspace(1.8,3.3,2000);
zgrid = linspace(-1,1,4000);

[r_mesh z_mesh] = meshgrid(rgrid,zgrid);

% B-field profile

br_q = 0*r_mesh;
bz_q = 0*r_mesh;
bphi_q = 0*r_mesh;

    for i=1:length(r)
        [in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
        if (length(find(in)) > 0 || length(find(on)) > 0)
        
        % figure; patch(r(:,i),z(:,i),data(i),'EdgeColor','k');
        % hold on;
        % scatter(rq,zq,'g')
        br_q(find(in)) = br(i);
        i
        end
    end
        
    for i=1:length(r)
        [in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
        if (length(find(in)) > 0 || length(find(on)) > 0)
        
        bz_q(find(in)) = bz(i);
        i
        end
    end
        
    for i=1:length(r)
        [in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
        if (length(find(in)) > 0 || length(find(on)) > 0)
        
        bphi_q(find(in)) = bphi(i);
        i
        end
    end

figure
h = imagesc(rgrid,zgrid,br_q);
set(gca,'YDir','normal')
% h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','linear')
xlim([1.8 3.2])
ylim

figure
h = imagesc(rgrid,zgrid,bz_q);
set(gca,'YDir','normal')
% h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','linear')
xlim([1.8 3.2])
ylim

figure
h = imagesc(rgrid,zgrid,bphi_q);
set(gca,'YDir','normal')
% h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','linear')
xlim([1.8 3.2])
ylim

% Electron profile
% ----------------
ne_q = 0*r_mesh;
te_q = 0*r_mesh;
vez_q = 0*r_mesh;

for i=1:length(r)
[in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
if (length(find(in)) > 0 || length(find(on)) > 0)

% figure; patch(r(:,i),z(:,i),data(i),'EdgeColor','k');
% hold on;
% scatter(rq,zq,'g')
ne_q(find(in)) = ne(i);
i
end
end

for i=1:length(r)
[in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
if (length(find(in)) > 0 || length(find(on)) > 0)

te_q(find(in)) = te(i);
i
end
end

for i=1:length(r)
[in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
if (length(find(in)) > 0 || length(find(on)) > 0)

vez_q(find(in)) = vez(i);
i
end
end

figure
h = imagesc(rgrid,zgrid,ne_q);
set(gca,'YDir','normal')
% h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','log')
xlim([1.8 3.2])
ylim

figure
h = imagesc(rgrid,zgrid,te_q);
set(gca,'YDir','normal')
% h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','log')
xlim([1.8 3.2])
ylim

figure
h = imagesc(rgrid,zgrid,vez_q);
set(gca,'YDir','normal')
% h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','linear')
xlim([1.8 3.2])
ylim

%% D+ Ion profiles

% D+ Ion density
% ----------------
[r_mesh z_mesh] = meshgrid(rgrid,zgrid);
ni_q = 0*r_mesh;
ti_q = 0*r_mesh;
viz_q = 0*r_mesh;

    for i=1:length(r)
        [in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
        if (length(find(in)) > 0 || length(find(on)) > 0)
        
        % figure; patch(r(:,i),z(:,i),data(i),'EdgeColor','k');
        % hold on;
        % scatter(rq,zq,'g')
        ni_q(find(in)) = ni(i);
        i
        end
    end
        
    for i=1:length(r)
        [in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
        if (length(find(in)) > 0 || length(find(on)) > 0)
        
        ti_q(find(in)) = ti(i);
        i
        end
    end
        
    for i=1:length(r)
        [in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
        if (length(find(in)) > 0 || length(find(on)) > 0)
        
        viz_q(find(in)) = viz(i);
        i
        end
    end

figure
h = imagesc(rgrid,zgrid,ni_q);
set(gca,'YDir','normal')
% h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','log')
xlim([1.8 3.2])
ylim


figure
h = imagesc(rgrid,zgrid,ti_q);
set(gca,'YDir','normal')
% h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','log')
xlim([1.8 3.2])
ylim

figure
h = imagesc(rgrid,zgrid,viz_q);
set(gca,'YDir','normal')
% h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','linear')
xlim([1.8 3.2])
ylim


%% Oxygen profile

% O8+ density
% ----------------

numOxygenSpecies=2;

for ss=1:numOxygenSpecies
% O+1 profiles

speciesDensity = ['/o_',int2str(ss),'/dens'];
speciesParrFlow = ['/o_',int2str(ss),'/parr_flow'];
speciesTemp = ['/o_',int2str(ss),'/temp'];


n_o{ss} = h5read(fileName,speciesDensity)';
t_o{ss} = h5read(fileName,speciesParrFlow)';
vz_o{ss} = h5read(fileName,speciesTemp)';

no_q{ss} = 0*r_mesh;
to_q{ss} = 0*r_mesh;
vzo_q{ss} = 0*r_mesh;
    for i=1:length(r)
        [in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
        if (length(find(in)) > 0 || length(find(on)) > 0)
    
        no_q{ss}(find(in)) = n_o{ss}(i);
        i
        end
    end

    for i=1:length(r)
        [in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
        if (length(find(in)) > 0 || length(find(on)) > 0)
        
        to_q{ss}(find(in)) = t_o{ss}(i);
        i
        end
    end


    for i=1:length(r)
        [in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
        if (length(find(in)) > 0 || length(find(on)) > 0)
        
        vzo_q{ss}(find(in)) = vz_o{ss}(i);
        i
        end
    end

end

figure
h = imagesc(rgrid,zgrid,no_q{1});
set(gca,'YDir','normal')
% h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','log')
xlim([1.8 3.2])
ylim

figure
h = imagesc(rgrid,zgrid,to_q{1});
set(gca,'YDir','normal')
% h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','log')
xlim([1.8 3.2])
ylim

figure
h = imagesc(rgrid,zgrid,vzo_q{1});
set(gca,'YDir','normal')
% h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','linear')
xlim([1.8 3.2])
ylim


disp('Saving .mat file for this interpolation')
fileName = fileName(1:strfind(fileName,'.')-1);
save(fileName)

return 


%% Values at the wall
val_wall = interpn(zgrid,rgrid,ne_q,zW,rW)

soledge_wall = 0*rW;
for i=1:length(r)
[in,on] = inpolygon(rW,zW,r(:,i),z(:,i));
if (length(find(in)) > 0 || length(find(on)) > 0)


% figure; patch(r(:,i),z(:,i),data(i),'EdgeColor','k');
% hold on;
% scatter(rq,zq,'g')
soledge_wall(find(in)) = data(i);
i
end
end

figure
plot(val_wall)
hold on
plot(soledge_wall)