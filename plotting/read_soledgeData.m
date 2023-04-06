close all;
clear all;
clc;

r=h5read('3mw_test_3.h5','/solps_like/r');
z=h5read('3mw_test_3.h5','/solps_like/z');

% rCentroid=sum(r,1)./4;
% zCentroid=sum(z,1)./4;



data=h5read('3mw_test_3.h5','/n_e/dens');
data=data';
r=reshape(r,4,165*146);
z=reshape(z,4,165*146);



zW=h5read('3mw_test_3.h5','/z_wall_points');
rW=h5read('3mw_test_3.h5','/r_wall_points');

figure(1); patch(r,z,data(:),'EdgeColor','k');
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

%%
rq = 2.2;
zq = -0.5;
rgrid = linspace(1.8,3.3,100);
zgrid = linspace(-1,1,150);

[r_mesh z_mesh] = meshgrid(rgrid,zgrid);
data_q = 0*r_mesh;
for i=1:length(r)
[in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
if (length(find(in)) > 0 || length(find(on)) > 0)


% figure; patch(r(:,i),z(:,i),data(i),'EdgeColor','k');
% hold on;
% scatter(rq,zq,'g')
data_q(find(in)) = data(i);
i
end
end

figure
h = pcolor(rgrid,zgrid,data_q);
h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','log')

val_wall = interpn(rgrid,zgrid,data_q,rW,zW)

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