% PIVsuite CC Algorithm Modified

clear all
close all
clc

%% Load Images

clear;

im1 = '../Sample Images/vortex_pair_particles_1.tif';
im2 = '../Sample Images/vortex_pair_particles_2.tif';

  
[pivData1] = pivAnalyzeImagePair(im1,im2);


% Save Velociy Field
ux_cc = pivData1.U;
uy_cc = pivData1.V;

u_mag_cc = (ux_cc.^2+uy_cc.^2).^0.5;
u_mag_cc = (ux_cc.^2+uy_cc.^2).^0.5;
u_max_cc = max(max(u_mag_cc));
u_mag_cc = u_mag_cc/u_max_cc;
figure(1)
vlims=[0, 1];
imagesc(u_mag_cc,vlims);
xlabel('x (pixels)');
ylabel('y (pixels)');
axis image;
set(gca,'YDir','reverse');
colorbar;
colormap jet

% cd('/Users/DavidMoussa/Desktop/CC-OF Hybrid Code/Sample Displacement Field')
% dlmwrite('ux_cc.dat',ux_cc)
% dlmwrite('uy_cc.dat',uy_cc)


