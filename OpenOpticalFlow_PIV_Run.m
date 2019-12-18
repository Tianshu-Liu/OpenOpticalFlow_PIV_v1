%% This Matlab program integrats the optical flow method with the cross-correlation method
%% for extraction of high-resolution velocity fields from particle images.  
%% This hybrid method provides an additional tool to process PIV images, 
%% which combines the advantages of the optical flow method and cross-correlation method.  

close all
clear all
clc

%% Load Images

% Im1=imread('I1_vortexpair_dt0p03.tif');
% Im2=imread('I2_vortexpair_dt0p03.tif');

% Im1=imread('I1_vortexpair_dt0p05.tif');
% Im2=imread('I2_vortexpair_dt0p05.tif');

% Im1=imread('wall_jet_1.tif');
% Im2=imread('wall_jet_2.tif');

% Im1=imread('twin_jet_1.tif');
% Im2=imread('twin_jet_2.tif');

Im1=imread('White_Oval_1.tif');
Im2=imread('White_Oval_2.tif');



%% Selete region of interest, "0" for processing the whole image, "1" for processing a selected region
index_region=0; 

Im1=double(Im1);
Im2=double(Im2);

if (index_region == 1)
    imagesc(uint8(Im1));
    colormap(gray);
    axis image;

    xy=ginput(2);
    x1=floor(min(xy(:,1)));
    x2=floor(max(xy(:,1)));
    y1=floor(min(xy(:,2)));
    y2=floor(max(xy(:,2)));
    I1=double(Im1(y1:y2,x1:x2)); 
    I2=double(Im2(y1:y2,x1:x2));
elseif (index_region == 0)
    I1=Im1;
    I2=Im2;
end

I1_original=I1;
I2_original=I2;


%% Set the Parameters for Optical Flow Computation

% Set the lagrange multipleirs in optical computation 
lambda_1=20;  % the Horn_schunck estimator for initial field
lambda_2=2000; % the Liu-Shen estimator for refined estimation

%% Number of iterations in the coarse-to-fine iterative process from
%% initial estimation, "1" means 1 iteration
no_iteration=1; % fixed

%% Initial coarse field estimation in the coarse-to-fine iterative process,
%% scale_im is a scale factor for down-sizing of images
scale_im=1; 
%% For Image Pre-Processing

%% For local illumination intensity adjustment, To bypass it, set size_average = 0
size_average=0; % in pixels

%% Gausian filter size for removing random noise in images
size_filter=6; % in pixels

%% correcting the global and local intensity change in images
[m1,n1]=size(I1);
window_shifting=[1;n1;1;m1]; % [x1,x2,y1,y2] deines a rectangular window for global intensity correction
[I1,I2]=correction_illumination(I1,I2,window_shifting,size_average);


%% cleaning the left and upper edges since some data near the edges are corrupted due to interperlation
edge_width=1; % in pixels


%% pre-processing for reducing random noise,
%% and downsampling images if displacements are large
[I1,I2] = pre_processing_a(I1,I2,scale_im,size_filter);

I_region1=I1;
I_region2=I2;


%% initial correlation calculation for a coarse-grained velocity field (ux0,uy0)
% ux is the velocity (pixels/unit time) in the image x-coordinate (from the left-up corner to right)
% uy is the velocity (pixels/unit time) in the image y-coordinate (from the left-up corner to bottom)


%% run FFT cross-correlation algorithm
Im1=I1;
Im2=I2;

pivPar.iaSizeX = [64 16 8];     % size of interrogation area in X 
pivPar.iaStepX = [32 8 4];       % grid spacing of velocity vectors in X 

pivPar.ccMethod = 'fft'; 

[pivData1] = pivAnalyzeImagePair(Im1,Im2,pivPar);

ux0=pivData1.U;
uy0=pivData1.V;


%% re-size the initial velocity dield (u0, v0)
[n0,m0]=size(ux0);
[n1,m1]=size(Im1);

scale=round((n1*m1/(n0*m0))^0.5);

ux0=imresize(ux0,scale);
uy0=imresize(uy0,scale);


%% generate the shifted image from Im1 based on the initial coarse-grained velocity field (ux0, uy0),
%% and then calculate velocity difference for iterative correction


%% estimate the displacement vector and make correction in iterations

ux=ux0;
uy=uy0;

k=1;
while k<=no_iteration
   [Im1_shift,uxI,uyI]=shift_image_fun_refine_1(ux,uy,Im1,Im2);
    
    I1=double(Im1_shift);
    I2=double(Im2);
               
    % calculation of correction of the optical flow 
    [dux,duy,vor,dux_horn,duy_horn,error2]=OpticalFlowPhysics_fun(I1,I2,lambda_1,lambda_2);

    % refined optical flow
    ux_corr=uxI+dux;
    uy_corr=uyI+duy;
    
    
    k=k+1;
end

%% refined velocity field
ux = ux_corr;    %%%%%
uy = uy_corr;    %%%%%



%% clean up the edges
ux(:,1:edge_width)=ux(:,(edge_width+1):(2*edge_width));
uy(:,1:edge_width)=uy(:,(edge_width+1):(2*edge_width));

ux(1:edge_width,:)=ux((edge_width+1):(2*edge_width),:);
uy(1:edge_width,:)=uy((edge_width+1):(2*edge_width),:);



%% show the images and processed results
%% plot the images, velocity vector, and streamlines in the initail and
%% refined estimations
 plots_set_1;

%% plot the fields of velocity magnitude, vorticity and the second invariant Q
 plots_set_2;
 
 
% save Ux_vortexpair_hybrid_dt0p2.dat ux -ascii;
% save Uy_vortexpair_hybrid_dt0p2.dat uy -ascii;
%  
% save Ux_vortexpair_corr_dt0p2.dat ux0 -ascii;
% save Uy_vortexpair_corr_dt0p2.dat uy0 -ascii;
% %  
 
 
 
 