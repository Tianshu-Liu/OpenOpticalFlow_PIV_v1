clear all;
close all;


Ux_exact_1=load('Ux_vortexpair_true_dt0p01.dat');
Uy_exact_1=load('Uy_vortexpair_true_dt0p01.dat');

Ux_exact_2=load('Ux_vortexpair_true_dt0p02.dat');
Uy_exact_2=load('Uy_vortexpair_true_dt0p02.dat');

Ux_exact_3=load('Ux_vortexpair_true_dt0p03.dat');
Uy_exact_3=load('Uy_vortexpair_true_dt0p03.dat');

Ux_exact_4=load('Ux_vortexpair_true_dt0p05.dat');
Uy_exact_4=load('Uy_vortexpair_true_dt0p05.dat');

Ux_exact_5=load('Ux_vortexpair_true_dt0p1.dat');
Uy_exact_5=load('Uy_vortexpair_true_dt0p1.dat');

Ux_exact_6=load('Ux_vortexpair_true_dt0p2.dat');
Uy_exact_6=load('Uy_vortexpair_true_dt0p2.dat');


Ux_liu_1=load('Ux_vortexpair_hybrid_dt0p01.dat');
Uy_liu_1=load('Uy_vortexpair_hybrid_dt0p01.dat');

Ux_liu_2=load('Ux_vortexpair_hybrid_dt0p02.dat');
Uy_liu_2=load('Uy_vortexpair_hybrid_dt0p02.dat');

Ux_liu_3=load('Ux_vortexpair_hybrid_dt0p03.dat');
Uy_liu_3=load('Uy_vortexpair_hybrid_dt0p03.dat');

Ux_liu_4=load('Ux_vortexpair_hybrid_dt0p05.dat');
Uy_liu_4=load('Uy_vortexpair_hybrid_dt0p05.dat');

Ux_liu_5=load('Ux_vortexpair_hybrid_dt0p1.dat');
Uy_liu_5=load('Uy_vortexpair_hybrid_dt0p1.dat');

Ux_liu_6=load('Ux_vortexpair_hybrid_dt0p2.dat');
Uy_liu_6=load('Uy_vortexpair_hybrid_dt0p2.dat');


Ux_corr_1=load('Ux_vortexpair_corr_dt0p01.dat');
Uy_corr_1=load('Uy_vortexpair_corr_dt0p01.dat');

Ux_corr_2=load('Ux_vortexpair_corr_dt0p02.dat');
Uy_corr_2=load('Uy_vortexpair_corr_dt0p02.dat');

Ux_corr_3=load('Ux_vortexpair_corr_dt0p03.dat');
Uy_corr_3=load('Uy_vortexpair_corr_dt0p03.dat');

Ux_corr_4=load('Ux_vortexpair_corr_dt0p05.dat');
Uy_corr_4=load('Uy_vortexpair_corr_dt0p05.dat');

Ux_corr_5=load('Ux_vortexpair_corr_dt0p1.dat');
Uy_corr_5=load('Uy_vortexpair_corr_dt0p1.dat');

Ux_corr_6=load('Ux_vortexpair_corr_dt0p2.dat');
Uy_corr_6=load('Uy_vortexpair_corr_dt0p2.dat');


% 

[m,n]=size(Ux_exact_2);

Ux_corr_a_1 = imresize(Ux_corr_1, [m n]); 
Uy_corr_a_1 = imresize(Uy_corr_1, [m n]); 

Ux_corr_a_2 = imresize(Ux_corr_2, [m n]); 
Uy_corr_a_2 = imresize(Uy_corr_2, [m n]); 

Ux_corr_a_3 = imresize(Ux_corr_3, [m n]); 
Uy_corr_a_3 = imresize(Uy_corr_3, [m n]); 

Ux_corr_a_4 = imresize(Ux_corr_4, [m n]); 
Uy_corr_a_4 = imresize(Uy_corr_4, [m n]); 

Ux_corr_a_5 = imresize(Ux_corr_5, [m n]); 
Uy_corr_a_5 = imresize(Uy_corr_5, [m n]); 

Ux_corr_a_6 = imresize(Ux_corr_6, [m n]); 
Uy_corr_a_6 = imresize(Uy_corr_6, [m n]); 



%dt=0.02;
% U_mag_exact=(Ux_exact_2.^2+Uy_exact_2.^2).^0.5;
% U_mag_liu=(Ux_liu_2.^2+Uy_liu_2.^2).^0.5;
% U_mag_corr_a=(Ux_corr_1a.^2+Uy_corr_1a.^2).^0.5;
% U_mag_corr_b=(Ux_corr_1b.^2+Uy_corr_1b.^2).^0.5;


%x0=floor(n/2);
x0=250;
y=[1:1:m];


figure(2);
plot(y,Ux_exact_4(1:m,x0),'.k',y,Ux_liu_4(:,x0),'-k',y(1:8:end),Ux_corr_a_4(1:8:end,x0),'+k');
grid;
axis([0 500 -3 5]);
xlabel('y (pixels)');
ylabel('u_x (pixels/unit time)');
legend('Truth', 'Hybrid Method', 'Correlation Method 1');


% plot velocity vector field
figure(3);
gx=30; offset=1;
h = vis_flow (Ux_liu_4, Uy_liu_4, gx, offset, 3, 'm');
set(h, 'Color', 'black');
xlabel('x (pixels)');
ylabel('y (pixels)');
axis image;
set(gca,'YDir','reverse');
title('Velocity Field (Hybrid Method)');


% plot streamlines
figure(4);
[m,n]=size(Ux_liu_6);
[x,y]=meshgrid(1:n,1:m);
dn=10;
dm=10;
[sx,sy]=meshgrid(1:dn:n,1:dm:m);
%h=streamline(x, y, ux, uy, sx, sy);
h=streamslice(x, y, Ux_liu_4, Uy_liu_4, 4);
set(h, 'Color', 'black');
xlabel('x (pixels)');
ylabel('y (pixels)');
axis image;
set(gca,'YDir','reverse');
title('Streamlines (Hybrid Method)');


% plot velocity vector field
figure(5);
gx=30; offset=1;
h = vis_flow (Ux_corr_a_4(1:8:end,1:8:end), Uy_corr_a_4(1:8:end,1:8:end), gx, offset, 3, 'm');
set(h, 'Color', 'black');
xlabel('x (regions)');
ylabel('y (regions)');
axis image;
set(gca,'YDir','reverse');
title('Velocity Field (Correlation Method 1)');


% plot streamlines
figure(6);
[m,n]=size(Ux_corr_a_6);
[x,y]=meshgrid(1:n,1:m);
dn=10;
dm=10;
[sx,sy]=meshgrid(1:dn:n,1:dm:m);
%h=streamline(x, y, ux, uy, sx, sy);
h=streamslice(x(1:8:end,1:8:end), y(1:8:end,1:8:end), Ux_corr_a_4(1:8:end,1:8:end), Uy_corr_a_4(1:8:end,1:8:end), 4);
set(h, 'Color', 'black');
xlabel('x (regions)');
ylabel('y (regions)');
axis image;
set(gca,'YDir','reverse');
title('Streamlines (Correlation Method 1)');




for i=2:6
    Ux_op_name=strcat('Ux_liu_',num2str(i));
    Uy_op_name=strcat('Uy_liu_',num2str(i));
    
    Ux_liu=eval(Ux_op_name);
    Uy_liu=eval(Uy_op_name);
    
    Ux_corr_name=strcat('Ux_corr_a_',num2str(i));
    Uy_corr_name=strcat('Uy_corr_a_',num2str(i));
    
    Ux_corr=eval(Ux_corr_name);
    Uy_corr=eval(Uy_corr_name);
    
    Ux_exact_name=strcat('Ux_exact_',num2str(i));
    Uy_exact_name=strcat('Uy_exact_',num2str(i));
    
    Ux_exact=eval(Ux_exact_name);
    Uy_exact=eval(Uy_exact_name);
    
    dU_liu=((Ux_liu(50:end,50:end)-Ux_exact(50:end,50:end)).^2+(Uy_liu(50:end,50:end)-Uy_exact(50:end,50:end)).^2).^0.5;
    dU_corr=((Ux_corr(1:8:end,1:8:end)-Ux_exact(1:8:end,1:8:end)).^2+(Uy_corr(1:8:end,1:8:end)-Uy_exact(1:8:end,1:8:end)).^2).^0.5;
    U_mag=(Ux_exact.^2+Uy_exact.^2).^0.5;   
    U_max(i)=max(max(U_mag));
    
    grad_mag=gradient(Ux_exact(50:n,50:n),Uy_exact(50:n,50:n));   
    grad_avg(i)=mean(mean(grad_mag));
        
    error_liu(i)=mean(mean(dU_liu));
    error_corr(i)=mean(mean(dU_corr));
end


% data_out=[U_max' error_liu'];
% save error_liu_(50_5000).dat data_out -ascii;



data_err=load('errors_OF_lavision.dat');

figure(20);
plot(U_max,error_liu,'o-k',data_err(:,1), data_err(:,2),'-sk',U_max,error_corr,'>-k',data_err(:,1), data_err(:,3),'-dk');
grid;
xlabel('Max Displacement (pixels)');
ylabel('RMS Error (pixels)');
legend('Hybrid Method','Optical Flow Method','Correlation Method 1','Correlation Method 2');


figure(21);
plot(U_max,error_liu./U_max,'o-',U_max,error_corr./U_max,'>-k');
grid;
xlabel('Max Displacement (pixels)');
ylabel('Relative RMS Error (pixels)');
legend('Hybrid Method','Correlation Method');


figure(22);
mesh(dU_liu);
axis([1 500 1 500 0 1.5]);
xlabel('x (pixels)');
ylabel('y (pixels)');
zlabel('Error (pixels/unit time)');
title('Hybrid Method');


figure(23);
mesh(dU_corr);
axis([1 62 1 62 0 1.5]);
xlabel('x (regions)');
ylabel('y (regions)');
zlabel('Error (pixels/unit time)');
title('Correlation Method');    
    
    
    
    
figure(30);
plot(grad_avg,error_liu,'o-k',grad_avg, data_err(:,2),'-sk',grad_avg,error_corr,'>-k',...
    grad_avg, data_err(:,3),'-dk');
grid;
xlabel('Velocity Gradient Magnitude (1/unit time)');
ylabel('RMS Error (pixels)');
legend('Hybrid Method','Optical Flow Method','Correlation Method 1','Correlation Method 2');   
    
    
    
    
    















