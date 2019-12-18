function [grad_mag]=gradient(Vx, Vy)

% Vx = imfilter(Vx, [1 1 1 1 1]'*[1 1 1 1 1]/25,'symmetric');
% Vy = imfilter(Vy, [1 1 1 1 1]'*[1,1 1 1,1]/25,'symmetric');

dx=1;
D = [0, -1, 0; 0,0,0; 0,1,0]/2; %%% partial derivative 
Vy_x = imfilter(Vy, D'/dx, 'symmetric',  'same'); 
Vx_y = imfilter(Vx, D/dx, 'symmetric',  'same');
grad_mag=(Vy_x.^2+Vx_y.^2).^0.5;




















