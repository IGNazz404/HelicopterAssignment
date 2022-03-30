% Filename: trim_H145D3
% Course: AE4314 Helicopter Performance
% Professor: Marilena Pavel
% Author: Chloë van Droogenbroeck, Ynias Prencipe 
% Student Number: xxxxxxx, 4777158
% Date of Delivery: xx-xx-xxxx

%% Helicopter Trim %%

clc; close all; clear all;

InducedVelocity_H145D3;

% Helicopter Parameters
R = 5.4;                %[m]
D = R*2;                %[m]             
c = 0.257;          %[m]
N = 5;       
omega = 38.94;
rho = 0.75145 ;           %[kg/m^3]
W = 3010*9.80665;                %[N]
CDS = 1.9;              % estimated from figure 5.11
Cla = rad2deg(0.1);         % NACA 0015 chosen since symmetrical

sigma = N*c*R/(pi*R^2); %rotor solidity

li_lst = vi_low./(omega*R); 
D_lst = CDS*0.5*rho*V_lst.^2;
mu_lst = V_lst/(omega*R);
CT_lst = 2*li_lst.*sqrt((mu_lst.*cos(D_lst/W)).^2+(mu_lst.*sin(D_lst/W)+li_lst).^2);

trim_values = zeros(2, length(V_lst));

for i = 1:length(V_lst)
    li = li_lst(i);
    D  = D_lst(i);
    mu = mu_lst(i);
    CT = CT_lst(i);

    A = [1+(3/2)*mu^2     -(8/3)*mu;
        -mu            (2/3)+mu^2];
    B = [-2*mu^2*(D/W)-2*mu*li;
         4/sigma*CT/Cla+mu*(D/W)+li];

    trim_values(:,i) = inv(A)*B;
end

a1 = trim_values(1,:);
theta0 = trim_values(2,:);

plot(V_lst, rad2deg(a1), V_lst, rad2deg(theta0))
legend('a1=\theta_c', '\theta_0')

save dumpfile a1 theta0 Cla sigma
clear
load dumpfile