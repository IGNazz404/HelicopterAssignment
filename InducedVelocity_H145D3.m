% Filename: InducedVelocity_H145D3
% Course: AE4314 Helicopter Performance
% Professor: Marilena Pavel
% Author: Chloë van Droogenbroeck, Ynias Prencipe 
% Student Number: xxxxxxx, 4777158
% Date of Delivery: xx-xx-xxxx

%% Induced Velocity Calculations

clc; close all; clear all;

m_hover = 3010;      
g = 9.80665;
D = 10.8    ;        
R = D/2     ;        
rho = 0.75145;       
Vmax = 80;

% Induced Velocity Hover
T = m_hover * g;
vi_hov = sqrt(T/(2*rho*pi*R^2));


% Induced Velocity Forward Flight
V_lst = [0:1:Vmax];
vi_lst = sqrt(-V_lst.^2/2 + sqrt(V_lst.^4/4+1));
vi_low = zeros(1, length(V_lst));
vi_high = zeros(1, length(V_lst));

for i = 1:length(V_lst)
    V = V_lst(i);
    
    vi_fwd_lowspeed = sqrt(-V^2/2 + sqrt((V^2/2)^2 + vi_hov^4));
    vi_fwd_hispeed = vi_hov^2/V;
    
    vi_low(:,i) =  vi_fwd_lowspeed;
    vi_high(:,i) =  vi_fwd_hispeed;
end


plot(V_lst, vi_low, V_lst, vi_high, V_lst, vi_lst)
axis([0 80 0 30]);
ylabel('Induced Velocity v_i [m/s]'); xlabel('Forward Speed V [m/s]');
legend("Induced Velocity Low Speed", "Induced Velocity High speed", "vi_{bar}");

save dumpfile vi_low Vmax V_lst vi_lst
clear
load dumpfile
