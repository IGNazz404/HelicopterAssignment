% Filename: EOM_H145D3
% Course: AE4314 Helicopter Performance
% Professor: Marilena Pavel
% Author: Chloë van Droogenbroeck, Ynias Prencipe 
% Student Number: xxxxxxx, 4777158
% Date of Delivery: xx-xx-xxxx

%% 

clc; close all; clear all; 

InducedVelocity_H145D3;
trim_H145D3;

m = 3010;
g = 9.80665;
W = m*g;
Iy = 3163;          % small compared to other file - 10 000, massa is 1000kg lichter
CDS = 1.9;
omega = 38.94;
R = 5.4;
rho = 0.75145;
lok = 6;        % check deze
y_cg = 2.1;         % called mast in other file? 1m


vtip = omega*R;
Arotor = pi*R^2;

tau = 0.1;

% INITIAL CONDITIONS STATES
t0 = 0;
u0 = 0;
w0 = 0;
pitch0 = deg2rad(0);
q0 = 0;
li0 = sqrt(m*g/(Arotor*2*rho))/vtip;

t(1) = t0;
u(1) = u0;
w(1) = w0;
pitch(1) = pitch0;
q(1) = q0;

%li(1) = li0;

dt = 0.01; T = 40; t = [t0:dt:T]; N = length(t);
li = zeros(N, 1);
li(1) = li0;

collective(1) = deg2rad(6);
cyclic(1) = deg2rad(0);

for i= 1:N
   
    % SIMULATE A 1 DEG CYCLIC INPUT BETWEEN 0.5s AND 1s
    if t(i)>= 0.5 & t(i) <= 1
        cyclic(i) = deg2rad(1);
    else
        cyclic(i) = deg2rad(0);
    end
    
    % COLLECTIVE REMAINS CONSTANT
    collective(i) = collective(1); 
    
    % CALCULATE DIMENSIONLESS PITCH RATE AND SPEED
    qdiml(i) = q(i)/omega;
    vdiml(i) = sqrt(u(i)^2+w(i)^2)/vtip;
    
    % CALCULATION OF phi
    if u(i)==0 	
        if w(i)>0 	
            phi(i)=pi/2;
        else phi(i)=-pi/2;
        end
    else
        phi(i)=atan(w(i)/u(i));
    end
    if u(i)<0
        phi(i) = phi(i)+pi;
    end
    alfa_c(i) = cyclic(i)-phi(i);     %(alfa_c)
    
    mu(i) = vdiml(i)*cos(alfa_c(i));  
    lc(i) = vdiml(i)*sin(alfa_c(i)); 

    % CALCULATE a1
    numerator_a1(i) = -16/lok*qdiml(i)+8/3*mu(i)*collective(i)-...
                      2*mu(i)*(lc(i)+li(i));
    a1(i) = numerator_a1(i)/(1-0.5*mu(i)^2);   % van waar komt deze formule?
    
    % TRHUST COEFFICIENTS
        % THRUST COEFFICIENT BEM
        ctelem(i) = Cla*sigma/4*(2/3*collective(i)*(1+1.5*mu(i)^2)-...
               (lc(i)+li(i)));
        % THRUST COEFFICIENT GLAUERT
        alfd(i) = alfa_c(i)-a1(i);  % alfa_c - a1
        ctglau(i) = 2*li(i)*sqrt((vdiml(i)*cos(alfd(i)))^2+(vdiml(i)*...
                    sin(alfd(i))+li(i))^2);
    
    % EQUATIONS OF MOTION
    lidot(i) = ctelem(i);
    
    T(i) = lidot(i)*rho*vtip^2*Arotor;
    cyclic_min_a1(i) = cyclic(i)-a1(i);
    vv(i) = vdiml(i)*vtip;
    
    udot(i) = -g*sin(pitch(i)) - 0.5*rho*CDS/m*u(i)*vv(i)+...
              T(i)/m*sin(cyclic_min_a1(i)) - q(i)*w(i);
    wdot(i) = g*cos(pitch(i)) - CDS/m*.5*rho*w(i)*vv(i)-...
              T(i)/m*cos(cyclic_min_a1(i))+q(i)*u(i); 
    qdot(i) = -T(i)*y_cg/Iy*sin(cyclic_min_a1(i));
    pitchdot(i) = q(i);
    
    % EULER INTEGRATION FOR DISCRETE TIME-SERIES OF EOM
    u(i+1) = u(i) + dt*udot(i);
    w(i+1) = w(i) + dt*wdot(i);
    q(i+1) = q(i) + dt*qdot(i);
    pitch(i+1) =pitch(i) + dt*pitchdot(i);
    
end

u(end) = [];
pitch(end) = [];
w(end) = [];
q(end) = [];

plot(t, rad2deg(cyclic))
xlabel('t (s)')
ylabel('longit grd - cyclic angle')
grid 
pause


plot(t,u, t,rad2deg(pitch)),xlabel('t [s]'),ylabel('u [m/s]'), ylabel('pitch [deg]'),grid,pause;
%plot(t,x),xlabel('t (s)'),ylabel('x(m)'),grid,pause;
plot(t,w),xlabel('t (s)'),ylabel('w(m)'),grid,pause;
plot(t,q),xlabel('t (s)'),ylabel('q(m)'),grid,pause; 
