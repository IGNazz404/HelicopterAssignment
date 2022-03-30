%% Helicopter Trim %%

% Helicopter Parameters
R = 5.4;                %[m]
D = R*2;                %[m]             
c = 0.257;          %[m]
N = 5;       
omega = 39.4794;
rho = 0.75145 ;           %[kg/m^3]
W = 3010*9.80665;                %[N]
CDS = 1.9;              % estimated from figure 5.11
V_max = 80;     % max speed in [m/s]
Cla = 0.1;         % NACA 0015 chosen since symmetrical

sigma = N*c/(pi*R); %rotor solidity

V_lst = [0:0.1:V_max];
x_lst = zeros(2,length(V_lst));

for i = 1:length(V_lst)
V = V_lst(i);
D = CDS*rho*0.5*V^2;
T = sqrt(W^2+D^2);
CT = T/(rho*(omega*R)^2*pi*R^2);
mu = V/(omega*R);

coeff = [4  8*mu*sin(D/W)   4*mu^2 0 -CT^2];
r = roots(coeff);
r_real =r(r>=0);
lambda_i = r_real;

A = [1+3/2*mu^2     -8/3*mu;
     -mu            2/3+mu^2];

B = [-2*mu^2 * (D/W) - 2*mu*sol;
     4/sigma * CT/Cla + mu*(D/W) + sol];

x_lst(:,i) = (inv(A)*B);
end

a1 = x_lst(1,:);
theta0 = x_lst(2,:);

plot(V_lst, a1, V_lst, theta0)