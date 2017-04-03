%% ECE308 Project
%  Spring 2017, modified by Yue
clc; clear all; close all;

%% Design Specifications
tau_s = 0.2;    % settling time
os = 0.20;       % percent overshoot


%% Motor Parameters
K0 = 2;            
K = 40;       
tau = 0.3;       
Kp = 10/pi;        

%% Sampling Time
Ts = 0.001;     % sampling time for discretization
a = 1/tau;

%% Derivation of Motor TF and Root Locus (Uncompensated System)
z = tf('z',Ts);                     % setting up the 'z' variable
G_unc_s = tf(K0*K, [tau 1 0]);      % plant TF in s-domain 
G_unc = c2d(G_unc_s, Ts,'zoh');           % zero order hold
z_p = zero(G_unc);         
p_p = pole(G_unc);
figure(1), rlocus(G_unc);
title('Root Locus Plot of the Uncompensated System');

%% Desired Pole Calculation
zeta = sqrt((log(os)^2)/(pi^2+log(os)^2));   % desired damping ratio
omega_n = 4/(tau_s*zeta);                    % desired natural frequency
p_d_s = -omega_n*zeta +i*omega_n*sqrt(1-zeta^2);                % desired pole in s-domain
p_d_s_2=-omega_n*zeta -i*omega_n*sqrt(1-zeta^2);  
d_1 = exp(Ts*p_d_s);                          % desired pole in z-domain
d_2 = exp(Ts*p_d_s_2);

%% Lead Compensator Design
alpha_c=p_p(2); % Use compensator zero to cancel plant pole
% Multiple choices for zeros, as long as satisfying anagle deficiency

phi = (pi + angle(d_1-alpha_c) + angle(d_1 - z_p)...
    - angle(d_1-p_p(1)) - angle(d_1-p_p(2)))/pi*180; % Angle condition
beta_c = (imag(d_1)/tand(phi)) - real(d_1);
z = tf('z',Ts);
temp = G_unc*(z-alpha_c)/(z+beta_c);
b1=K0*K*(Ts+tau*(exp(-Ts/tau)-1));
b2=K0*K*(tau-(Ts+tau)*(exp(-Ts/tau)));
Kc = norm(p_d_s*(p_d_s+beta_c)*(p_d_s+1/tau)/(K*K0*(p_d_s-alpha_c)/tau)); % Magnitude condition

%% Compensated System (Analysis)

G_cc = (z-alpha_c)/(z+beta_c);
G_c = Kc*G_cc;              % Finish Compensator Design
figure(2), rlocus(series(G_c, G_unc));

figure(3)
rlocus(G_unc); hold on;
plot(d_1, 'r*');plot(d_2,'r*'); hold off;
title('Root Locus of the Uncompensated System');
axis([0.94 1.02 -0.06 0.06]);

figure(4) 
rlocus(series(G_c, G_unc)); hold on;
plot(d_1, 'r*');plot(d_2,'r*'); hold off;
title('Root Locus of the Compensated System');
axis([0.94 1.02 -0.06 0.06]);

figure(5)
stepplot(feedback(series(G_c, G_unc), 1));
axis([0 0.5 0 1.4])
stepinfo(feedback(series(G_c, G_unc), 1))

