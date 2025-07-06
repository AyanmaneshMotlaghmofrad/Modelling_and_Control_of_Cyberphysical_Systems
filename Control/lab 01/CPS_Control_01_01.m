clear variables
close all
clc
format compact

%%
s = tf('s');
z = tf('z',-1);

Gp = 100/(s^2 + 1.2*s + 1)

Gp_d = c2d(Gp,1,'zoh')

[num_d, den_d] = tfdata(Gp_d,'v')
theta = [den_d num_d]'
%% simulation
H = 100000;% in the noise free case 7 samplesa are enough

u = rand([H 1]);

w = lsim(Gp_d,u);

%% Shaping the matrix A and the vector b

% w(k) = den_d(2)*w(k-1) + den_d(3)*w(k-2) - num_d(2)*u(k-1) - num_d(3)*u(k-2) = 0 

% [w(3:H) ] = [-w(2:H-1) -w(1:H-2) u(3:H) u(2:H-1) u(1:H-2)]*[den_d(2) den_d(3) num_d(1) num_d(2) num_d(3)]'

b = w(3:H,1);

A = [-w(2:H-1) -w(1:H-2) u(3:H) u(2:H-1) u(1:H-2)];

theta_eq_er = A \ b

num_hat = theta_eq_er(3:end)'
den_hat = [1 theta_eq_er(1:2)']

[num_hat' num_d']

[den_hat' den_d']

%% equation error

eq_error = 5*randn([H 1]);

% corrupting the equations
y_eq = zeros(H,1);

for i = 3:H
    y_eq(i) = -theta(2)*y_eq(i-1) - theta(3)*y_eq(i-2) + theta(5)*u(i-1) + theta(6)*u(i-2) + eq_error(i);
end

A_eq = [-y_eq(2:H-1) -y_eq(1:H-2) u(3:H) u(2:H-1) u(1:H-2)];
b_eq_error = y_eq(3:end);

theta_eq_er = A_eq \ b_eq_error;

num_hat_eq = theta_eq_er(3:end)';

den_hat_eq = [1 theta_eq_er(1:2)'];

[num_hat_eq' num_d']

[den_hat_eq' den_d']

% it can be seen that LS is effective with 
% equation error, since consistency properties
% assumptions are satisfied.

%% output error

eta = 5*randn([H 1]); % gaussian noise assumption

y = w + eta;

b_OE = y(3:H);
A_OE = [-y(2:H-1) -y(1:H-2) u(3:H) u(2:H-1) u(1:H-2)];

theta_OE = A_OE \ b_OE

num_hat_OE = theta_OE(3:end)'

den_hat_OE = [1 theta_OE(1:2)']

[num_hat_OE' num_d']

[den_hat_OE' den_d']

% also for large values of H, the value of theta
% does not converge to the real parameters.