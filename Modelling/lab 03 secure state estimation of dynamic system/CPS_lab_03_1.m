%secure state estimation of a dynamical system
% in the case of A = cte
clc
clear variables
close all
format compact

%%
load dynamic_CPS_data.mat

[n,~] = size(A);
q = length(a);
G = [C eye(q)];


%% part 1
max_iteration = 1e5;

% SSO
nu = 0.99/(norm(G,2)^2);
lambda = 0.1;

% initialization
estimation_error_SSO = zeros(1,max_iteration);
attack_support_error_SSO = zeros(1,max_iteration);

x = [x0,zeros(n,max_iteration-1)];
x_hat_SSO = zeros(n,max_iteration);

a_hat_SSO = zeros(q,max_iteration);

for k = 1:max_iteration
    % system simulation and 
    y= C*x(:,k) + a;
    y_hat_SSO = C*x_hat_SSO(:,k) + a_hat_SSO(:,k);
    x_hat_SSO(:,k+1) = A*x_hat_SSO(:,k) - nu*A*C'*(y_hat_SSO - y);
    a_hat_SSO(:,k+1) = prox_l1(a_hat_SSO(:,k) - nu*(y_hat_SSO - y),nu*lambda);
    x(:,k+1) = A*x(:,k);
    
    % errors
    attack_support_error_SSO(k) = sum((a_hat_SSO(:,k+1)~=0) ~= (a~= 0));
    estimation_error_SSO(k) = norm(x_hat_SSO(:,k+1)-x(:,k+1))/norm(x(:,k+1));
    
end

% D-SSO
nu = 0.7;
lambda = 0.1;

desired_lambda = 0.01*abs(rand(1,n));

L = place(A',C',desired_lambda);

% initialization
estimation_error_DSSO = zeros(1,max_iteration);
attack_support_error_DSSO = zeros(1,max_iteration);

x = [x0,zeros(n,max_iteration-1)];
x_hat_DSSO = zeros(n,max_iteration);

a_hat_DSSO = zeros(q,max_iteration);

for k = 1:max_iteration
    % system simulation and 
    y= C*x(:,k) + a;
    y_hat_DSSO = C*x_hat_DSSO(:,k) + a_hat_DSSO(:,k);
    x_hat_DSSO(:,k+1) = A*x_hat_DSSO(:,k) - L'*(y_hat_DSSO - y);
    a_hat_DSSO(:,k+1) = prox_l1(a_hat_DSSO(:,k) - nu*(y_hat_DSSO - y),nu*lambda);
    x(:,k+1) = A*x(:,k);
    
    % errors
    attack_support_error_DSSO(k) = sum((a_hat_DSSO(:,k+1)~=0) ~= (a~= 0));
    estimation_error_DSSO(k) = norm(x_hat_DSSO(:,k+1)-x(:,k+1))/norm(x(:,k+1));
    
end
%%
figure(1)
loglog(estimation_error_SSO,'b')
hold on, grid on
loglog(estimation_error_DSSO,'r')
title('estimation error')
xlabel('the number of iterations')
ylabel('estimation error')

figure(2)
semilogx(attack_support_error_SSO,'b')
hold on, grid on
semilogx(attack_support_error_DSSO,'r')
title('attack support error')
xlabel('the number of iterations')
ylabel('attack support error')

%% part 2
max_iteration = 1e5;

% SSO
nu = 0.99/(norm(G,2)^2);
lambda = 0.1;

% initialization
estimation_error_SSO = zeros(1,max_iteration);
attack_support_error_SSO = zeros(1,max_iteration);

x = [x0,zeros(n,max_iteration-1)];
x_hat_SSO = zeros(n,max_iteration);

a_hat_SSO = zeros(q,max_iteration);

for k = 1:1e3
    % system simulation and 
    y= C*x(:,k) + a;
    y_hat_SSO = C*x_hat_SSO(:,k) + a_hat_SSO(:,k);
    x_hat_SSO(:,k+1) = A*x_hat_SSO(:,k) - nu*A*C'*(y_hat_SSO - y);
    a_hat_SSO(:,k+1) = prox_l1(a_hat_SSO(:,k) - nu*(y_hat_SSO - y),nu*lambda);
    x(:,k+1) = A*x(:,k);
    
    % errors
    attack_support_error_SSO(k) = sum((a_hat_SSO(:,k+1)~=0) ~= (a~= 0));
    estimation_error_SSO(k) = norm(x_hat_SSO(:,k+1)-x(:,k+1))/norm(x(:,k+1));
    
end
last_k = k;

for k = last_k+1:max_iteration
    % system simulation and 
    y= C*x(:,k) + a;
    y_hat_SSO = C*x_hat_SSO(:,k) + a_hat_SSO(:,last_k);
    x_hat_SSO(:,k+1) = A*x_hat_SSO(:,k) ;
    x(:,k+1) = A*x(:,k);
    
    % errors
    attack_support_error_SSO(k) = sum((a_hat_SSO(:,last_k)~=0) ~= (a~= 0));
    estimation_error_SSO(k) = norm(x_hat_SSO(:,k+1)-x(:,k+1))/norm(x(:,k+1));
    
end
% D-SSO
nu = 0.7;
lambda = 0.1;

desired_lambda = 0.01*abs(rand(1,n));

L = place(A',C',desired_lambda);

% initialization
estimation_error_DSSO = zeros(1,max_iteration);
attack_support_error_DSSO = zeros(1,max_iteration);

x = [x0,zeros(n,max_iteration-1)];
x_hat_DSSO = zeros(n,max_iteration);

a_hat_DSSO = zeros(q,max_iteration);

for k = 1:1e3
    % system simulation and 
    y= C*x(:,k) + a;
    y_hat_DSSO = C*x_hat_DSSO(:,k) + a_hat_DSSO(:,k);
    x_hat_DSSO(:,k+1) = A*x_hat_DSSO(:,k) - L'*(y_hat_DSSO - y);
    a_hat_DSSO(:,k+1) = prox_l1(a_hat_DSSO(:,k) - nu*(y_hat_DSSO - y),nu*lambda);
    x(:,k+1) = A*x(:,k);
    
    % errors
    attack_support_error_DSSO(k) = sum((a_hat_DSSO(:,k+1)~=0) ~= (a~= 0));
    estimation_error_DSSO(k) = norm(x_hat_DSSO(:,k+1)-x(:,k+1))/norm(x(:,k+1));
    
end

last_k = k;
for k = last_k:max_iteration
    % system simulation and 
    y= C*x(:,k) + a;
    y_hat_DSSO = C*x_hat_DSSO(:,k) + a_hat_DSSO(:,k);
    x_hat_DSSO(:,k+1) = A*x_hat_DSSO(:,k) ;
    a_hat_DSSO(:,k+1) = prox_l1(a_hat_DSSO(:,k) - nu*(y_hat_DSSO - y),nu*lambda);
    x(:,k+1) = A*x(:,k);
    
    % errors
    attack_support_error_DSSO(k) = sum((a_hat_DSSO(:,last_k)~=0) ~= (a~= 0));
    estimation_error_DSSO(k) = norm(x_hat_DSSO(:,k+1)-x(:,k+1))/norm(x(:,k+1));
    
end

%%
figure(3)
loglog(estimation_error_SSO,'b')
hold on, grid on
loglog(estimation_error_DSSO,'r')
title('estimation error')
xlabel('the number of iterations')
ylabel('estimation error')

figure(4)
semilogx(attack_support_error_SSO,'b')
hold on, grid on
semilogx(attack_support_error_DSSO,'r')
title('attack support error')
xlabel('the number of iterations')
ylabel('attack support error')
