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

%% checking observability
obsv_rank = rank(obsv(A,C))

% Observability matrix is full rank, so all the states are observable
% without attack

% With attack the system wouldn't be observable
rank(obsv([A zeros(n,q);zeros(q,n) eye(q)],G))
%% part 1
max_iteration = 1e5;

% SSO
nu = 0.99/(norm(G,2)^2)*2;
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
legend('SSO', 'DSSO')


figure(2)
semilogx(attack_support_error_SSO,'b')
hold on, grid on
semilogx(attack_support_error_DSSO,'r')
title('attack support error')
xlabel('the number of iterations')
ylabel('attack support error')
legend('SSO', 'DSSO')


%% part 2
max_iteration = 1e5;

%IDEAS to be implemented in part 3 and 4
% 1) once the position of the attacks are found put away those sensors
%    and go on with SSO and D-SSO

% 2) change the states, use SSO or D-SSO for the new states corresponding to
%  the dcgain mode, then use a Luenberg observer for the other modes. 

% Here the input of the observer is cut after some iterations, and 
% surprisingly, the precision of the estimation is enhanced.

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
legend('SSO', 'DSSO')


figure(4)
semilogx(attack_support_error_SSO,'b')
hold on, grid on
semilogx(attack_support_error_DSSO,'r')
title('attack support error')
xlabel('the number of iterations')
ylabel('attack support error')
legend('SSO', 'DSSO')

%% part 3
max_iteration = 1e5;

%IDEAS
% 1) once the position of the attacks are found put away those sensors
%    and go on with SSO and D-SSO

nu = 0.99/(norm(G,2)^2);
lambda = 0.1;

% initialization
delta = 1e-5;

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

    if (x_hat_SSO(:,k+1)-x_hat_SSO(:,k) < delta)
        break
    end

    
end
last_k = k;
C_new_SSO = C;
C_new_SSO(find(a_hat_SSO(:,k)),:) = [];

if (rank(obsv(A,C_new_SSO)) ~= n)
     %some of the states are not going to be observable,
     % unobservable states are going to be separated and
     % the observer should run on the observable states.
     fprintf('Not all the states are observable')
end

for k = last_k+1:max_iteration
    % system simulation and 
    y= C_new_SSO*x(:,k);
    y_hat_SSO = C_new_SSO*x_hat_SSO(:,k);
    if((y_hat_SSO - y)<1e-8)
        x_hat_SSO(:,k+1) = A*x_hat_SSO(:,k);
    else
        x_hat_SSO(:,k+1) = A*x_hat_SSO(:,k) - nu*A*C_new_SSO'*(y_hat_SSO - y);
    end
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
    
    if(x_hat_DSSO(:,k+1)-x_hat_DSSO(:,k) < delta)
        break
    end

end

last_k = k;
% removing attacked sensors
C_new_DSSO = C;
C_new_DSSO(find(a_hat_DSSO(:,k)),:) = [];

if (rank(obsv(A,C_new_DSSO)) ~= n)
     %some of the states are not going to be observable,
     % unobservable states are going to be separated and
     % the observer should run on the observable states.
     fprintf('Not all the states are observable')
end

% placing new eigen values
desired_lambda = 0.01*abs(rand(1,n));

L = place(A',C_new_DSSO',desired_lambda);

for k = last_k:max_iteration
    % system simulation and 
    y= C_new_DSSO*x(:,k);
    y_hat_DSSO = C_new_DSSO*x_hat_DSSO(:,k);
    if ((y_hat_DSSO - y)<1e-8)
        x_hat_DSSO(:,k+1) = A*x_hat_DSSO(:,k);
    else
        x_hat_DSSO(:,k+1) = A*x_hat_DSSO(:,k) - L'*(y_hat_DSSO - y);
    end
    x(:,k+1) = A*x(:,k);
    
    % errors
    attack_support_error_DSSO(k) = sum((a_hat_DSSO(:,last_k)~=0) ~= (a~= 0));
    estimation_error_DSSO(k) = norm(x_hat_DSSO(:,k+1)-x(:,k+1))/norm(x(:,k+1));
    
end

%%
figure(5)
loglog(estimation_error_SSO,'b')
hold on, grid on
loglog(estimation_error_DSSO,'r')
title('estimation error')
xlabel('the number of iterations')
ylabel('estimation error')
legend('SSO', 'DSSO')

figure(6)
semilogx(attack_support_error_SSO,'b')
hold on, grid on
semilogx(attack_support_error_DSSO,'r')
title('attack support error')
xlabel('the number of iterations')
ylabel('attack support error')
legend('SSO', 'DSSO')
