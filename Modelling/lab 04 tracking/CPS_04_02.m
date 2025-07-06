% tracking
clc
clear variables 
close all
format compact

%%

load tracking_data.mat

n = size(A,1);
q = size(D,1);
%%
G = normalize([D eye(q)]);

% checking observability
% without attack
rank(obsv(A,D))
% the system is observable full rank

% with constantattack
A_with_attack = [A,zeros(n,q);zeros(q,n) eye(q)];
rank(obsv(A_with_attack,G)) % the system with attack is not observable
D1 = G(:,1:n);
I = G(:,n+1:end);
% the system with attack is not observable so a classic Luenburge cannot


%% checking observability
% without attack
rank(obsv(A,D))
% the system is observable full rank

% with constantattack
A_with_attack = [A,zeros(n,q);zeros(q,n) eye(q)];
rank(obsv(A_with_attack,G))
%% SSO for tracking
max_iteration = 80;

% SSO
nu1 = 0.99/(norm(D1,2)^2);
nu2 = 0.99/(norm(I,2)^2);
lambda1 = 10;
lambda2 = 10;
% initialization
state_support_error_SSO = zeros(1,max_iteration);
attack_support_error_SSO = zeros(1,max_iteration);

x = [xtrue0,zeros(n,max_iteration-1)];
x_hat_SSO = zeros(n,max_iteration);
a_hat_SSO = zeros(q,max_iteration);

figure(1)
hold on
for k = 1:max_iteration
    % system simulation and 
    y= D1*x(:,k) + atrue;
    y_hat_SSO = D1*x_hat_SSO(:,k) + a_hat_SSO(:,k);
    x_hat_SSO(:,k+1) = prox_l1(A*x_hat_SSO(:,k) - nu1*A*D1'*(y_hat_SSO - y),nu1*lambda1);
    a_hat_SSO(:,k+1) = prox_l1(a_hat_SSO(:,k) - nu2*I'*(y_hat_SSO - y),nu2*lambda2);
    x(:,k+1) = A*x(:,k);
    % Finding position
    non_zero_x = find(x_hat_SSO(:,k+1));
    non_zero_a = find(a_hat_SSO(:,k+1));

    cleaning_threshold = 0.1;
    position = zeros(n,1);
    % the state of the cell can be {0,1}
    % so we have to clean the data
    for i = 1:length(non_zero_x)
        if(x_hat_SSO(non_zero_x(i),k+1)<cleaning_threshold)
            position(non_zero_x(i)) = 0;      
        else
            position(non_zero_x(i)) = 1;
        end

    end
    % k
    % find(position)
    % find(x(:,k+1))

    % errors
    attack_support_error_SSO(k) = sum((a_hat_SSO(:,k+1)~=0) ~= (atrue~= 0));
    state_support_error_SSO(k) = sum((position~=0) ~= (x(:,k+1)~= 0));
    
    if(k>10)
    %plot
           plot(k,find(x(:,k+1)),'+b')
           plot(k,find(position),'+r')
    end
end

grid on
title('Tracking Representation SSO')
xlabel('number of iterations')
ylabel('the number of the cells')
%%
figure(2)
plot(state_support_error_SSO,'b')
hold on, grid on
title('state support error SSO')
xlabel('the number of iterations')
ylabel('state support error')


figure(3)
plot(attack_support_error_SSO,'b')
hold on, grid on
title('attack support error SSO')
xlabel('the number of iterations')
ylabel('attack support error')



%% Try with DSSO?

% whether system without attack is observable?
% Since the attacks are constant we can consider the dynamic of the attacks
% that is encoded by an identity matrix

% reagrding observability we can do some manual analysis
% like the one for luenburge

% We can try with deadbeat

%% DSSO for tracking
max_iteration = 230;

% DSSO
nu1 = 0.05;
lambda1 = 0.01;

nu2 = 0.05;
lambda2 = 0.01;

% initialization
state_support_error_DSSO = zeros(1,max_iteration);
attack_support_error_DSSO = zeros(1,max_iteration);

x = [xtrue0,zeros(n,max_iteration-1)];
x_hat_DSSO = zeros(n,max_iteration);
a_hat_DSSO = zeros(q,max_iteration);

desired_lambda = 0.01*abs(rand(1,n));
L = place(A',D',desired_lambda)';

figure(4)
hold on
for k = 1:max_iteration
    % system simulation and 
    y= D*x(:,k) + atrue;
    y_hat_DSSO = D*x_hat_DSSO(:,k) + a_hat_DSSO(:,k);
    x_hat_DSSO(:,k+1) = prox_l1(A*x_hat_DSSO(:,k) - L*(y_hat_DSSO - y),nu1*lambda1);
    a_hat_DSSO(:,k+1) = prox_l1(a_hat_DSSO(:,k) - nu2*(y_hat_DSSO - y),nu2*lambda2);
    x(:,k+1) = A*x(:,k);
    
    
    % Finding position
    non_zero_x = find(x_hat_DSSO(:,k+1));
    non_zero_a = find(a_hat_DSSO(:,k+1));

    cleaning_threshold = 0.6;
    position = zeros(n,1);
    % the state of the cell can be {0,1}
    % so we have to clean the data
    for i = 1:length(non_zero_x)
        if(x_hat_DSSO(non_zero_x(i),k+1)<cleaning_threshold)
            position(non_zero_x(i)) = 0;      
        else
            position(non_zero_x(i)) = 1;
        end

    end
    k;
    find(position);
    find(x(:,k+1));
    
    % errors
    attack_support_error_DSSO(k) = sum((a_hat_DSSO(:,k+1)~=0) ~= (atrue~= 0));
    state_support_error_DSSO(k) = sum((position~=0) ~= (x(:,k+1)~= 0));
    
    if(k>180)
        %plot
        plot(k,find(x(:,k+1)),'+b')
        plot(k,find(position),'+r')
    end
end

grid on
title('Tracking Representation DSSO')
xlabel('number of iterations')
ylabel('the number of the cells')
%%
figure(5)
plot(state_support_error_DSSO,'r')
hold on, grid on
title('state support error DSSO')
xlabel('the number of iterations')
ylabel('state support error')

figure(6)
plot(attack_support_error_DSSO,'r')
hold on, grid on
title('attack support error DSSO')
xlabel('the number of iterations')
ylabel('attack support error')


