clc
clear variables
close all
format compact

%% loading the data

load localization_data.mat

size = size(D);
q = size(1);
n = size(2);

%solution
target_position = 37;
x_tilde = zeros(n,1);
x_tilde(target_position) = 1;

a_support_idx = [1,10,14,16,17];
a_tilde_support = zeros(q,1);
a_tilde_support(a_support_idx) = 1;


%% solving the problem
% solver parameters

delta = 1e-10;
delta_my_ISTA = 1e-6;
max_iteration = 1e4;

x_hat_ISTA = zeros(n,max_iteration);
x_hat_my_ISTA = zeros(n,max_iteration);


% to avoid numerical problems component of 
% D have large values compared to I.
%normalizing D
D = normalize(D);
G = [D, eye(q)];

                       
%ISTA for localization full LASSO +++++++++++++++++++++++++++++++++++++++++
nu = 0.99/(sqrt(norm(G,2)));
lambda_1 = 10;
lambda_2 = 10;

%initialization
x_hat_ISTA = zeros(n,max_iteration);
a_hat_ISTA = zeros(q,max_iteration);

for k = 1:max_iteration
    
    x_hat_ISTA(:,k+1) = prox_l1(x_hat_ISTA(:,k) - ... 
                        nu*D'*(D*x_hat_ISTA(:,k) + a_hat_ISTA(:,k)- y),nu*lambda_1);
  
    a_hat_ISTA(:,k+1) = prox_l1(a_hat_ISTA(:,k) - nu*(D*x_hat_ISTA(:,k) ...
                         +a_hat_ISTA(:,k)-y),nu*lambda_2);
    if(norm(x_hat_ISTA(:,k+1)-x_hat_ISTA(:,k))^2<delta)  % condition for maximum iteration number
        break
    end
end
last_col_ISTA = find(any(x_hat_ISTA, 1), 1, 'last'); %finding the last non-zero column

[[1:n]' x_tilde x_hat_ISTA(:,last_col_ISTA)]

[[1:q]' a_tilde_support a_hat_ISTA(:,last_col_ISTA)]

return
figure(1)  
for j = 1:length(x_tilde)
    plot(x_hat_ISTA(j,1:last_col_ISTA))
    hold on,grid on
end                         


