clc
clear variables
close all
format compact

%% loading the data

load localization_data.mat

size_D = size(D);
q = size_D(1);
n = size_D(2);

%solution
target_position = 37;
x_tilde = zeros(n,1);
x_tilde(target_position) = 1;

a_support_idx = [1,10,14,16,17];
a_tilde_support = zeros(q,1);
a_tilde_support(a_support_idx) = 1;


% solving the problem
% solver parameters

delta = 1e-10;
delta_my_ISTA = 1e-6;
max_iteration = 1e4;


% to avoid numerical problems component of 
% D have large values compared to I.
%normalizing D

G = normalize([D, eye(q)]);

                       
%ISTA for localization full LASSO +++++++++++++++++++++++++++++++++++++++++
nu = 2.33/(norm(G,2)^2);
lambda_1 = 14;
lambda_2 = 11.4;

% lambda_1 = 10;
% lambda_2 = 10;
%initialization
z_hat_ISTA = zeros(n + q,max_iteration);


for k = 1:max_iteration
    
    z_hat_ISTA(:,k+1) = prox_l1_vec(z_hat_ISTA(:,k) - ... 
                        nu*G'*(G*z_hat_ISTA(:,k)-y),[nu*lambda_1;nu*lambda_2],n,q);
  
    if(norm(z_hat_ISTA(:,k+1)-z_hat_ISTA(:,k))^2<delta)  % condition for maximum iteration number
        break
    end
end
last_col_ISTA = find(any(z_hat_ISTA, 1), 1, 'last'); %finding the last non-zero column

% cleaning the data %the x_i  =  {0, 1}
x_hat_ISTA = z_hat_ISTA(1:n,last_col_ISTA);
a_hat_ISTA = z_hat_ISTA(n+1:end,last_col_ISTA);

non_zero_x = find(x_hat_ISTA);
non_zero_a = find(a_hat_ISTA);

cleaning_threshold = 0.4;

% the state of the cell can be {0,1}
% so we have to clean the data
for i = 1:length(non_zero_x)
    if(x_hat_ISTA(non_zero_x(i))<cleaning_threshold)
        x_hat_ISTA(non_zero_x(i)) = 0;
    
    else
        x_hat_ISTA(non_zero_x(i)) = 1;
    end

end
target_position = find(x_hat_ISTA)
% the check whether position of the attacks are detected correctly.
% [a_support_idx' non_zero_a]

% identifying the attack
a_support_matrix = zeros(q);

for i = 1:length(non_zero_a)
    a_support_matrix(non_zero_a(i),non_zero_a(i)) = 1;
end

a_hat = pinv(a_support_matrix)*(y-D*x_hat_ISTA);

% just to check the index of the attacks
[find(a_hat) find(a_hat_ISTA)]



% now that the position of the target is clear, we can solve the 
% problem for identifying the attack

%% plotting

figure(1)
bar(1:length(find(a_hat)),a_hat(find(a_hat)),'g')

xticks([1, 2, 3, 4, 5]);  % Specify the x-tick positions
xticklabels({'1', '10', '14', '16','17'});  % Set the labels for the x-ticks
xlabel('attacked sensor')
ylabel('the value of the attack')
title('The value of the attack corrupting the measurement')                    

%% finding the position of the sensors

[~,idx] = max(D,[],1)
