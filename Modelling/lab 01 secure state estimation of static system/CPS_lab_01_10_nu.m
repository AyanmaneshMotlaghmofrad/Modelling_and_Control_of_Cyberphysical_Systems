%secure state estimation
% in the case of A = cte
clc
clear variables
close all
format compact

%% setting the parameters

% System dimension information

n = 15; % the number of the states
q = 30; % the length of the output vector
h = 2; % two sensors are under attack.


% solver parameters

delta = 1e-10;
delta_my_IJAM = 1e-10;
delta_my_ISTA = 1e-10;
max_iteration = 1e5;

lambda_vec = [0.1, 0.4, 0.7, 1,3,5]; % Lambda values

% Create legend labels for lambda values
legendLabels = arrayfun(@(x) sprintf('\\lambda = %.1f', x), lambda_vec, 'UniformOutput', false);

% simulation part
run = 50; % the number of simulations
for lam = 1:length(lambda_vec)
    lambda = lambda_vec(lam);
    % performance matrices initialization
    last_col_IJAM = zeros(1,run);    
    last_col_ISTA = zeros(1,run);
    last_col_my_IJAM = zeros(1,run);
    last_col_my_ISTA = zeros(1,run);
    
    state_estimation_error_IJAM = zeros(run,max_iteration);
    state_estimation_error_ISTA = zeros(run,max_iteration);
    state_estimation_error_my_IJAM = zeros(run,max_iteration); 
    state_estimation_error_my_ISTA = zeros(run,max_iteration); 
    
    attack_support_error_IJAM = zeros(run,max_iteration);
    attack_support_error_ISTA = zeros(run,max_iteration);
    attack_support_error_my_IJAM = zeros(run,max_iteration);
    attack_support_error_my_ISTA = zeros(run,max_iteration);
    
    
    for i = 1:run  
    
        %Setting the system++++++++++++++++++++++++++++++++++++++++++++++++++
        % making C according to normal distribution
    
        C = randn(q,n);
        G = [C,eye(q)];
        % support of the attack vector a: uniform distribution
        % and "real" attack 
    
        a_tilde = zeros(q,1);
        h_count = h;
        
        while 1
            
            idx = randi([1,q]);
        
            if (a_tilde(idx) == 0)
                
                %chossing the number
                side = randi([1,2]);
                if (side == 1)
                    a_tilde(idx,1) = rand + 2;
                    
                    h_count = h_count-1;
                    if h_count== 0
                        break
                    end
        
                else
                    a_tilde(idx,1) = rand - 5;
                
                    h_count = h_count-1;
                    if h_count== 0
                        break
                    end
        
                end
                % the number is chosen
            end
        
        end
        % the real attack is defined here
    
        % "real" state
        x_tilde = zeros(n,1);
        
        for j = 1:n
            
            side = randi([1,2]);
            if (side == 1)
                x_tilde(j,1) = rand + 2;
            else
                x_tilde(j,1) = rand - 3;
            end
        
        end
        % real state is defined here
        
        % measurement noise 
       
        eta = 10^-2*randn(q,1);
    
        
        % y, the measurement vector 
        % corrupted by noise and attack
    
        y = C*x_tilde + a_tilde + eta;
        
        %checking the condition for h attack correctability
        % since the system is stable this can be checked
        zero_norm = length(find(C*x_tilde));
    
        if( zero_norm < 2*h+1)
            break
        end
    
    
        
        %--------------------------------------------------------------------
    
        % SOLVERS ***********************************************************
    
        %IJAM ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % parameter setting
        % lambda = 0.1;
        nu = 0.7;
        % initialization
        x_hat_IJAM = zeros(n,max_iteration);
        a_hat_IJAM = zeros(q,max_iteration);
    
        for k = 1:max_iteration
            
            x_hat_IJAM(:,k+1) = pinv(C)*(y-a_hat_IJAM(:,k));
    
            a_hat_IJAM(:,k+1) = prox_l1(a_hat_IJAM(:,k) - nu*(C*x_hat_IJAM(:,k) ...
                                +a_hat_IJAM(:,k)-y),lambda*nu);
            
            
            %calculating errors
            state_estimation_error_IJAM(i,k+1) = norm(x_hat_IJAM(:,k+1)-x_tilde)/norm(x_tilde); 
            attack_support_error_IJAM(i,k+1) = sum((a_hat_IJAM(:, k+1) ~= 0) ~= (a_tilde ~= 0));
    
            if(norm(x_hat_IJAM(:,k+1)-x_hat_IJAM(:,k))^2<delta)  % condition for maximum iteration number
                break
            end
    
        end
     
        last_col_IJAM(i) = k;
        attack_support_error_IJAM(i,k+2:end)=attack_support_error_IJAM(i,k+1);
        state_estimation_error_IJAM(i,k+2:end) = state_estimation_error_IJAM(i,k+1);
        
        %ISTA ++++++++++++++++++++++++++++++++++++++++++++++++++
        nu = 0.99/(norm(G,2)^2);
        
        %initialization
        x_hat_ISTA = zeros(n,max_iteration);
        a_hat_ISTA = zeros(q,max_iteration);
    
        for k = 1:max_iteration
            
            x_hat_ISTA(:,k+1) = x_hat_ISTA(:,k) - ... 
                                nu*C'*(C*x_hat_ISTA(:,k) + a_hat_ISTA(:,k)- y);
          
            a_hat_ISTA(:,k+1) = prox_l1(a_hat_ISTA(:,k) - nu*(C*x_hat_ISTA(:,k) ...
                                 +a_hat_ISTA(:,k)-y),lambda*nu);
            
            %calculating errors
            state_estimation_error_ISTA(i,k+1) = norm(x_hat_ISTA(:,k+1)-x_tilde)/norm(x_tilde); 
            attack_support_error_ISTA(i,k+1) = sum((a_hat_ISTA(:, k+1) ~= 0) ~= (a_tilde ~= 0));
    
            if(norm(x_hat_ISTA(:,k+1)-x_hat_ISTA(:,k))^2<delta)  % condition for maximum iteration number
                break
            end
        end
        last_col_ISTA(i) = k; %finding the last non-zero column
        attack_support_error_ISTA(i,k+2:end)=attack_support_error_ISTA(i,k+1);
        state_estimation_error_ISTA(i,k+2:end) = state_estimation_error_ISTA(i,k+1);
        
        % [[1:n]' x_tilde x_hat_ISTA(:,last_col_ISTA(i))]
        % [[1:q]' a_tilde a_hat_ISTA(:,last_col_ISTA(i))]
        
        
    
        % MIX my IJAM +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        % delta_my_IJAM = 1e-5; we don't need a very exact convergence
        % lambda = 0.1;
        nu = 0.7;
        
        %initialization
        x_hat_my_IJAM = zeros(n,max_iteration);
        a_hat_my_IJAM = zeros(q,max_iteration);
        %try initializing x with pinv(G)*C
        
        for k = 1:max_iteration
            
            x_hat_my_IJAM(:,k+1) = pinv(C)*(y-a_hat_my_IJAM(:,k));
    
            a_hat_my_IJAM(:,k+1) = prox_l1(a_hat_my_IJAM(:,k) - nu*(C*x_hat_my_IJAM(:,k) ...
                                +a_hat_my_IJAM(:,k)-y),lambda*nu);
            
            %calculating errors
            state_estimation_error_my_IJAM(i,k+1) = norm(x_hat_my_IJAM(:,k+1)-x_tilde)/norm(x_tilde); 
            attack_support_error_my_IJAM(i,k+1) = sum((a_hat_my_IJAM(:,k+1)~=0) ~= (a_tilde~= 0));
    
            if(norm(x_hat_my_IJAM(:,k+1)-x_hat_my_IJAM(:,k))^2<delta_my_IJAM)  % condition for maximum iteration number
                break
            end
            
        end
        
        % ATTENTION!!
        % still one step needs to be done so these numbers are incremented by 1
        % at the end of the process.
        last_col_my_IJAM(i) = k;
           
    
        
        non_zero_idx = find(a_hat_my_IJAM(:,last_col_my_IJAM(i)));
        
        
        % attack support matrix
        Ga = zeros(q); 
        for j = 1:length(non_zero_idx)
            Ga(non_zero_idx(j),non_zero_idx(j)) = 1;
        end
        
        G_sparse = [C, Ga]; %the position of the attacked is known now.
        z_sparse = pinv(G_sparse)*y;
        
        % the last time k+1 was filled in the table
        x_hat_my_IJAM(:,k+2) = z_sparse(1:n);
        
        % attack needs to be cleaned due to numerical errors
        attack = z_sparse(n+non_zero_idx);
        for j = 1:length(non_zero_idx)
           a_hat_my_IJAM(non_zero_idx(j),k+2) = attack(j);
        end
        % cleaning 'a' due to numerical error
        
        % since one step is done outside of the for cycle
        last_col_my_IJAM(i) = last_col_my_IJAM(i) + 1;
        
        state_estimation_error_my_IJAM(i,k+2) = norm(x_hat_my_IJAM(:,k+2)-x_tilde)/norm(x_tilde); 
        state_estimation_error_my_IJAM(i,k+3:end) = state_estimation_error_my_IJAM(i,k+2); 
    
        
        attack_support_error_my_IJAM(i,k+2) = sum((a_hat_my_IJAM(:,k+2)~=0) ~= (a_tilde~= 0));
        attack_support_error_my_IJAM(i,k+3:end)=attack_support_error_my_IJAM(i,k+2);
    
        
        % semilogx(state_estimation_error_my_IJAM(i,:),'r')
        
        % states comparison (just for checking if it works)
         % [[1:n]' x_tilde x_hat_my_IJAM(:,last_col_my_IJAM(i)) x_hat_IJAM(:,last_col_IJAM(i))]
         % [[1:q]' a_tilde a_hat_my_IJAM(:,last_col_my_IJAM(i)) a_hat_IJAM(:,last_col_IJAM(i))]
      
        % MIX my ISTA +++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % delta_my_ISTA = 1e-5; we don't need a very exact convergence
       
        % lambda = 0.1;
        nu = 0.99/(norm(G,2)^2);
        %initialization
        x_hat_my_ISTA = zeros(n,max_iteration);
        a_hat_my_ISTA = zeros(q,max_iteration);
        %try initializing x with pinv(G)*C
        
        for k = 1:max_iteration
            
            x_hat_my_ISTA(:,k+1) = x_hat_my_ISTA(:,k) - ... 
                                nu*C'*(C*x_hat_my_ISTA(:,k) + a_hat_my_ISTA(:,k)- y);
          
            a_hat_my_ISTA(:,k+1) = prox_l1(a_hat_my_ISTA(:,k) - nu*(C*x_hat_my_ISTA(:,k) ...
                                 +a_hat_my_ISTA(:,k)-y),lambda*nu);
    
            state_estimation_error_my_ISTA(i,k+1) = norm(x_hat_my_ISTA(:,k+1)-x_tilde)/norm(x_tilde); 
            attack_support_error_my_ISTA(i,k+1) = sum((a_hat_my_ISTA(:, k+1) ~= 0) ~= (a_tilde ~= 0));
    
            if(norm(x_hat_my_ISTA(:,k+1)-x_hat_my_ISTA(:,k))^2<delta_my_ISTA)  % condition for maximum iteration number
                break
            end
            
        end
        
        % ATTENTION!!
        % still one step needs to be done so these numbers are incremented by 1
        % at the end of the process.
        last_col_my_ISTA(i) = k;
           
        non_zero_idx = find(a_hat_my_ISTA(:,last_col_my_ISTA(:,i)));
        
        
        % attack support matrix
        Ga = zeros(q); 
        for j = 1:length(non_zero_idx)
            Ga(non_zero_idx(j),non_zero_idx(j)) = 1;
        end
        
        G_sparse = [C, Ga]; %the position of the attacked is known now.
        z_sparse = pinv(G_sparse)*y;
        
        % the last time k+1 was filled in the table
        x_hat_my_ISTA(:,k+2) = z_sparse(1:n);
    
        % attack needs to be cleaned due to numerical errors
        attack = z_sparse(n+non_zero_idx);
        for j = 1:length(non_zero_idx)
           a_hat_my_ISTA(non_zero_idx(j),k+2) = attack(j);
        end
        % cleaning a due to numerical error
        
        % since one step is done outside of the for cycle
        last_col_my_ISTA(i) = last_col_my_ISTA(i) + 1;
       
        state_estimation_error_my_ISTA(i,k+2) = norm(x_hat_my_ISTA(:,k+2)-x_tilde)/norm(x_tilde); 
        state_estimation_error_my_ISTA(i,k+3:end) = state_estimation_error_my_ISTA(i,k+2); 
        
        attack_support_error_my_ISTA(i,k+2) = sum((a_hat_my_ISTA(:,k+2)~=0) ~= (a_tilde~= 0));
        attack_support_error_my_ISTA(i,k+3:end)=attack_support_error_my_ISTA(i,k+2);
    
        
     
        % states comparison (just for checking if it works)
        % [[1:n]' x_tilde x_hat_my_ISTA(:,last_col_my_ISTA(i)) x_hat_ISTA(:,last_col_ISTA(i))]
        % [[1:q]' a_tilde a_hat_my_ISTA(:,last_col_my_ISTA(i)) a_hat_ISTA(:,last_col_ISTA(i))]
        
    
        %---------------------------------------------------------------------
     
    end
    
    
    % Making average of the performance matrices
    
    %average convergance comparison
    mean_iteration_IJAM(lam) = mean(last_col_IJAM);
    mean_iteration_ISTA(lam)  = mean(last_col_ISTA);
    mean_iteration_my_IJAM(lam)  = mean(last_col_my_IJAM);
    mean_iteration_my_ISTA(lam)  = mean(last_col_my_ISTA);
    
    % average attack support error
    mean_attack_support_error_IJAM = mean(attack_support_error_IJAM,1);
    mean_attack_support_error_ISTA = mean(attack_support_error_ISTA,1);
    mean_attack_support_error_my_IJAM = mean(attack_support_error_my_IJAM,1);
    mean_attack_support_error_my_ISTA = mean(attack_support_error_my_ISTA,1);
    
    % average estimation error
    mean_state_estimation_error_IJAM = mean(state_estimation_error_IJAM,1);
    mean_state_estimation_error_ISTA = mean(state_estimation_error_ISTA,1);
    mean_state_estimation_error_my_IJAM = mean(state_estimation_error_my_IJAM,1);
    mean_state_estimation_error_my_ISTA = mean(state_estimation_error_my_ISTA,1);
    
    % % precision comparison
    % last_non_zero_IJAM = find(mean_state_estimation_error_IJAM ~= 0, 1, 'last');
    % last_non_zero_ISTA = find(mean_state_estimation_error_ISTA ~= 0, 1, 'last');
    % last_non_zero_my_IJAM = find(mean_state_estimation_error_my_IJAM ~= 0, 1, 'last');
    % last_non_zero_my_ISTA = find(mean_state_estimation_error_my_ISTA ~= 0, 1, 'last');
    
    
    % error_my_IJAM(counter) = mean_attack_support_error_my_IJAM(last_non_zero_my_IJAM)
    error_IJAM(lam) = mean_attack_support_error_IJAM(end)
    % error_my_ISTA(counter) = mean_attack_support_error_my_ISTA(last_non_zero_my_ISTA)
    error_ISTA(lam) = mean_attack_support_error_ISTA(end)
    
    
    
end
%% plotting
% 
% figure(1)  
% semilogx(mean_state_estimation_error_IJAM,'m')
% hold on,grid on
% semilogx(mean_state_estimation_error_ISTA,'g')
% semilogx(mean_state_estimation_error_my_IJAM,'b')
% semilogx(mean_state_estimation_error_my_ISTA,'r')
% 
% title('Estate Estimation Error')
% ylabel( 'mean relative ℓ₂ norm error')
% xlabel('iteration number')
% legend('IJAM','ISTA','my IJAM', 'my ISTA')
% 
% figure(2)  
% semilogx(mean_attack_support_error_IJAM,'m')
% hold on,grid on
% semilogx(mean_attack_support_error_ISTA,'g')
% % semilogx(mean_attack_support_error_my_IJAM,'k')
% % semilogx(mean_attack_support_error_my_ISTA,'r')
% 
% title('Support Attack Error')
% ylabel( 'mean support attack error')
% xlabel('iteration number')
% legend('IJAM','ISTA')


%% plotting for resilience toward attack

% step 1: counter should be set to 1
% step 2: algorithm should be run for different values of attack
% step 3: these lines should be used in order to plot 


mean_iterations = [error_IJAM', error_ISTA'];

% Labels for methods
methods = {'IJAM', 'ISTA'};

% Plot grouped bar chart
bar_handle = bar(1:length(lambda_vec), mean_iterations, 'grouped');

colors = [0.2 0.6 0.8;  % Soft Cyan 
          %0.8 0.4 0.2;  % Warm Orange
         % 0.4 0.2 0.6;  % Purple 
          0.2 0.8 0.4]; % Green 

% Apply colors to each bar group
for k = 1:length(bar_handle)
    bar_handle(k).FaceColor = colors(k, :);
end

xticks(1:length(lambda_vec));                % Set tick positions
xticklabels(string(lambda_vec));              % Set tick labels to lambda values

xlabel('\lambda') 
ylabel('mean attack support error') 
title('Mean Attak Support Error Vs. \lambda')
legend(methods, 'Location', 'northwest')
grid on

%% iteration

mean_iterations = [mean_iteration_my_IJAM', mean_iteration_my_ISTA'];

% Labels for methods
methods = {'IJAM', 'ISTA'};

% Plot grouped bar chart
bar_handle = bar(1:length(lambda_vec), mean_iterations, 'grouped');

colors = [0.2 0.6 0.8;  % Soft Cyan 
          %0.8 0.4 0.2;  % Warm Orange
         % 0.4 0.2 0.6;  % Purple 
          0.2 0.8 0.4]; % Green 

% Apply colors to each bar group
for k = 1:length(bar_handle)
    bar_handle(k).FaceColor = colors(k, :);
end

xticks(1:length(lambda_vec));                % Set tick positions
xticklabels(string(lambda_vec));              % Set tick labels to lambda values

xlabel('\lambda') 
ylabel('mean iteration') 
title('Number of iterations Vs. \lambda')
legend(methods, 'Location', 'northwest')
grid on
