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
delta_my = 1e-6;
max_iteration = 1e4;

% simulation part
run = 5; % the number of simulations

% performance matrices initialization
last_col_IJAM = zeros(1,run);    
last_col_ISTA = zeros(1,run);
last_col_my = zeros(1,run);

state_estimation_error_IJAM = zeros(1,run);
state_estimation_error_ISTA = zeros(1,run);
state_estimation_error_my = zeros(1,run); 

attack_support_error_IJAM = zeros(1,run);
attack_support_error_ISTA = zeros(1,run);
attack_support_error_my = zeros(1,run);


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

    %--------------------------------------------------------------------

    % SOLVERS ***********************************************************

    %IJAM ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % parameter setting
    lambda = 0.1;
    nu = 0.7;
    % initialization
    x_hat_IJAM = zeros(n,max_iteration);
    a_hat_IJAM = zeros(q,max_iteration);

    for k = 1:max_iteration
        
        x_hat_IJAM(:,k+1) = pinv(C)*(y-a_hat_IJAM(:,k));

        a_hat_IJAM(:,k+1) = prox_l1(a_hat_IJAM(:,k) - nu*(C*x_hat_IJAM(:,k) ...
                            +a_hat_IJAM(:,k)-y),lambda);
        
        if(norm(x_hat_IJAM(:,k+1)-x_hat_IJAM(:,k))^2<delta)  % condition for maximum iteration number
            break
        end

    end
 
    last_col_IJAM(i) = find(any(x_hat_IJAM, 1), 1, 'last');

    figure(i)  
    subplot(3,1,1);
    for j = 1:length(x_tilde)
        plot(x_hat_IJAM(j,1:last_col_IJAM(i)))
        hold on,grid on
        plot(x_tilde(j)*ones(1,last_col_IJAM(i)),'--k')
        title('x hat IJAM')
    end

    
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

        if(norm(x_hat_ISTA(:,k+1)-x_hat_ISTA(:,k))^2<delta)  % condition for maximum iteration number
            break
        end
    end
    last_col_ISTA(i) = find(any(x_hat_ISTA, 1), 1, 'last'); %finding the last non-zero column

    % [[1:n]' x_tilde x_hat_ISTA(:,last_col_ISTA(i))]
    % [[1:q]' a_tilde a_hat_ISTA(:,last_col_ISTA(i))]
    
    
    figure(i)  
    subplot(3,1,2);
    for j = 1:length(x_tilde)
        plot(x_hat_ISTA(j,1:last_col_ISTA(i)))
        hold on,grid on
        plot(x_tilde(j)*ones(1,last_col_ISTA(i)),'--k')
        title('x hat ISTA')
    end
    
    % MIX +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % After a_hat became h-sparse, we simply make a support matrix 
    % for the attack, and then we solve the problem in one step
    % by least-square

    % delta_my = 1e-5; we don't need a very exact convergence
    lambda = 0.1;
    nu = 0.7;
    
    %initialization
    x_hat_my = zeros(n,max_iteration);
    a_hat_my = zeros(q,max_iteration);
    %try initializing x with pinv(G)*C
    
    for k = 1:max_iteration
        
        x_hat_my(:,k+1) = pinv(C)*(y-a_hat_my(:,k));

        a_hat_my(:,k+1) = prox_l1(a_hat_my(:,k) - nu*(C*x_hat_my(:,k) ...
                            +a_hat_my(:,k)-y),lambda);
        
        if(norm(x_hat_my(:,k+1)-x_hat_my(:,k))^2<delta_my)  % condition for maximum iteration number
            break
        end
        
    end
    
    % ATTENTION!!
    % still one step needs to be done so these numbers are incremented by 1
    % at the end of the process.
    last_col_my(i) = find(any(x_hat_my, 1), 1, 'last');
       

    
    non_zero_idx = find(a_hat_my(:,last_col_my(:,i)));
    
    
    % attack support matrix
    Ga = zeros(q); 
    for j = 1:length(non_zero_idx)
        Ga(non_zero_idx(j),non_zero_idx(j)) = 1;
    end
    
    G_sparse = [C, Ga]; %the position of the attacked is known now.
    z_sparse = pinv(G_sparse)*y;
    
    % the last time k+1 was filled in the table
    x_hat_my(:,k+2) = z_sparse(1:n);

    % attack needs to be cleaned due to numerical errors
    attack = z_sparse(n+non_zero_idx);
    for j = 1:length(non_zero_idx)
       a_hat_my(non_zero_idx(j),k+2) = attack(j);
    end
    % cleaning a due to numerical error
    
    % since one step is done outside of the for cycle
    last_col_my(i) = last_col_my(i) + 1;
    
    figure(i)  
    subplot(3,1,3);
    for j = 1:length(x_tilde)
        plot(x_hat_my(j,1:last_col_my(i)))
        hold on,grid on
        plot(x_tilde(j)*ones(1,last_col_my(i)),'--k')
        title('x hat my idea')
    end
    % states comparison (just for checking if it works)
    % [[1:n]' x_tilde x_hat_my(:,last_col_my(i)) x_hat_IJAM(:,last_col_IJAM(i))]
    % [[1:q]' a_tilde a_hat_my(:,last_col_my(i)) a_hat_IJAM(:,last_col_IJAM(i))]
    

    %---------------------------------------------------------------------
 
    % Calculating performance matrices +++++++++++++++++++++++++++++++++++
    state_estimation_error_IJAM(i) = norm(x_hat_IJAM(:,last_col_IJAM(i))-x_tilde)/norm(x_tilde);
    state_estimation_error_ISTA(i) = norm(x_hat_ISTA(:,last_col_ISTA(i))-x_tilde)/norm(x_tilde);
    state_estimation_error_my(i) = norm(x_hat_my(:,last_col_my(i))-x_tilde)/norm(x_tilde);
    

    attack_idx = find(a_tilde)
    attack_hat_IJAM_idx = find((a_hat_IJAM(:,last_col_IJAM(i))))
    attack_hat_ISTA_idx = find((a_hat_ISTA(:,last_col_ISTA(i)))) %%seems weired
    attack_hat_my_idx = find(abs(a_hat_my(:,last_col_my(i))))

    
    
    %IJAM SUPPORT ERROR
    for j = 1:h          
        if attack_idx(j) ~= attack_hat_IJAM_idx(j)
            attack_support_error_IJAM(i) =attack_support_error_IJAM(i) +1;  
        end                                      
    end
    
    %ISTA SUPPORT ERROR
    flag = 0;
    if(isempty(attack_hat_ISTA_idx))
            attack_support_error_ISTA(i)= attack_support_error_ISTA(i) +2;
            flag = 1;
    end

    for j = 1:h
    
        if flag == 1
            break
        end

        if(length(attack_idx) ==length(attack_hat_ISTA_idx))
            if attack_idx(j) ~= attack_hat_ISTA_idx(j)
                attack_support_error_ISTA(i) =attack_support_error_ISTA(i) +1;  
            end   
            
        else
            if attack_idx(j) ~= attack_hat_ISTA_idx
                attack_support_error_ISTA(i) =attack_support_error_ISTA(i) +1;  
            end
        end
        
    end

    %my SUPPORT ERROR
    for j = 1:h          
        if attack_idx(j) ~= attack_hat_my_idx(j)
            attack_support_error_my(i) =attack_support_error_my(i) +1;  
        end                                      
    end
    
    % take into account last_col_IJAM and last_col_ISTA
end


% Making average of the performance matrices

%average convergance comparison
mean_iteration_IJAM = mean(last_col_IJAM)
mean_iteration_ISTA = mean(last_col_ISTA)
mean_iteration_my = mean(last_col_my)

% average attack support error
mean_attack_support_error_IJAM = mean(attack_support_error_IJAM)
mean_attack_support_error_ISTA = mean(attack_support_error_ISTA)
mean_attack_support_error_my = mean(attack_support_error_my)

% average esetimation error
mean_state_estimation_error_IJAM = mean(state_estimation_error_IJAM)
mean_state_estimation_error_ISTA = mean(state_estimation_error_ISTA)
mean_state_estimation_error_my = mean(state_estimation_error_my)

my_precision_vs_IJAM = (mean_state_estimation_error_IJAM-mean_state_estimation_error_my)...
                        /mean_state_estimation_error_IJAM*100


