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
max_iteration = 1e4;

% simulation part
run = 20; % the number of simulations

% performance matrices initialization
last_col_x_IJAM = zeros(1,run);    
last_col_x_ISTA = zeros(1,run);

last_col_a_IJAM = zeros(1,run);    
last_col_a_ISTA = zeros(1,run);

state_estimation_error_IJAM = zeros(1,run);
state_estimation_error_ISTA = zeros(1,run);

attack_support_error_IJAM = zeros(1,run);
attack_support_error_ISTA = zeros(1,run);
    
for i = 1:run  

    %Setting the system+++++++++++++++++++++++++++++
    % making C according to normal distribution

    C = randn(q,n);
    G = [C,eye(q)];
    % support of the attack vector a: uniform distribution
    % and "real" attack +++++++++++++++++++++++++++++++++++++++++++++++++

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
    % the real attack is defined here-------------------------------

    % "real" state++++++++++++++++++++++++++++++++++++++++++++++++++++
    x_tilde = zeros(n,1);
    
    for j = 1:n
        
        side = randi([1,2]);
        if (side == 1)
            x_tilde(j,1) = rand + 2;
        else
            x_tilde(j,1) = rand - 3;
        end
    
    end
    % real state is defined here -------------------------------------
    
    % measurement noise  +++++++++++++++++++++++++++++++++++++++++++++++
   
    eta = 10^-2*randn(q,1);

    %-------------------------------------------------------------------
    
    % y, the measurement vector ++++++++++++++++++++++++++++++++++++++++
    % corrupted by noise and attack

    y = C*x_tilde + a_tilde + eta;

    %--------------------------------------------------------------------
    % SOLVERS ++++++++++++++++++++++++++++++++++++++++
    % INITIALIZING THE VARIABLES FOR THE SOLVER 

    
    
    
    %initialization
    x_hat_IJAM = zeros(n,max_iteration);
    x_hat_ISTA = zeros(n,max_iteration);
    
    a_hat_IJAM = zeros(q,max_iteration);
    a_hat_ISTA = zeros(q,max_iteration);
    


    %IJAM ++++++++++++++++++++++++++++++++++++++++++++++++++
    lambda = 0.1;
    nu = 0.7;
    
    for k = 1:max_iteration
        
        x_hat_IJAM(:,k+1) = pinv(C)*(y-a_hat_IJAM(:,k));

        a_hat_IJAM(:,k+1) = prox_l1(a_hat_IJAM(:,k) - nu*(C*x_hat_IJAM(:,k) ...
                            +a_hat_IJAM(:,k)-y),lambda);
        
        if(norm(x_hat_IJAM(:,k+1)-x_hat_IJAM(:,k))^2<delta)  % condition for maximum iteration number
            break
        end

    end
    last_col_x_IJAM(i) = find(any(x_hat_IJAM, 1), 1, 'last');
    last_col_a_IJAM(i) = find(any(a_hat_IJAM, 1), 1, 'last');
    
    [a_hat_IJAM(:,last_col_a_IJAM(i)) a_tilde]
    
    y_hat_unattacked = y-a_hat_IJAM(:,last_col_a_IJAM(i));

    x_hat_idea = pinv(C)*y_hat_unattacked;

    [[1:n]' x_tilde x_hat_idea x_hat_IJAM(:,last_col_x_IJAM(i))]
    return
    figure(i)  
    subplot(2,1,1);
    for j = 1:length(x_tilde)
        plot(x_hat_IJAM(j,1:last_col_x_IJAM(i)))
        hold on,grid on
        plot(x_tilde(j)*ones(1,last_col_x_IJAM(i)),'--k')
        title('x hat IJAM')
    end

    
    %ISTA ++++++++++++++++++++++++++++++++++++++++++++++++++
    nu = 0.99/(norm(G,2)^2);

    for k = 1:max_iteration
        
        x_hat_ISTA(:,k+1) = x_hat_ISTA(:,k) - ... 
                            nu*C'*(C*x_hat_ISTA(:,k) + a_hat_ISTA(:,k)- y);
      
        a_hat_ISTA(:,k+1) = prox_l1(a_hat_ISTA(:,k) - nu*(C*x_hat_ISTA(:,k) ...
                             +a_hat_ISTA(:,k)-y),lambda);
       
        if(norm(x_hat_ISTA(:,k+1)-x_hat_ISTA(:,k))^2<delta)  % condition for maximum iteration number
            break
        end
    end
    last_col_x_ISTA(i) = find(any(x_hat_ISTA, 1), 1, 'last'); %finding the last non-zero column
    last_col_a_ISTA(i) = find(any(a_hat_ISTA, 1), 1, 'last');


    figure(i)  
    subplot(2,1,2);
    for j = 1:length(x_tilde)
        plot(x_hat_ISTA(j,1:last_col_x_ISTA(i)))
        hold on,grid on
        plot(x_tilde(j)*ones(1,last_col_x_ISTA(i)),'--k')
        title('x hat ISTA')
    end

    % Calculating performance matrices ++++++++++++++++++
    state_estimation_error_IJAM(i) = norm(x_hat_IJAM(:,last_col_x_IJAM(i))-x_tilde)/norm(x_tilde);
    state_estimation_error_ISTA(i) = norm(x_hat_ISTA(:,last_col_x_ISTA(i))-x_tilde)/norm(x_tilde);
    
    attack_idx = find(a_tilde)
    attack_hat_IJAM_idx = find((a_hat_IJAM(:,last_col_a_IJAM(i))))
    attack_hat_ISTA_idx = find((a_hat_ISTA(:,last_col_a_ISTA(i)))) %%seems weired
    
    %IJAM SUPPORT ERROR
    for j = 1:h          
        if attack_idx(j) ~= attack_hat_IJAM_idx(j)
            attack_support_error_IJAM(i) =attack_support_error_IJAM(i) +1;  
        end                                      
    end
    
    %ISTA SUPPORT ERROR
    for j = 1:h
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


    
    % take into account last_col_IJAM and last_col_ISTA
end


% Making average of the performance matrices

%average convergance comparison
mean_last_col_x_IJAM = mean(last_col_x_IJAM)
mean_last_col_x_ISTA = mean(last_col_x_ISTA)

% average attack support 
mean_attack_support_error_IJAM = mean(attack_support_error_IJAM)
mean_attack_support_error_ISTA = mean(attack_support_error_ISTA)

% 
mean_state_estimation_error_IJAM = mean(state_estimation_error_IJAM)
mean_state_estimation_error_ISTA = mean(state_estimation_error_ISTA)



