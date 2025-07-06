clc
clear variables 
close all
format compact

%% System dimension information

n = 15; % the number of the states
q = 30; % the length of the output vector
h = 2; % two sensors are under attack.
run = 1000;

count = 0;

for i = 1:run
    C = randn(q,n);
    G = [C,eye(q)];
    % support of the attack vector a: uniform distribution
    % and "real" attack +++++++++++++++++++++++++++++++++++++++++++++++++
    
    a_tilde = zeros(q,1);
    h_count = h;
    
    while 1
        
        idx = randi([1,15]);
    
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
    
    % attack detection
    
    
    z = pinv(G)*y;
    
    % sparsifying
    a_non_sparse = z(n+1:end);
    
    a_sparse = zeros(q,1);
    
    [~, idx] = maxk(abs(a_non_sparse), h); % Find indices of the h largest values
    
    a_sparse(idx) = a_non_sparse(idx);
    
    [a_tilde a_sparse a_non_sparse];
    
    % checking the performance ++++++++++++++++++++++++++++++++++++++++++
    non_zero_idx_tilde = find(a_tilde);
    
    non_zero_idx_sparse = find(a_sparse);  
    

    if all(ismember(non_zero_idx_tilde, non_zero_idx_tilde)) && all(ismember(non_zero_idx_sparse, non_zero_idx_tilde))
        count = count + 1;
    end

    % identification+++++++++++++++++++++++++++++++++++++++++++++++++++++
    % first method (does not work)
    
    y_unattacked = y;
    y_unattacked(non_zero_idx_sparse) = [];

    y_attacked = y(non_zero_idx_sparse);
    

    C_unattacked = C;
    C_unattacked(:,non_zero_idx_sparse) = [];
    C_unattacked(non_zero_idx_sparse,:) = [];
    C_attacked = C(non_zero_idx_sparse,:);
    
    % identification of the unattacked states
    x_hat_unattacked = pinv(C_unattacked)*y_unattacked;   
    
    % identificaiton of the attacked states
    G_attacked = [C_attacked eye(h)];
    z_hat_attacked = pinv(G_attacked)*y_attacked;
    
    x_hat_attacked = z_hat_attacked(1:h);
    a_hat_attacked = z_hat_attacked(h+1:2*h);

    %shaping the estimated vectors
    x_hat_1 = zeros(n,1);
    idx_attack = 1; %for counting
    idx_unattack = 1;

    for idx_hat = 1:n %maybe does not work when the first element is attacked       
        if(non_zero_idx_sparse(idx_attack) ~=idx_hat)
            x_hat_1(idx_hat) = x_hat_unattacked(idx_unattack);
            idx_unattack = idx_unattack+1;

        else
            x_hat_1(idx_hat) = x_hat_attacked(idx_attack);
            idx_attack = idx_attack + 1;
            if(idx_attack > h)
                break
            end
        end
    end
 
    % second method
    % it seems like it works well, at least for x_hat
    Ga = zeros(q);
    
    for i = non_zero_idx_sparse
        Ga(i,i) = 1;
    end
      
    
    G_sparse = [C, Ga];
    z_sparse = pinv(G_sparse)*y;

    x_hat_2 = z_sparse(1:n);
    a_hat = z_sparse(n+1:end);

    % non_zero_idx_sparse
    % [[1:n]' x_hat_1 x_hat_2 x_tilde]
    % [a_hat a_sparse a_tilde]
    % 
    detection_performance =count/run*100 %delete this line with return
    
 
end

% detection performance
detection_performance =count/run*100


