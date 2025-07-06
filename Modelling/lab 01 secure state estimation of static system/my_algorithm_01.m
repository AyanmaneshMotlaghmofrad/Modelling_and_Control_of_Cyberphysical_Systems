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
    
    % attack detection+++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    z = pinv(G)*y;
    
    % sparsifying
    a_non_sparse = z(n+1:end);
    
    a_sparse = zeros(q,1);
    
    [~, idx] = maxk(abs(a_non_sparse), h); % Find indices of the h largest values
    
    a_sparse(idx) = a_non_sparse(idx);
    
    [a_tilde a_sparse a_non_sparse]
    
    
    % checking detection performance +++++++++++++++++++++++++++++++
    non_zero_idx_tilde = find(a_tilde);
    
    non_zero_idx_sparse = find(a_sparse);
    
    
    if all(ismember(non_zero_idx_tilde, non_zero_idx_tilde)) && all(ismember(non_zero_idx_sparse, non_zero_idx_tilde))
        count = count + 1;
    end
    %          ----------------------------------
    
    % identification+++++++++++++++++++++++++++++++++++++++++++++++++++++
    % shaping Ga for G_sparse = (C Ga)
    Ga = zeros(q);
    
    for i = 1:length(non_zero_idx_sparse)
        Ga(non_zero_idx_sparse(i),non_zero_idx_sparse(i)) = 1;
    end
    
    
    G_sparse = [C, Ga];
    z_sparse = pinv(G_sparse)*y;

    x_hat = z_sparse(1:n);
    a_hat = z_sparse(n+1:end);

    non_zero_idx_sparse
    [[1:15]' x_hat x_tilde]
    [a_sparse a_hat a_tilde]
    
    detection_performance =count/run*100 %delete this line with return
    
 
end

% detection performance
detection_performance =count/run*100


