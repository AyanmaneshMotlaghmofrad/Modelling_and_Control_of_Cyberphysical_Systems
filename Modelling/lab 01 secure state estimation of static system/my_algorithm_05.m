clc
clear variables 
close all
format compact

h = 1; % resilient up to 8 attack
%% System dimension information

n = 15; % the number of the states
q = 30; % the length of the output vector
run = 1000;

count = 0;

State_estimation_error = zeros(1,run);

for i = 1:run
    C = 0.01*randn(q, n);
    % C = 10*eye(q,n)
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
    
    % attack detection+++++++++++++++++++++++++++++++++++++++++++++++++++
    G = [C,eye(q)];
    
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
    
    if all(sum((a_sparse ~= 0) ~= (a_tilde ~= 0)) == 0)
        count = count + 1;
    end
    %          ----------------------------------
    
    % 
    % identification+++++++++++++++++++++++++++++++++++++++++++++++++++++
    % shaping Ga for G_sparse = (C Ga)
    Ga = zeros(q);
    
    for j = 1:length(non_zero_idx_sparse)
        Ga(non_zero_idx_sparse(j),non_zero_idx_sparse(j)) = 1;
    end
    
    
    G_sparse = [C, Ga];
    z_sparse = pinv(G_sparse)*y;

    x_hat = z_sparse(1:n);
    a_hat = z_sparse(n+1:end);


    
    % estimation performance
    state_estimation_error(i) = norm(x_hat-x_tilde)/norm(x_tilde);
    
    % attack detection


    % non_zero_idx_sparse
    % [[1:n]' x_hat x_tilde]
    % [[1:q]' a_non_sparse a_sparse a_hat a_tilde]
    count %if count =1 the detection of the attack is exact
          %it does not detect fully in half of the times

    state_estimation_error(i) % in case of correct detection
                              % one order of magnitude more exact than IJAM
    
end

% exact detection percentage
detection_performance(h) =count/run*100
mean_state_estimation_error(h) = mean(state_estimation_error)
h = h +1;


%% plotting
figure(1)
bar(1:length(detection_performance),detection_performance,'g')
xlabel('h [the number of attacks]')
ylabel('detection_precision %')
title('Detection versus the number of attacks')


figure(2)

bar(1:length(mean_state_estimation_error),mean_state_estimation_error)
xlabel('number of attacks') 
ylabel('mean estimation error') 
title('Mean estimation error Vs. number of attacks')
grid on
