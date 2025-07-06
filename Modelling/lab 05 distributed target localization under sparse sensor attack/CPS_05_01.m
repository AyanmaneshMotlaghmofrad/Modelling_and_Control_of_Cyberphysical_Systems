clc
clear variables
close all
format compact
%%
load localization_data.mat

n = size(D,2)
q = size(D,1)

D_non_normalized = D;
G = normalize([D_non_normalized eye(q)]);
D = G(:,1:n);

%% ring

load Q_ring.mat
figure(1)
plot(digraph(Q))


max_iteration = 1e4;
delta = 1e-8;

lambda1 = 1;
lambda2 = 1;
nu1 = 0.01/norm(G(:,1:n),2)^2;
nu2 = 0.01/norm(G(:,n+1:end),2)^2;

z_hat_DISTA = zeros(n+q,max_iteration,q);
attack_tilde_support = [1,10,14,16,17];


for k = 1:max_iteration   
    for node = 1:q
        % calculating the average term
        average_term = zeros(n+q,1);

        for neighbor = find(Q(node,:))
            average_term = average_term +Q(node,neighbor)*z_hat_DISTA(:,k,neighbor);
        end
        
        % updating the state
        z_hat_DISTA(:,k+1,node) = prox_l1(average_term + nu1*G(node,:)'*(y(node)-G(node,:)*z_hat_DISTA(:,k,node)) ,nu1*lambda1);
    end

    if(norm(z_hat_DISTA(:,k+1,node)-z_hat_DISTA(:,k,node))^2<delta)  % condition for maximum iteration number
        break
    end
end

last_index = k+1;


for node = 1:q
    %position
    x_hat_DISTA(:,node) = z_hat_DISTA(1:n,k+1,node);
    node
    [~,idx] = max(x_hat_DISTA(:,node))

end

for node = 1:q
    %position
    a_hat_DISTA(:,node) = z_hat_DISTA(n+1:end,k+1,node);
    node
    [~,idx] = maxk(a_hat_DISTA(:,node),5);
    sort(idx)
end
% for node = 1:q
%     x_hat_DISTA(:,node) = z_hat_DISTA(1:n,k+1,node);
%     a_hat_DISTA(:,node) = z_hat_DISTA(n+1:end,last_index,node);
% 
%     non_zero_x = find(x_hat_DISTA);
%     non_zero_a = find(a_hat_DISTA);
% 
%     cleaning_threshold = 0.29;
% 
%     % the state of the cell can be {0,1}
%     % so we have to clean the data
%     for i = 1:length(non_zero_x)
%         if(x_hat_DISTA(non_zero_x(i))<cleaning_threshold)
%             x_hat_DISTA(non_zero_x(i)) = 0;
% 
%         else
%             x_hat_DISTA(non_zero_x(i)) = 1;
%         end
% 
%     end
%     node
%     target_position = find(x_hat_DISTA)
% 
% 
% end
% 

%% star

load Q_star.mat
figure(2)
plot(digraph(Q))
%% STAR configuration

max_iteration = 1e4;
delta = 1e-8;

lambda1 = 1;
lambda2 = 1;
nu1 = 0.01/norm(G(:,1:n),2)^2;
nu2 = 0.01/norm(G(:,n+1:end),2)^2;

z_hat_DISTA = zeros(n+q,max_iteration,q);
attack_tilde_support = [1,10,14,16,17];


for k = 1:max_iteration   
     % calculating the average term
    average_term0 = zeros(n+q,1);

    for neighbor = find(Q(1,:))
        average_term0 = average_term0 + Q(1,neighbor)*z_hat_DISTA(:,k,neighbor-1);
    end

    for node = 1:q
        % updating the state
        z_hat_DISTA(:,k+1,node) = prox_l1(average_term0 + nu1*G(node,:)'*(y(node)-G(node,:)*z_hat_DISTA(:,k,node)) ,nu1*lambda1);
    end

    if(norm(z_hat_DISTA(:,k+1,node)-z_hat_DISTA(:,k,node))^2<delta)  % condition for maximum iteration number
        break
    end
end

last_index = k+1;


for node = 1:q
    %position
    x_hat_DISTA(:,node) = z_hat_DISTA(1:n,k+1,node);
    node
    [~,idx] = max(x_hat_DISTA(:,node))

end

for node = 1:q
    %position
    a_hat_DISTA(:,node) = z_hat_DISTA(n+1:end,k+1,node);
    node
    [~,idx] = maxk(a_hat_DISTA(:,node),5);
    sort(idx)
end
% for node = 1:q
%     x_hat_DISTA(:,node) = z_hat_DISTA(1:n,k+1,node);
%     a_hat_DISTA(:,node) = z_hat_DISTA(n+1:end,last_index,node);
% 
%     non_zero_x = find(x_hat_DISTA);
%     non_zero_a = find(a_hat_DISTA);
% 
%     cleaning_threshold = 0.29;
% 
%     % the state of the cell can be {0,1}
%     % so we have to clean the data
%     for i = 1:length(non_zero_x)
%         if(x_hat_DISTA(non_zero_x(i))<cleaning_threshold)
%             x_hat_DISTA(non_zero_x(i)) = 0;
% 
%         else
%             x_hat_DISTA(non_zero_x(i)) = 1;
%         end
% 
%     end
%     node
%     target_position = find(x_hat_DISTA)
% 
% 
% end
% 
