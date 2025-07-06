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
title('Ring Network Topology Graph')

%%

max_iteration = 1e4;
delta = 1e-8;

lambda1 = 1*2;
lambda2 = 1*2;
nu1 = 0.015/norm(G(:,1:n),2)^2;
nu2 = 0.015/norm(G(:,n+1:end),2)^2;

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

last_index_ring = k+1;

% cleaning the data since we only have one target
x_hat_DISTA_star = zeros(n,q);
a_hat_DISTA_star = zeros(q,q);

% retrieving position
for node = 1:q
    %position
    [~,idx] = max(z_hat_DISTA(1:n,k+1,node));
    x_hat_DISTA_star(idx,node) = 1;
end

% ploting the heat map of the position
figure(2)
imagesc(x_hat_DISTA_star);                 % 100 x q matrix
colormap([0 0 0; 1 1 1]);                  % Black for 0, White for 1
colorbar('Ticks',[0,1], 'TickLabels',{'0','1'});

xlabel('Node Index');
ylabel('Cell Index');
title('state support heat map after the consensus - ring structure');

% Highlight only cell 37 on y-axis
yticks(37);
yticklabels({'37'});

% retrieving the attack vector
for node = 1:q
    idx = find(z_hat_DISTA(n+1:end,k+1,node)>0.07);
    a_hat_DISTA_star(idx,node) = 1;
end

% attack heat map
figure(3)
imagesc(a_hat_DISTA_star);                 % 100 x q matrix
colormap([0 0 0; 1 1 1]);                  % Black for 0, White for 1
colorbar('Ticks',[0,1], 'TickLabels',{'0','1'});

xlabel('Node Index');
ylabel('Cell Index');
title('attack support heat map after the consensus - ring structure');

% Highlight only cell 37 on y-axis
yticks(attack_tilde_support);
yticklabels({'1','10','14','16','17'});


%% star

load Q_star.mat
figure(3)
plot(digraph(Q))
title('Ring Network Topology Graph')

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

last_index_star = k+1;

% cleaning the data since we only have one target
x_hat_DISTA_star = zeros(n,q);
a_hat_DISTA_star = zeros(q,q);

% retrieving position
for node = 1:q
    %position
    [~,idx] = max(z_hat_DISTA(1:n,k+1,node));
    x_hat_DISTA_star(idx,node) = 1;
end

% ploting the heat map of the position
figure(4)
imagesc(x_hat_DISTA_star);                 % 100 x q matrix
colormap([0 0 0; 1 1 1]);                  % Black for 0, White for 1
colorbar('Ticks',[0,1], 'TickLabels',{'0','1'});

xlabel('Node Index');
ylabel('Cell Index');
title('state support heat map after the consensus - star structure');

% Highlight only cell 37 on y-axis
yticks(37);
yticklabels({'37'});

% retrieving the attack vector
for node = 1:q
    idx = find(z_hat_DISTA(n+1:end,k+1,node)>0.07);
    a_hat_DISTA_star(idx,node) = 1;
end

% attack heat map
figure(5)
imagesc(a_hat_DISTA_star);                 % 100 x q matrix
colormap([0 0 0; 1 1 1]);                  % Black for 0, White for 1
colorbar('Ticks',[0,1], 'TickLabels',{'0','1'});

xlabel('Node Index');
ylabel('Cell Index');
title('attack support heat map after the consensus - star structure');

% Highlight only cell 37 on y-axis
yticks(attack_tilde_support);
yticklabels({'1','10','14','16','17'});
%%

figure(6)
bar([1,2],[last_index_ring, last_index_star])
xticks([1 2]);
xticklabels({'Ring','Star'})
yticks([last_index_ring, last_index_star]);
ylabel('number of iterations')
title('Number of iterations to reach consensus')

