clc
clear variables
close all
format compact

% 1) design a local controller to stabilize the system
% since the system is unstable but since we need also the observer
% once we need to do it with local and once distributed observer.

% 2) test with noise

%% parameters of the system and the leader referece
n_sys = 2; % the number of the states
N = 6; % the number of the followers

% unstable system
A_unstable = [0       1
     880.87  0];
B_unstable = [0 -9.9453]';
C_unstable = [708.27 0];
D_unstable = 0;

rank(obsv(A_unstable,C_unstable));
rank(ctrb(A_unstable,B_unstable));

% local observer and controller
lambda_obs = [-10 -15];
L = place(A_unstable',C_unstable',lambda_obs)';

lambda_ctrl = [0 -5]; 

s= tf('s');
K_loc = acker(A_unstable,B_unstable,lambda_ctrl);
N_loc = 1/dcgain(s*tf(ss(A_unstable-B_unstable*K_loc,B_unstable,C_unstable,D_unstable)));
% N_loc = 0.0114*2.5/100;

A = [A_unstable -B_unstable*K_loc;L*C_unstable A_unstable-B_unstable*K_loc-L*C_unstable];
B = [B_unstable;B_unstable]*N_loc;
C = [C_unstable zeros(1,n_sys)];
D = 0;

%dynamics after the controller
sys0 = ss(A,B,C,D);
A_node = A_unstable-B_unstable*K_loc;
B_node = B_unstable;
C_node = C_unstable;
D_node = D_unstable;
sys_node = ss(A_unstable-B_unstable*K_loc,B_unstable*N_loc,C_unstable,D_unstable);

t = linspace(0,80,5000);
u = zeros(1,length(t));
Ts = t(2) - t(1);
u(1) = 0.025/Ts;

% f = 0.02; %frequency of square wave, in Hz
% u = square(2*pi*f*t);


x0 = [0 0 0 0]';



[y0_transpose,~,x0_transpose] = lsim(sys0,u,t,x0);


plot(t,y0_transpose)
ylabel('x(1) [m]')
title('sinusoidal reference generated by the leader')
xlabel('time [sec]')

x0_bar = repmat(x0_transpose',[N 1]);

%% First Network topology 
% tree topology
close all
n = 4;
G_cal = diag([1 0 0 0 0 0]);
A_cal = [0 0 0 0 0 0;
         1 0 0 0 0 0;
         0 1 0 0 0 0;
         0 0 1 0 0 0;
         0 0 0 1 0 0;
         0 0 0 0 1 0];
% later use different A_cal
% for the controller and the observer

D_cal= diag(sum(A_cal,2));
L_cal = D_cal-A_cal;
lambda = eig(L_cal+G_cal);


% Distributed control protocol
Q_c = 1e4*eye(4);
R_c = 1;

P_c = are(A,B*inv(R_c)*B',Q_c);
K = inv(R_c)*B'*P_c;

c_c = zeros(N); %coupling gain for the controller

for i = 1:6
    c_c(i,i) = (real(lambda(i)))/2+1;
end

Ac = kron(eye(N),A) - kron(c_c*(L_cal+G_cal),B*K);
Bc = kron(c_c*(L_cal+G_cal),B*K);
Cc = eye(N*n);
Dc = zeros(N*n);



sys_global = ss(Ac,Bc,Cc,Dc);
[yr,~,xr]= lsim(sys_global, x0_bar, t, zeros(4*N,1)); 

for i = 1:2
    if(i == 1 || i ==3)
         figure(i)
         plot(t,xr(:,i:4:end)*708.2700)
         hold on
         plot(t,y0_transpose)
         legend({'Node 01','Node 02','Node 03','Node 04','Node 05','Node 06','Node 00'})
         title('node sinusoidal responses in tree topology with local observer')
         ylabel('x(1) [m]')
         xlabel('time [sec]')
    else
        figure(i)
         plot(t,xr(:,i:4:end))
         legend
    end
end

% plotting the graph
% Number of agents and leader
n_agents = size(A_cal, 1);
n_total = n_agents + 1;

A_aug =[
     0     0     0     0     0     0     1
     1     0     0     0     0     0     0
     0     1     0     0     0     0     0
     0     0     1     0     0     0     0
     0     0     0     1     0     0     0
     0     0     0     0     1     0     0
     0     0     0     0     0     0     0]
       
     
% Create directed graph
G = digraph(A_aug');

% Color: blue = agents, green = leader
node_colors = repmat([0 0 1], n_total, 1); % blue
node_colors(end, :) = [0 0.6 0];           % green for leader

% Labels
labels = arrayfun(@(x) sprintf('Node %d', x), 1:n_agents, 'UniformOutput', false);
labels{end+1} = 'Leader';

% Plot
figure
plot(G, ...
     'NodeColor', node_colors, ...
     'MarkerSize', 7, ...
     'LineWidth', 1.5, ...
     'NodeLabel', labels);
title('Directed Network: Tree Topology')

%% first structure with NOISE

% xdot = Ax + Bu+ Gw, w:n*q q = 2 source of noise
% y = Cx + Du + v  

q = 12; % two source of noise for each node;
m = 1; % measurement noise

G = 0.00005*randn(24,12); %process noise gain G:n*q
w = 0.0001*randn(q,length(t));
v = 0.001*randn(24,length(t));
% v = zeros(24,length(t));

G_cal = diag([1 0 0 0 0 0]);
A_cal = [0 0 0 0 0 0;
         1 0 0 0 0 0;
         0 1 0 0 0 0;
         0 0 1 0 0 0;
         0 0 0 1 0 0;
         0 0 0 0 1 0]*1;
% later use different A_cal
% for the controller and the observer

D_cal= diag(sum(A_cal,2));
L_cal = D_cal-A_cal;
lambda = eig(L_cal+G_cal);


% Distributed control protocol
Q_c = 1e4*eye(4);
R_c = 1;

P_c = are(A,B*inv(R_c)*B',Q_c);
K = inv(R_c)*B'*P_c;

c_c = zeros(N); %coupling gain for the controller

for i = 1:6
    c_c(i,i) = (real(lambda(i)))/2+100;
end

Ac = kron(eye(N),A) - kron(c_c*(L_cal+G_cal),B*K);
Bc = kron(c_c*(L_cal+G_cal),B*K);
Cc = eye(N*n);
Dc = zeros(N*n);


sys_global = ss(Ac,Bc,Cc,Dc);
xr = zeros(24,length(t));
yr = zeros(24,length(t));

for k = 1:length(t)-1
    dx = Ac * xr(:,k) + Bc * x0_bar(:,k);       % deterministic dynamics
    xr(:,k+1) = xr(:,k) + Ts * dx + G * w(:,k);  % Euler integration + process noise
    yr(:,k+1) = Cc * xr(:,k+1) + Dc * x0_bar(:,k) + v(:,k);  % output + measurement noise
end


for i = 1:2
    if(i == 1 || i ==3)
         figure(i)
         plot(t,xr(i:4:end,:)*708.2700)
         hold on
         plot(t,y0_transpose)
         legend({'Node 01','Node 02','Node 03','Node 04','Node 05','Node 06','Node 00'})
         title('node noisy step responses in tree topology with local observer')
         ylabel('x(1) [m]')
         xlabel('time [sec]')
    else
        figure(i)
         plot(t,xr(i:4:end,:))
         legend
    end
end

%% second Network topology 
% ring topology
n = 4;
G_cal = diag([1 0 0 0 0 0]);
A_cal = [0 1 0 0 0 1;
         1 0 1 0 0 0;
         0 1 0 1 0 0;
         0 0 1 0 1 0;
         0 0 0 1 0 1;
         1 0 0 0 1 0];

% later use different A_cal
% for the controller and the observer

D_cal= diag(sum(A_cal,2));
L_cal = D_cal-A_cal;
lambda = eig(L_cal+G_cal);


% Distributed control protocol
Q_c = 1e4*eye(4);
R_c = 1;

P_c = are(A,B*inv(R_c)*B',Q_c);
K = inv(R_c)*B'*P_c;

c_c = zeros(N); %coupling gain for the controller

for i = 1:6
    c_c(i,i) = (real(lambda(i)))/2+1;
end

Ac = kron(eye(N),A) - kron(c_c*(L_cal+G_cal),B*K);
Bc = kron(c_c*(L_cal+G_cal),B*K);
% C_single = [
%     708.27     0     0       0  ;
%       0        0   708.27    0  ;
%       0        0     0       1  ;
% ];
% Cc = kron(eye(N), C_single);
% Dc = zeros(N*3,N*4);

Cc = eye(N*n);
Dc = zeros(N*n);



sys_global = ss(Ac,Bc,Cc,Dc);
[yr,~,xr]= lsim(sys_global, x0_bar, t, zeros(4*N,1)); 



for i = 3:4
    if(i-2 == 1 || i-2 ==3)
         figure(i)
         plot(t,xr(:,i-2:4:end)*708.2700)
         hold on
         plot(t,y0_transpose)
         legend({'Node 01','Node 02','Node 03','Node 04','Node 05','Node 06','Node 00'})
         title('node sinusoidal responses in ring topology with local observer')
         ylabel('x(1) [m]')
         xlabel('time [sec]')
    else
        figure(i)
         plot(t,xr(:,i-2:4:end))
         legend
    end
end

% plotting the graph
% Number of agents and leader
n_agents = size(A_cal, 1);
n_total = n_agents + 1;

A_aug =[
     0     1     0     0     0     1     1
     1     0     1     0     0     0     0
     0     1     0     1     0     0     0
     0     0     1     0     1     0     0
     0     0     0     1     0     1     0
     1     0     0     0     1     0     0
     0     0     0     0     0     0     0]

% Create directed graph
G = digraph(A_aug');

% Color: blue = agents, green = leader
node_colors = repmat([0 0 1], n_total, 1); % blue
node_colors(end, :) = [0 0.6 0];           % green for leader

% Labels
labels = arrayfun(@(x) sprintf('Node %d', x), 1:n_agents, 'UniformOutput', false);
labels{end+1} = 'Leader';

% Plot
figure(1)
plot(G, ...
     'NodeColor', node_colors, ...
     'MarkerSize', 7, ...
     'LineWidth', 1.5, ...
     'NodeLabel', labels);
title('Directed Network: Ring Topology')

%% Forth Network topology 
% fully connected-nodes topology
n = 4;
G_cal = diag([1 0 0 0 0 0]);
A_cal = ones(6)- diag(ones(6,1));

% later use different A_cal
% for the controller and the observer
D_cal= diag(sum(A_cal,2));
L_cal = D_cal-A_cal;
lambda = eig(L_cal+G_cal);


% Distributed control protocol
Q_c = 1e4*eye(4);
R_c = 1;

P_c = are(A,B*inv(R_c)*B',Q_c);
K = inv(R_c)*B'*P_c;

c_c = zeros(N); %coupling gain for the controller

for i = 1:6
    c_c(i,i) = (real(lambda(i)))/2+1;
end

Ac = kron(eye(N),A) - kron(c_c*(L_cal+G_cal),B*K);
Bc = kron(c_c*(L_cal+G_cal),B*K);
% C_single = [
%     708.27     0     0       0  ;
%       0        0   708.27    0  ;
%       0        0     0       1  ;
% ];
% Cc = kron(eye(N), C_single);
% Dc = zeros(N*3,N*4);

Cc = eye(N*n);
Dc = zeros(N*n);



sys_global = ss(Ac,Bc,Cc,Dc);
[yr,~,xr]= lsim(sys_global, x0_bar, t, zeros(4*N,1)); 



for i = 7:8
    if(i-6 == 1 || i-6 ==3)
         figure(i)
         plot(t,xr(:,i-6:4:end)*708.2700)
         hold on
         plot(t,y0_transpose)
         legend({'Node 01','Node 02','Node 03','Node 04','Node 05','Node 06','Node 00'})
         title('node sinusoidal responses in full topology with local observer')
         ylabel('x(1) [m]')
         xlabel('time [sec]')   
    else
         figure(i)
         plot(t,xr(:,i-6:4:end))
         legend
    end
end

%% fully connected
n = 4;
G_cal = diag([1 0 0 0 0 0]);
A_cal = ones(6)- diag(ones(1,6));

% later use different A_cal
% for the controller and the observer

D_cal= diag(sum(A_cal,2));
L_cal = D_cal-A_cal;
lambda = eig(L_cal+G_cal);


% Distributed control protocol
Q_c = 1e4*eye(4);
R_c = 1;

P_c = are(A,B*inv(R_c)*B',Q_c);
K = inv(R_c)*B'*P_c;

c_c = zeros(N); %coupling gain for the controller

for i = 1:6
    c_c(i,i) = (real(lambda(i)))/2+1;
end

Ac = kron(eye(N),A) - kron(c_c*(L_cal+G_cal),B*K);
Bc = kron(c_c*(L_cal+G_cal),B*K);
% C_single = [
%     708.27     0     0       0  ;
%       0        0   708.27    0  ;
%       0        0     0       1  ;
% ];
% Cc = kron(eye(N), C_single);
% Dc = zeros(N*3,N*4);

Cc = eye(N*n);
Dc = zeros(N*n);



sys_global = ss(Ac,Bc,Cc,Dc);
[yr,~,xr]= lsim(sys_global, x0_bar, t, zeros(4*N,1)); 



for i = 7:8
    if(i-6 == 1 || i-6 ==3)
         figure(i)
         plot(t,xr(:,i-6:4:end)*708.2700)
         legend
    else
        figure(i)
         plot(t,xr(:,i-6:4:end))
         legend
    end
end


% plotting the graph
% Number of agents and leader
n_agents = size(A_cal, 1);
n_total = n_agents + 1;

A_aug =[
     0     1     1     1     1     1     1
     1     0     1     1     1     1     0
     1     1     0     1     1     1     0
     1     1     1     0     1     1     0
     1     1     1     1     0     1     0
     1     1     1     1     1     0     0
     0     0     0     0     0     0     0]

% Create directed graph
G = digraph(A_aug');

% Color: blue = agents, green = leader
node_colors = repmat([0 0 1], n_total, 1); % blue
node_colors(end, :) = [0 0.6 0];           % green for leader

% Labels
labels = arrayfun(@(x) sprintf('Node %d', x), 1:n_agents, 'UniformOutput', false);
labels{end+1} = 'Leader';

% Plot
figure(1)
plot(G, ...
     'NodeColor', node_colors, ...
     'MarkerSize', 7, ...
     'LineWidth', 1.5, ...
     'NodeLabel', labels);
title('Directed Network: Fully Connected Topology')
%% Fifth Network topology 
% star nodes topology
n = 4;
G_cal = diag([1 0 0 0 0 0]);
A_cal = [0 1 1 1 1 1;
         1 0 0 0 0 0;
         1 0 0 0 0 0;
         1 0 0 0 0 0;
         1 0 0 0 0 0;
         1 0 0 0 0 0];
% later use different A_cal
% for the controller and the observer

D_cal= diag(sum(A_cal,2));
L_cal = D_cal-A_cal;
lambda = eig(L_cal+G_cal);


% Distributed control protocol
Q_c = 1e4*eye(4);
R_c = 1;

P_c = are(A,B*inv(R_c)*B',Q_c);
K = inv(R_c)*B'*P_c;

c_c = zeros(N); %coupling gain for the controller

for i = 1:6
    c_c(i,i) = (real(lambda(i)))/2+1;
end

Ac = kron(eye(N),A) - kron(c_c*(L_cal+G_cal),B*K);
Bc = kron(c_c*(L_cal+G_cal),B*K);
% C_single = [
%     708.27     0     0       0  ;
%       0        0   708.27    0  ;
%       0        0     0       1  ;
% ];
% Cc = kron(eye(N), C_single);
% Dc = zeros(N*3,N*4);

Cc = eye(N*n);
Dc = zeros(N*n);



sys_global = ss(Ac,Bc,Cc,Dc);
[yr,~,xr]= lsim(sys_global, x0_bar, t, zeros(4*N,1)); 



for i = 9:10
    if(i-8 == 1 || i-8 ==3)
         figure(i)
         plot(t,xr(:,i-8:4:end)*708.2700)
         hold on
         plot(t,y0_transpose)
         legend({'Node 01','Node 02','Node 03','Node 04','Node 05','Node 06','Node 00'})
         title('node sinusoidal responses in star topology with local observer')
         ylabel('x(1) [m]')
         xlabel('time [sec]') 
    else
        figure(i)
         plot(t,xr(:,i-8:4:end))
         legend
    end
end


% plotting the graph
% Number of agents and leader
n_agents = size(A_cal, 1);
n_total = n_agents + 1;

A_aug =[
     0     1     1     1     1     1     1
     1     0     0     0     0     0     0
     1     0     0     0     0     0     0
     1     0     0     0     0     0     0
     1     0     0     0     0     0     0
     1     0     0     0     0     0     0
     0     0     0     0     0     0     0]

% Create directed graph
G = digraph(A_aug');

% Color: blue = agents, green = leader
node_colors = repmat([0 0 1], n_total, 1); % blue
node_colors(end, :) = [0 0.6 0];           % green for leader

% Labels
labels = arrayfun(@(x) sprintf('Node %d', x), 1:n_agents, 'UniformOutput', false);
labels{end+1} = 'Leader';

% Plot
figure(1)
plot(G, ...
     'NodeColor', node_colors, ...
     'MarkerSize', 7, ...
     'LineWidth', 1.5, ...
     'NodeLabel', labels);
title('Directed Network: Star Topology')

%% sixth Network topology 
% unidirected chain topology
n = 4;
A_cal = [0 1 0 0 0 0;
         1 0 1 0 0 0;
         0 1 0 1 0 0;
         0 0 1 0 1 0;
         0 0 0 1 0 1;
         0 0 0 0 1 0];
G_cal = diag([1 0 0 1 0 0]);
% later use different A_cal
% for the controller and the observer

D_cal= diag(sum(A_cal,2));
L_cal = D_cal-A_cal;
lambda = eig(L_cal+G_cal);


% Distributed control protocol
Q_c = 1e4*eye(4);
R_c = 1;

P_c = are(A,B*inv(R_c)*B',Q_c);
K = inv(R_c)*B'*P_c;

c_c = zeros(N); %coupling gain for the controller

for i = 1:6
    c_c(i,i) = (real(lambda(i)))/2+1;
end

Ac = kron(eye(N),A) - kron(c_c*(L_cal+G_cal),B*K);
Bc = kron(c_c*(L_cal+G_cal),B*K);
% C_single = [
%     708.27     0     0       0  ;
%       0        0   708.27    0  ;
%       0        0     0       1  ;
% ];
% Cc = kron(eye(N), C_single);
% Dc = zeros(N*3,N*4);

Cc = eye(N*n);
Dc = zeros(N*n);



sys_global = ss(Ac,Bc,Cc,Dc);
[yr,~,xr]= lsim(sys_global, x0_bar, t, zeros(4*N,1)); 



for i = 11:12
    if(i-10 == 1 || i-10 ==3)
         figure(i)
         plot(t,xr(:,i-10:4:end)*708.2700)
         legend
    else
        figure(i)
         plot(t,xr(:,i-10:4:end))
         legend
    end
end

%% Considering global observer
% parameters of the system and the leader referece
n_sys = 2; % the number of the states
N = 6; % the number of the followers

% local observer and controller
lambda_obs = [-10 -15];
L = place(A_unstable',C_unstable',lambda_obs)';

lambda_ctrl = [1*j -1*j]; 

s= tf('s');
K_loc = acker(A_unstable,B_unstable,lambda_ctrl);
% N_loc1 = 1/dcgain(tf(ss(A_unstable-B_unstable*K_loc,B_unstable,C_unstable,D_unstable)));
N_loc = 0.0114*2.5/100;

A = [A_unstable -B_unstable*K_loc;L*C_unstable A_unstable-B_unstable*K_loc-L*C_unstable];
B = [B_unstable;B_unstable]*N_loc;
C = [C_unstable zeros(1,n_sys)];
D = 0;

%dynamics after the controller
sys0 = ss(A,B,C,D);

t = linspace(0,7,5000);
u = zeros(1,length(t));
Ts = t(2) - t(1);
u(1) = 0.0125/Ts;


[y0_transpose,~,x0_transpose] = lsim(sys0,u,t,x0);
plot(t,y0_transpose)

x0_bar = repmat(x0_transpose(:,1:2)',[N 1]);

%% First Network topology global observer
% Chain topology
close all
n = 2;

%tree
G_cal = diag([1 0 0 0 0 0]);
A_cal = [0 0 0 0 0 0;
         1 0 0 0 0 0;
         0 1 0 0 0 0;
         0 0 1 0 0 0;
         0 0 0 1 0 0;
         0 0 0 0 1 0]; 

% %ring
% G_cal = diag([1 0 0 0 0 0]);
% A_cal = [0 1 0 0 0 1;
%          1 0 1 0 0 0;
%          0 1 0 1 0 0;
%          0 0 1 0 1 0;
%          0 0 0 1 0 1;
%          1 0 0 0 1 0];
% 
% %full
% G_cal = diag([1 0 0 0 0 0]);
% A_cal = ones(6)- diag(ones(1,6));

% % star
% G_cal = diag([1 0 0 0 0 0]);
% A_cal = [0 1 1 1 1 1;
%          1 0 0 0 0 0;
%          1 0 0 0 0 0;
%          1 0 0 0 0 0;
%          1 0 0 0 0 0;
%          1 0 0 0 0 0];

D_cal= diag(sum(A_cal,2));
L_cal = D_cal-A_cal;
lambda = eig(L_cal+G_cal);


% Distributed control protocol
Q_c = eye(n);
R_c = 1;

P_c = are(A_unstable,B_unstable*inv(R_c)*B_unstable',Q_c);
K = inv(R_c)*B_unstable'*P_c;

c_c = zeros(N); %coupling gain for the controller

for i = 1:6
    c_c(i,i) = (max(real(lambda)))/2+1;
end

% global observer
Q_o = eye(n); % two states
R_o = 1;

P_o = are(A_unstable',C_unstable'*inv(R_o)*C_unstable,Q_o);
F = P_o*C_unstable'*inv(R_o);

c_o = zeros(N); %coupling gain for the controller

for i = 1:6
    c_o(i,i) = (max(real(lambda)))/2 + 1;
end


% global matrices
% the first time run without the correction for the dcgain
Ao = kron(eye(N),A_node)-c_o(1,1)*kron((L_cal+G_cal),F*C_node);
Ac = [kron(eye(N),A_node), -c_c(1,1)*kron((L_cal+G_cal),B_node*K);
      c_o(1,1)*kron((L_cal+G_cal),F*C_node), Ao-c_c(1,1)*kron((L_cal+G_cal),B_node*K)];
Bc = [kron(c_c*(L_cal +G_cal),B_node*K);c_c(1,1)*kron((L_cal +G_cal),B_node*K)];
Cc = eye(N*2*n);
Dc = zeros(N*2*n,N*n);

sys_global = ss(Ac,Bc,Cc,Dc);


sys_global = ss(Ac,Bc,Cc,Dc);
[yr,~,xr]= lsim(sys_global, x0_bar, t, zeros(4*N,1)); 

for i = 1:2
    if(i == 1)
         figure(i)
         plot(t,xr(:,i:2:24)*708.2700)
         hold on
         plot(t,y0_transpose)
         legend({'Node 01','Node 02','Node 03','Node 04','Node 05','Node 06','Node 00'})
         title('node sinusoidal responses in tree topology with distributed observer')
         ylabel('x(1) [m]')
         xlabel('time [sec]') 

    else
        figure(i)
         plot(t,xr(:,i:2:12))
         legend
    end
end


%% local algorithm simulation
% tree topology
G_cal = diag([1 0 0 0 0 0]);
A_cal = [0 0 0 0 0 0;
         1 0 0 0 0 0;
         0 1 0 0 0 0;
         0 0 1 0 0 0;
         0 0 0 1 0 0;
         0 0 0 0 1 0]; 

% later use different A_cal
% for the controller and the observer

D_cal= diag(sum(A_cal,2));
L_cal = D_cal-A_cal;
lambda = eig(L_cal+G_cal);
Ts = 0.5;

sys_d = c2d(sys_node,Ts,'zoh');
Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;
Dd = sys_d.D;


% discrete simulation of the leader

t = 0:Ts:100;
u = zeros(1,length(t));
u(1) = 0.025/Ts;
x0 = [0 0 0 0]';
[y,~,x0_d_transpose] = lsim(c2d(sys0,Ts,'zoh'),u,t,x0);
x0_d = x0_d_transpose';

% plot(t,y)


% initialization the states and the input signal
x_d = zeros(4,length(t),N); %[x;x_hat]
x_d(:,1,:) = 0.001*rand([4 1 6]);
u_d = zeros(1,length(t),N);


% Distributed control protocol
n = 2;
Q_c = [1 0
        0  1];
R_c = 1;

K = dlqr(Ad, Bd, Q_c, R_c);

c_c = zeros(N); %coupling gain for the controller

for i = 1:6
    c_c(i,i) = (max(real(lambda)))/2+100;
end

% global observer
Q_o = [1 0
        0  1]; % two states
R_o = 1;

F = dlqr(Ad', Cd', Q_o, R_o)';

c_o = zeros(N); %coupling gain for the controller

for i = 1:6
    c_o(i,i) = (max(real(lambda)))/2 + 1;
end

epsilon_hat = zeros(2,6);
xsi = zeros(1,6);




for k = 1:length(t)-1
    %dynamic update
      % 2 stands for the number of states of  the observer and 6 for the number of nodes
    
    %calculating epsilon_hat for the plant
    % error calculations for command input
    epsilon_hat = zeros(2,6);
    xsi = zeros(1,6);
    for i = 1:N  
        for j = 1:N
            epsilon_hat(:,i) = epsilon_hat(:,i) + A_cal(i,j)*(x_d(3:4,k,j) - x_d(3:4,k,i));
        end
        epsilon_hat(:,i) =  epsilon_hat(:,i) + G_cal(i,i)*(x0_d(3:4,k) - x_d(3:4,k,i));
        

    end
    
    % calculating epsilon for the observer ( action on y_tilde)
    for i = 1:N  
        for j = 1:N
            xsi(:,i) = xsi(:,i) + A_cal(i,j)*(Cd*(x_d(1:2,k,j)-x_d(3:4,k,j)) - Cd*(x_d(1:2,k,i)-x_d(3:4,k,i)));
        end
        xsi(:,i) =  xsi(:,i) + G_cal(i,i)*(Cd*(x0_d(1:2,k)-x0_d(3:4,k)) - Cd*(x_d(1:2,k,i)-x_d(3:4,k,i)));
    
    end

    for i = 1:N
        u_d(1,k,i) = c_c(i,i)*K*epsilon_hat(:,i);
        
        % actual state equation
        x_d(1:2,k+1,i) = Ad*x_d(1:2,k,i) + Bd*u_d(1,k,i);        
        x_d(3:4,k+1,i) = Ad*x_d(3:4,k,i) + Bd*u_d(1,k,i) - c_o(i,i)*F*xsi(:,i);

    end

end



% plotting
close all
for i = 1:2 % each state
    figure(i)
    hold on
    for j = 1:6 % for all the nodes
        if(i == 1 || i == 3)
             plot(t,x_d(i,:,j)*708.27)
             hold on
             plot(t,y)
             legend({'Node 01','Node 02','Node 03','Node 04','Node 05','Node 06','Node 00'})
             title('distirbuted observer tree topology with Ts = 0.5')
        else
            plot(t,x_d(i,:,j))
            hold on
            % plot(x_d(i+2,:,j))
            legend
        end
    end
    
end



%% not needed to add one integrator for a step reference

% local controller with integrator
% A_aug = [0 -C_unstable;zeros(2,1) A_unstable];
% B_aug = [0;B_unstable];
% C_aug = [0 C_unstable];
% D_aug = 0;
% 
% lambda_ctrl = [-30 -31 -32]; % [q;x]
% 
% K_loc = place(A_aug,B_aug,lambda_ctrl);
% K_q = K_loc(1);
% K_x = K_loc(2:end);
% 
% local observer to make the system controllable
% lambda_obs = [-35 -37];
% L = place(A_unstable',C_unstable',lambda_obs)';
% 
% state space x = [q x x_hat]
% A = [  0              -C_unstable        zeros(1, n_sys);
%         -B_unstable*K_q  A_unstable         -B_unstable*K_x;
%         -B_unstable*K_q  L*C_unstable       A_unstable-L*C_unstable-B_unstable*K_x ];
% 
% B = [ 1;
%          zeros(n_sys,1);
%          zeros(n_sys,1) ];
% 
% C = [ 0  C_unstable  zeros(1,2) ];   % output is still y = C*x
% D = 0;
% 
% sys_augmented = ss(A, B, C, D);
% 
% step(sys_augmented)












