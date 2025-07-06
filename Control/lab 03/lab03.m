clc
clear variables 
close all
format compact

%% problem setting
N = 6 % the number of nodes
k1 = 1.5;
k2 = 1;
m1 = 1.1;
m2 = 0.9;

A = [0 0 1 0
     0 0 0 1
     -(k1+k2)/m1 k2/m1 0 0
     k2/m2 -k2/m2 0 0];

B = [0 0 1/m1 0]';

C = eye(4);
D = zeros(4,1);

sys = ss(A,B,C,D);


% topology of the network
G_cal = diag([1 0 0 0 0 0]);
A_cal = [0 0 0 0 0 0;
         2 0 0 0 0 0;
         0 6 0 0 0 0;
         0 0 1 0 0 0;
         0 0 0 1 0 0;
         0 0 0 0 3 0]; 
D_cal= diag(sum(A_cal,2));

L_cal = D_cal-A_cal;

lambda = eig(L_cal+G_cal);


%% Gain matrix and controller design
Q = eye(4);
R = 1;

P = are(A,B*inv(R)*B',Q);

K = inv(R)*B'*P;

c = zeros(N);

for i = 1:6
    c(i,i) = (real(lambda(i)))/2 + 1
end

%% global closed-loop system dynamics
% simulation for a node for the leader
sys = ss(A,B,C,D)
t = linspace(0,20,1000);
u = zeros(1,length(t));
x0 = [10 10 0 0];
x0_transponse = lsim(sys,u,t,x0)


x0_bar = repmat(x0_transponse',[N 1]);

% global system matrices
Ac = kron(eye(N),A) - kron(c*(L_cal+G_cal),B*K)
Bc = kron(c*(L_cal+G_cal),B*K);
Cc = eye(N*4);
Dc = zeros(N*4);

sys_global = ss(Ac,Bc,Cc,Dc);
yr = lsim(sys_global, x0_bar, t, zeros(4*N,1)); 

for i = 1:4
    figure(i)
    plot(t,yr(:,i:4:end))
    hold on
    plot(t,x0_transponse(:,i))
end

%% designing a local controller for each node so that the leader can
% generate a step reference

rank(ctrb(A,B))
% the system is fully controllable
lambda_des = [0,-5,-10,-15];
% 0 is placed since the aim is to generate a step response,
% put two zeros in order to create a ramp reference

k_local = place(A,B,lambda_des)

N_gain = 1/dcgain(ss(A-B*k_local,B,C,D));

%% new system
A_new = A-B*k_local;
B_new = B*N;
C_new = C;
D_new = D;

%% new weights
Q = eye(4);
R = 1;

P = are(A_new,B_new*inv(R)*B_new',Q);

K = inv(R)*B_new'*P;

c_new = zeros(N);

for i = 1:6
    c_new(i,i) = (real(lambda(i)))/2 + 1;
end

%% global closed-loop modified system dynamics
% simulation for a node for the leader
sys_new = ss(A_new,B_new,C_new,D_new)
t = linspace(0,1000,1000);
u = zeros(1,length(t));
x0 = [1 1 0 0];
x0_transponse = lsim(sys_new,u,t,x0)


x0_bar = repmat(x0_transponse',[N 1]);
Ac = kron(eye(N),A_new) - kron(c_new*(L_cal+G_cal),B_new*K)

Bc = kron(c_new*(L_cal+G_cal),B_new*K);

Cc = eye(N*4);

Dc = zeros(N*4);

sys_global = ss(Ac,Bc,Cc,Dc);

yr = lsim(sys_global, x0_bar, t, zeros(4*N,1)); 

for i = 5:8
    figure(i)
    plot(t,yr(:,i:4:end))
    hold on
    plot(t,x0_transponse(:,i-4))
end

%% local algorithm simulation
Ts = 0.05
sys_d = c2d(sys,Ts,'zoh');
Ad = sys_d.A;
Bd = sys_d.B;

t = 0:Ts:20;

x_d = zeros(4,length(t)+1,N);
u_d = zeros(1,length(t),N);



% discrete simulation of the leader
u = zeros(1,length(t));
x0_d_transpose = lsim(sys_d,u',t,x0);
x0_d = x0_d_transpose';


for i = 1:6
    c_new(i,i) = (real(lambda(i)))/2 +1;
end

K = dlqr(Ad, Bd, Q, R);

for k = 1:length(t)
    % error calculations for command input
    epsilon = zeros(4,6);
    for i = 1:N  
        for j = 1:N
            epsilon(:,i) = epsilon(:,i) + A_cal(i,j)*(x_d(:,k,j) - x_d(:,k,i));
        end
        epsilon(:,i) =  epsilon(:,i) + G_cal(i,i)*(x0_d(:,k) - x_d(:,k,i));
    end
    return
    %dynamic update
    for i = 1:N
        u_d(1,k,i) = c(i,i)*K*epsilon(:,i);
        x_d(:,k+1,i) = Ad*x_d(:,k,i) + Bd*u_d(:,k,i);
    end
end


% plotting

for i = 9:12
    figure(i)
    plot(x_d(i-8,:,1))
    hold on
    plot(x_d(i-8,:,2))
    plot(x_d(i-8,:,3))
    plot(x_d(i-8,:,4))
    plot(x_d(i-8,:,5))
    plot(x_d(i-8,:,6))
    plot(x0_d(i-8,:),'LineWidth',1)

end

