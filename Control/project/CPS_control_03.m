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
lambda_obs = [-15 -20];
L = place(A_unstable',C_unstable',lambda_obs)';

lambda_ctrl = [0 -15]; 

s= tf('s');
K_loc = place(A_unstable,B_unstable,lambda_ctrl);
N_loc = 1/dcgain(s*tf(ss(A_unstable-B_unstable*K_loc,B_unstable,C_unstable,D_unstable)));

A = [A_unstable -B_unstable*K_loc;L*C_unstable A_unstable-B_unstable*K_loc-L*C_unstable];
B = [B_unstable;B_unstable]*N_loc;
C = [C_unstable zeros(1,n_sys)];
D = 0;

sys_node = ss(A,B,C,D);
t = linspace(0,100,5000);
u = zeros(1,length(t));
Ts = t(2) - t(1);
u(1) = 1/Ts;
x0 = [0 0 0 0]';



[y0_transpose,~,x0_transpose] = lsim(sys_node,u,t,x0);
plot(t,y0_transpose)

x0_bar = repmat(x0_transpose',[N 1]);

%% First Network topology 
% Chain topology
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

for i = 1:2
    if(i == 1 || i ==3)
         figure(i)
         plot(t,xr(:,i:4:end)*708.2700)
         legend
    else
        figure(i)
         plot(t,xr(:,i:4:end))
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
         legend
    else
        figure(i)
         plot(t,xr(:,i-2:4:end))
         legend
    end
end

%% Third Network topology 
% leader star topology
n = 4;
G_cal = diag([1 1 1 1 1 1]);
A_cal = zeros(6);

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



for i = 5:6
    if(i-4 == 1 || i-4 ==3)
         figure(i)
         plot(t,xr(:,i-4:4:end)*708.2700)
         legend
    else
        figure(i)
         plot(t,xr(:,i-4:4:end))
         legend
    end
end

%% Forth Network topology 
% fully connected-nodes topology
n = 4;
G_cal = diag([1 0 0 0 0 0]);
A_cal = ones(6);

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

n = 4;
G_cal = diag([1 0 0 0 0 0]);
A_cal = ones(6);

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
         legend
    else
        figure(i)
         plot(t,xr(:,i-8:4:end))
         legend
    end
end

%% sixth Network topology 
% unidirected chain topology
n = 4;
A_cal = [0 1 0 0 0 0;
         1 0 1 0 0 0;
         0 1 0 1 0 0;
         0 0 1 0 1 0;
         0 0 0 1 0 1;
         0 0 0 0 1 0];
G_cal = diag([1 0 0 0 0 0]);
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

% unstable system
A_unstable = [0       1
     880.87  0];
B_unstable = [0 -9.9453]';
C_unstable = [708.27 0];
D_unstable = 0;

sys_node = ss(A_unstable,B_unstable,C_unstable,D_unstable);

rank(obsv(A_unstable,C_unstable));
rank(ctrb(A_unstable,B_unstable));

% local observer and controller
lambda_obs = [-15 -20];
L = place(A_unstable',C_unstable',lambda_obs)';

lambda_ctrl = [0 -15]; 

s= tf('s');
K_loc = place(A_unstable,B_unstable,lambda_ctrl);
N_loc = 1/dcgain(s*tf(ss(A_unstable-B_unstable*K_loc,B_unstable,C_unstable,D_unstable)));

% leader must have the observer separately
A0 = [A_unstable -B_unstable*K_loc;L*C_unstable A_unstable-B_unstable*K_loc-L*C_unstable];
B0 = [B_unstable;B_unstable]*N_loc;
C0 = [C_unstable zeros(1,n_sys)];
D0 = 0;

sys0 = ss(A0,B0,C0,D0);

t = linspace(0,10,5000);
u = zeros(1,length(t));
Ts = t(2) - t(1);
u(1) = 1/Ts;
x0 = [0 0 0 0]';



[y0_transpose,~,x0_transpose] = lsim(sys0,u,t,x0);
plot(t,y0_transpose)

x0_bar = repmat(x0_transpose(:,1:2)',[N 1]);

%% First Network topology global observer
% Chain topology
close all
n = 2;

G_cal = diag([10 0 0 0 0 0]);
A_cal = [0 0 0 0 0 0;
         1 0 0 0 0 0;
         0 1 0 0 0 0;
         0 0 1 0 0 0;
         0 0 0 1 0 0;
         0 0 0 0 1 0]*10;
% later use different A_cal
% for the controller and the observer

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
    c_c(i,i) = (real(lambda(i)))/2+1;
end

% global observer
Q_o = eye(n); % two states
R_o = 1;

P_o = are(A_unstable',C_unstable'*inv(R_o)*C_unstable,Q_o);
F = P_o*C_unstable'*inv(R_o);

c_o = zeros(N); %coupling gain for the controller

for i = 1:6
    c_o(i,i) = (real(lambda(i)))/2 + 1;
end


% global matrices
Ao = kron(eye(N),A_unstable)-kron(c_o*(L_cal+G_cal),F*C_unstable);
Ac = [kron(eye(N),A_unstable), -kron(c_c*(L_cal+G_cal),B_unstable*K);
      kron(c_o*(L_cal+G_cal),F*C_unstable), Ao-kron(c_c*(L_cal+G_cal),B_unstable*K)];
Bc = [kron(c_c*(L_cal +G_cal),B_unstable*K);kron(c_c*(L_cal +G_cal),B_unstable*K)];
Cc = eye(N*2*n);
Dc = zeros(N*2*n,N*n);

sys_global = ss(Ac,Bc,Cc,Dc);
[yr,~,xr]= lsim(sys_global, x0_bar, t, zeros(4*N,1)); 

for i = 1:2
    if(i == 1)
         figure(i)
         plot(t,xr(:,i:2:end)*708.2700)
         legend
    else
        figure(i)
         plot(t,xr(:,i:2:end))
         legend
    end
end


%% local algorithm simulation
% tree topology
G_cal = diag([100 0 0 0 0 0]);
A_cal = [0 0 0 0 0 0;
         1 0 0 0 0 0;
         0 1 0 0 0 0;
         0 0 1 0 0 0;
         0 0 0 1 0 0;
         0 0 0 0 1 0]*100;
% later use different A_cal
% for the controller and the observer

D_cal= diag(sum(A_cal,2));
L_cal = D_cal-A_cal;
lambda = eig(L_cal+G_cal);
Ts = 0.01;
sys_d = c2d(sys_node,Ts,'zoh');
Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;
Dd = sys_d.D;


% discrete simulation of the leader

t = 0:Ts:1;
u = zeros(1,length(t));
Ts = t(2) - t(1);
u(1) = 1/Ts;
x0 = [0 0 0 0]';
[y,~,x0_d_transpose] = lsim(sys0,u,t,x0);
x0_d = x0_d_transpose';


% initialization the states and the input signal
x_d = zeros(4,length(t)+1,N); %[x;x_hat]
u_d = zeros(1,length(t),N);


% Distributed control protocol
Q_c = eye(n);
R_c = 1;

K = dlqr(Ad, Bd, Q_c, R_c);

c_c = zeros(N); %coupling gain for the controller

for i = 1:6
    c_c(i,i) = (real(lambda(i)))/2+1;
end

% global observer
Q_o = eye(n); % two states
R_o = 1;

F = dlqr(Ad', Cd', Q_o, R_o)';

c_o = zeros(N); %coupling gain for the controller

for i = 1:6
    c_o(i,i) = (real(lambda(i)))/2 + 1;
end


for k = 1:length(t)
    % error calculations for command input
    
    epsilon_hat = zeros(2,6);
    epsilon = zeros(1,6);

    % 2 stands for the number of states of  the observer and 6 for the number of nodes
    
    %calculating epsilon_hat for the plant
    for i = 1:N  
        for j = 1:N
            epsilon_hat(:,i) = epsilon_hat(:,i) + A_cal(i,j)*(x_d(3:4,k,j) - x_d(3:4,k,i));
        end
        epsilon_hat(:,i) =  epsilon_hat(:,i) + G_cal(i,i)*(x0_d(3:4,k) - x_d(3:4,k,i));
    end
    
    % calculating epsilon for the observer ( action on y_tilde)
    for i = 1:N  
        for j = 1:N
            epsilon(:,i) = epsilon(:,i) + A_cal(i,j)*(Cd*x_d(1:2,k,j) - Cd*x_d(1:2,k,i));
        end
        epsilon(:,i) =  epsilon(:,i) + G_cal(i,i)*(Cd*x0_d(1:2,k) - Cd*x_d(1:2,k,i));
    end

    %dynamic update
    for i = 1:N
        u_d(1,k,i) = c_c(i,i)*K*epsilon_hat(:,i);
        
        % actual state equation
        x_d(1:2,k+1,i) = Ad*x_d(1:2,k,i) + Bd*u_d(1,k,i);        
        x_d(3:4,k+1,i) = Ad*x_d(3:4,k,i) + Bd*u_d(1,k,i) - c_o(i,i)*F*epsilon(:,i);
    
    end    
end



% plotting

for i = 1:4 % each state
    figure(i)
    hold on
    for j = 1:6 % for all the nodes
        if(i == 1 || i == 3)
             plot(x_d(i,:,j)*708.2700)
        else
            plot(x_d(i,:,j))
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












