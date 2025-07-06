clc
clear variables  
close all
format compact
%%
% load exp_data.mat
% 
% N = length(u)
% size_theta = 5;
% n = 2
%%
s = tf('s');
z = tf('z',-1);

Gp = 100/(s^2 + 1.2*s + 1)

Gp_d = c2d(Gp,1,'zoh')

[num_d, den_d] = tfdata(Gp_d,'v')
theta_real = [den_d num_d]'
%% simulation
H = 50;% in the noise free case 7 samplesa are enough

u = rand([H 1]);

w = lsim(Gp_d,u);

eta = 5*randn([H 1]); % gaussian noise assumption

y_tilde = w + eta;

N = H;
size_theta = 5;
%%
for c = 1:N-2

    % equality support matrix
    %initialization
    supp_eq = zeros(9, size_theta+ N);

    %theta rows
    supp_eq(2:6, 1:size_theta) = eye(size_theta); 
    supp_eq(8:end, 1:2) = eye(2);
    
    % eta rows
    supp_eq(7:9, size_theta+c:size_theta+c+2) = flip(eye(3));

% % xi rows
% supp_eq(10:12, size_theta(2)+c+N:size_theta(2)+c+2+N) = flip(eye(3));
% 
    %equality coefficient vector
    coeffs_eq = [y_tilde(c+2), y_tilde(c+1), y_tilde(c),...
        -u(c+2), -u(c+1), -u(c), -1, -1, -1]';
    
    % creating the structures
    ineqPolySys{c}.noTerms = 9;
    ineqPolySys{c}.degree = 2;
    ineqPolySys{c}.dimVar = size_theta + N;
    ineqPolySys{c}.typeCone = -1; % equality
    ineqPolySys{c}.supports = supp_eq;
    ineqPolySys{c}.coef = coeffs_eq;
end


%% setting parameters
delta_eta = 5;
% lower bound
lbd = [-1e10*ones(size_theta,1); -delta_eta*ones(N,1)];
% upper bound
ubd = -lbd;

param.relaxOrder = 1;
param.POPsolver = 'active-set';

%% upper bounds

for i = 1:size_theta

    support = zeros(1,size_theta + N);
    support(i) = 1;

    objPoly.typeCone = 1; % always 1
    objPoly.dimVar = size_theta + N;
    objPoly.degree = 1;
    objPoly.noTerms = 1;
    objPoly.supports = support;
    objPoly.coef = 1;


    [param,SDPobjValue,POP,elapsedTime,SDPsolverInfo,SDPinfo] = sparsePOP(objPoly,ineqPolySys,lbd,ubd,param)
    
    PUI_lbd(i) = POP.xVectL(i);
end

for i = 1:size_theta

    support = zeros(1,size_theta + N);
    support(i) = 1;

    objPoly.typeCone = 1; % always 1
    objPoly.dimVar = size_theta + N;
    objPoly.degree = 1;
    objPoly.noTerms = 1;
    objPoly.supports = support;
    objPoly.coef = -1;


    [param,SDPobjValue,POP,elapsedTime,SDPsolverInfo,SDPinfo] = sparsePOP(objPoly,ineqPolySys,lbd,ubd,param)
    
    PUI_ubd(i) = POP.xVectL(i);
end

%% central model
[[1 PUI_lbd]', theta_real, [1 PUI_ubd]']

theta = mean([PUI_lbd',PUI_ubd'],2)

num = theta(3:end)';
den = [1 theta(1:2)']

G_d = tf(num,den,-1)
dcgain(G_d)
%% forcing dcgain 118


% for c = N-1

    supp_eq = zeros(6,size_theta + N);
    supp_eq(2:end, 1:size_theta) = eye(5);

    coef_eq = [-118 -118 -118 1 1 1]'; 
    ineqPolySys{c+1}.noTerms = 6;
    ineqPolySys{c+1}.degree = 1;
    ineqPolySys{c+1}.dimVar = size_theta + N;
    ineqPolySys{c+1}.typeCone = -1; % equality
    ineqPolySys{c+1}.supports = supp_eq;
    ineqPolySys{c+1}.coef = coef_eq;

%%
for i = 1:size_theta

    support = zeros(1,size_theta + N);
    support(i) = 1;

    objPoly.typeCone = 1; % always 1
    objPoly.dimVar = size_theta + N;
    objPoly.degree = 1;
    objPoly.noTerms = 1;
    objPoly.supports = support;
    objPoly.coef = 1;


    [param,SDPobjValue,POP,elapsedTime,SDPsolverInfo,SDPinfo] = sparsePOP(objPoly,ineqPolySys,lbd,ubd,param)
    
    PUI_lbd_dcgain(i) = POP.xVectL(i);
end

for i = 1:size_theta

    support = zeros(1,size_theta + N);
    support(i) = 1;

    objPoly.typeCone = 1; % always 1
    objPoly.dimVar = size_theta + N;
    objPoly.degree = 1;
    objPoly.noTerms = 1;
    objPoly.supports = support;
    objPoly.coef = -1;


    [param,SDPobjValue,POP,elapsedTime,SDPsolverInfo,SDPinfo] = sparsePOP(objPoly,ineqPolySys,lbd,ubd,param)
    
    PUI_ubd_dcgain(i) = POP.xVectL(i);
end

%% 
[PUI_lbd_dcgain' PUI_ubd_dcgain']

theta = mean([PUI_lbd_dcgain', PUI_ubd_dcgain'],2)

num = theta(3:end)';
den = [1 theta(1:2)']

G_d = tf(num,den,-1)
dcgain(G_d)