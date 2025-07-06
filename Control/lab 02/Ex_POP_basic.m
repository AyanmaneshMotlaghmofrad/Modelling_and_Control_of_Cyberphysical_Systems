% Example of polynomial optimization in MATLAB (with SparsePOP and Sedumi)

clear all
close all

% Definition of objective function

support = [1 0 0;0 1 0;0 0 1];

coef = [-2 3 -2]';

objPoly.typeCone = 1; % always 1
objPoly.dimVar = 3;
objPoly.degree = 1;
objPoly.noTerms = 3;
objPoly.supports = support;
objPoly.coef = coef;

% Definition of constraints

% Fist constraint 19-17*x1+8*x2-14*x3+6*x1^2+3*x2^2-2*x2*x3+3*x3^2 >= 0

support1 = [0 0 0;1 0 0;0 1 0;0 0 1;2 0 0;0 2 0;0 1 1;0 0 2];

coef1 = [19 -17 8 -14 6 3 -2 3]';

ineqPolySys{1}.typeCone = 1; % 1 for inequality
ineqPolySys{1}.dimVar = 3;
ineqPolySys{1}.degree = 2;
ineqPolySys{1}.noTerms = 8;
ineqPolySys{1}.supports = support1;
ineqPolySys{1}.coef = coef1;

% Second constraint 5-x1-2*x2-x3 >= 0

support2 = [0 0 0;1 0 0;0 1 0;0 0 1];

coef2 = [5 -1 -2 -1]';

ineqPolySys{2}.typeCone = 1; % 1 for inequality
ineqPolySys{2}.dimVar = 3;
ineqPolySys{2}.degree = 1;
ineqPolySys{2}.noTerms = 4;
ineqPolySys{2}.supports = support2;
ineqPolySys{2}.coef = coef2;

% Third constraint 7-5*x2-2*x3 >= 0

support3 = [0 0 0;0 1 0;0 0 1];

coef3 = [7 -5 -2]';

ineqPolySys{3}.typeCone = -1; % -1 for equality
ineqPolySys{3}.dimVar = 3;
ineqPolySys{3}.degree = 1;
ineqPolySys{3}.noTerms = 3;
ineqPolySys{3}.supports = support3;
ineqPolySys{3}.coef = coef3;

% Bounds on variables

ubd = [2 1 1e10];  % upper bounds
lbd = [0 0 -1e10]; % lower bounds

% Set parameter for the optimization 

param.relaxOrder = 3; % Rlaxation order
param.POPsolver = 'active-set' % Refinement method

% Optimization via SparsePOP + Sedumi

[param,SDPobjValue,POP,elapsedTime,SDPsolverInfo,SDPinfo] = sparsePOP(objPoly,ineqPolySys,lbd,ubd,param)

SDPobjValue % print optimum value of the cost function (computed relaxed approximation of global minimum)
POP.objValueL % print optimum value of the cost function after refinement (computed approximation of global minimum)
POP.xVect % print optimizer (value of x = [x1 x2 x3] which provides the minimum)
POP.xVectL % print optimizer after refinement (value of x = [x1 x2 x3] which provides the minimum)