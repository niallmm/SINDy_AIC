% Copyright 2016, All Rights Reserved
% Code by Niall Mangan
% For Paper, "Model selection for dynamical systems via sparse regression
% and information criteria"
% by N. M. Mangan, J. N. Kutz, S. L. Brunton, and J. L. Proctor

% This code:
% 1) runs a model to produce data for a simple 1-state variable system
% 2) creates a data library (for SINDy)and screening matrix (for
% combinatorial least squares regression)
% 3) performs combinatorial L2 regression and alternating L2 and
% thresholding
% 4) calls functions to validate and calculate IC for all produced models


addpath('utils')
addpath('models')
changeplot
clear all, close all, clc

filename = '1Dpolynomial';

plottag = 2; % for plotting measurements and validation plottag = 1
% for output plots only plotag =2
% generate Data
polyorder = 5;  % search space up to sixth order polynomials
usesine = 0;    % no trig functions
n = 1;          % 1D system

eps = 0.001; % low noise

% number of crossvalidations
numvalidation = 100;
tspan_t = 0:0.25:5;

% model
Xitrue = [0 1 0 -0.2 -0.1 0]'; % define coefficients up to 5th order
rhs = @(y)(Xitrue(1) + Xitrue(2)*y +Xitrue(3)*y^2+ Xitrue(4)*y.^3 ...
    +Xitrue(5)*y.^4 + Xitrue(6)*y.^5);   % ODE right hand side

% data aquisition.
tspan=0:0.25:5;   % time span
x0 = [0.1, 0.5, 5];   % initial conditions
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));
t = [];
x = [];
for nn = 1:length(x0)
    [t1,x1]=ode45(@(t,x)rhs(x),tspan,x0(nn),options);  % integrate
    t = [t; t1];
    x = [x; x1];
end

% compute Derivative
for i=1:length(x)
    dx(i,:) = rhs(x(i,:));
end

rng(10);
% add noise
noisevec = eps*randn(size(x));
x = x + noisevec;

% simulate time series for validation
x0test = logspace(-3,0.5, numvalidation);
for nn = 1:numvalidation
    [tA, xA]=ode45(@(t,x)rhs(x),tspan_t,x0test(:,nn),options);  % integrate
    savetA{nn} = tA;
    savexA{nn} = xA + eps*randn(size(xA));
end
% save into a structure to use in functions
val.x0 = x0test;
val.tA = savetA;
val.xA = savexA;
val.options= options;

%% % pool Data  (i.e., build library of nonlinear time series)
Theta = poolData(x,n,polyorder,usesine,0);
m = size(Theta,2);

Thetalib.Theta = Theta;
Thetalib.normTheta = 0;
Thetalib.dx = dx;
Thetalib.polyorder = polyorder;
Thetalib.usesine = usesine;

%% screening permutation matrix
permM = []; numterms=[];
for nn = 1:m
    permM = [permM; permpos(nn,m)];
    numterms = [numterms; nn*ones(nchoosek(m,nn),1)];
end
[tf index] = ismember(permM, [0 1 0 1 1 0], 'rows');
correct = find(index);
[numperms, junk] = size(permM);

% do least squares on all permuations 
 for ll = 1:numperms
    Thetascreen = Theta*diag(permM(ll,:));
     Xiperm(:,ll) = Thetascreen\dx;
 end
% least square permuted validation
for mm = 1:numperms
[abserror, RMSE, savetB, savexB] = validateXi(Xiperm(:,mm), Thetalib, val, plottag);
LSPermute_abserror{mm} = abserror;
IC = ICcalculations(abserror', numterms(mm), numvalidation);
LSPermuteaic_c(mm) = IC.aic_c;
end
LSPermuteaic_c_rel = LSPermuteaic_c-min(LSPermuteaic_c);
%% plot permuations
if plottag>1
    fulAIC_rel = figure(3);
    plot(numterms, LSPermuteaic_c_rel, 'ok')
    hold on
    zoomAIC_rel = figure(4);
    plot(numterms, LSPermuteaic_c_rel, 'ok')
    hold on
end

%% compute Sparse regression

lambdavals.numlambda = 100;
lambdavals.lambdastart = -6;
lambdavals.lambdaend = 1;
% find coefficient vectors
[Xicomb, numcoeff, lambdavec] = multiD_Xilib(Thetalib, lambdavals);

%% Perform validation and calculate aic for each model
clear abserror RMSE tB xB IC
for nn = 1:length(Xicomb)
    Xi = Xicomb{nn};
    clear error RMSE1 savetB savexB
    [error, RMSE1, savetB, savexB] = validateXi(Xi, Thetalib, val, plottag);
    ICtemp = ICcalculations(error', numcoeff(nn), numvalidation);
    abserror(:,nn) = error';
    RMSE(:,nn) = RMSE1;
    tB{nn} =savetB;
    xB{nn} = savexB;
    IC(nn) = ICtemp;
end


AIC_rel =cell2mat({IC.aic_c})-min(cell2mat({IC.aic_c}));
% plot numterms vs AIC plot
AnalyzeOutput

rmpath('utils')
rmpath('models')
