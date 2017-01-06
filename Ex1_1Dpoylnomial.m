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
% 4) calls functions to cross-validate and calculate IC for all produced models

changeplot
clear all, close all, clc

filename = '1Dpolynomial';

saveplotstag = 1; % save yes =1
savetag = 1;
plottag = 2; % for plotting measurements and cross validation plottag = 1
% for output plots only plotag =2
% generate Data
polyorder = 5;  % search space up to sixth order polynomials
usesine = 0;    % no trig functions
n = 1;          % 1D system

eps = 0.001; % low noise

% number of crossvalidations
numvalidation = 100;
tspan_t = 0:0.1:5;

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
for kk = 1:length(x0)
    [t1,x1]=ode45(@(t,x)rhs(x),tspan,x0(kk),options);  % integrate
    t = [t; t1];
    x = [x; x1];
end

% compute Derivative
for i=1:length(x)
    dx(i,:) = rhs(x(i,:));
end
% add noise
noisevec = eps*randn(size(x));
x = x + noisevec;

if plottag == 1

end
% simulate time series for cross validation
x0test = logspace(-2,2, numvalidation);
for kk = 1:numvalidation
    [tA, xA]=ode45(@(t,x)rhs(x),tspan_t,x0test(:,kk),options);  % integrate
    savetA{kk} = tA;
    savexA{kk} = xA + eps*randn(size(xA));
end
 %% % pool Data  (i.e., build library of nonlinear time series)
Theta = poolData(x,n,polyorder,usesine,0);
m = size(Theta,2);
%% screening permutation matrix
permM = []; numterms=[];
for kk = 1:m
    permM = [permM; permpos(kk,m)];
    numterms = [numterms; kk*ones(nchoosek(m,kk),1)];
end
[tf index] = ismember(permM, [0 1 0 1 1 0], 'rows');
correct = find(index);
[numperms, junk] = size(permM);

% do least squares on all permuations 
 for ll = 1:numperms
    Thetascreen = Theta*diag(permM(ll,:));
     Xiperm(:,ll) = Thetascreen\dx;
 end
% least square permuted Cross-validation
for mm = 1:numperms
[LSPermute_abserror{mm}, RMSE, savetB, savexB] = ...
crossvalidateXi(Xiperm(:,mm), rhs, tspan_t, x0test, polyorder, usesine, ...
options, plottag, savetA, savexA);
IC = ICcalculations(LSPermute_abserror{mm}, numterms(mm), numvalidation);
LSPermuteaic_c(mm) = IC.aic_c;
end
LSPermuteaic_c_rel = LSPermuteaic_c-min(LSPermuteaic_c);
%%
fulAIC_rel = figure(3);
plot(numterms, LSPermuteaic_c_rel, 'ok')
hold on
zoomAIC_rel = figure(4);
plot(numterms, LSPermuteaic_c_rel, 'ok')
hold on

%% compute Sparse regression
Thetalib.Theta = Theta;
Thetalib.normTheta = 0;
Thetalib.dx = dx;
Thetalib.polyorder = polyorder;
Thetalib.usesine = usesine;

crossval.x0 = x0test;
crossval.tvec = tspan_t;
crossval.xA = savexA;
crossval.options= options;

lassovals.numlasso = 100;
lassovals.lassostart = -6;
lassovals.lassoend = 1;
[Lasso_CV] = multiD_pareto_crossn(Thetalib,crossval,lassovals, plottag);

savetB=Lasso_CV.tB;
savexB =Lasso_CV.xB;
abserror=Lasso_CV.abserror;
RMSE = Lasso_CV.RSME;
IC=Lasso_CV.IC;
numcoeff = Lasso_CV.numcoeff;
Xicomb = Lasso_CV.Xi;
lambdavec = Lasso_CV.lambda;

AIC_rel =cell2mat({IC.aic})-min(cell2mat({IC.aic}));
AnalyzeOutput
