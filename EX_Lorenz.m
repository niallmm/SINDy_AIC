% Copyright 2016, All Rights Reserved
% Code by Niall Mangan
% For Paper, "Model selection for dynamical systems via sparse regression
% and information criteria"
% by N. M. Mangan, J. N. Kutz, S. L. Brunton, and J. L. Proctor


changeplot
clear all, close all, clc

filename = 'Lorenz';
saveplotstag = 1; % save yes =1
savetag = 1;
plottag = 2; % for plotting measurements and cross validation plottag = 1
% for output plots only plotag =2

ncross = 100; % number of crossvalidation experiments
eps = 0.001; % noise
%% generate Data

sigma = 10;  % Lorenz's parameters (chaotic)
beta = 8/3;
rho = 28;

n = 3;

x0=[-8; 8; 27];  % Initial condition

% Integrate
tspan=[.001:.001:100];
N = length(tspan);
options = odeset('RelTol',1e-12,'AbsTol',1e-5*ones(1,n));
[t,x]=ode45(@(t,x) lorenz(t,x,sigma,beta,rho),tspan,x0,options);

%% compute Derivative

for i=1:length(x)
    dx(i,:) = lorenz(0,x(i,:),sigma,beta,rho);
end
% dx = dx + eps*randn(size(dx));
x = x+eps*randn(size(x));

%% pool Data  (i.e., build library of nonlinear time series)
polyorder = 2;  % search space up to 3rd order polynomials
usesine = 0;    % no trig functions
laurent = 0;
Theta = poolData(x,n,polyorder,usesine, laurent);
m = size(Theta,2);
%%
% calculate cros validation data for new intial conditions.

tspanx=[.001:.001:5];
x0cross = 20*randn(n,ncross);
for ii = 1:ncross
    [t2,x2]=ode45(@(t,x) lorenz(t,x,sigma,beta,rho),tspanx,x0cross(:,ii),options);
    xA{ii}= x2 + eps*rand(size(x2));
end
    
%% compute Sparse regression
Thetalib.Theta = Theta;
Thetalib.normTheta = 0;
Thetalib.dx = dx;
Thetalib.polyorder = polyorder;
Thetalib.usesine = usesine;

crossval.x0 = x0cross;
crossval.tvec = tspanx;
crossval.xA = xA;
crossval.options= options;

lassovals.numlasso = 30;
lassovals.lassostart = -6;
lassovals.lassoend = 2;

Lasso_CV =multiD_pareto_crossn(Thetalib,crossval,lassovals, plottag);
                     
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


