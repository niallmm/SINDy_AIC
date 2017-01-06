% Copyright 2016, All Rights Reserved
% Code by Niall Mangan
% For Paper, "Model selection for dynamical systems via sparse regression
% and information criteria"
% by N. M. Mangan, J. N. Kutz, S. L. Brunton, and J. L. Proctor


changeplot
clear all, close all, clc

filename ='Cubic2D';
saveplotstag = 1; % save yes =1
savetag = 1;
plottag = 2; % for plotting measurements and cross validation plottag = 1
% for output plots only plotag =2

numvalidation = 100;
eps = 0.001; % noise
%% generate Data
polyorder = 5;  % search space up to fifth order polynomials
usesine = 0;    % no trig functions
laurent = 0;

n = 2;          % 2D system
    A = [-.1 2; -2 -.1];
rhs = @(y)A*y.^3;   % ODE right hand side
tspan=[0:.01:25];   % time span
x0 = [2; 0];        % initial conditions
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));
[t,x]=ode45(@(t,x)rhs(x),tspan,x0,options);  % integrate


%% compute Derivative 
eps = 0.001;      % noise strength
for i=1:length(x)
    dx(i,:) = A*(x(i,:).^3)';
end
dx = dx + eps*randn(size(dx));   % add noise

%% pool Data  (i.e., build library of nonlinear time series)
Theta = poolData(x,n,polyorder,usesine, laurent);
m = size(Theta,2);

%% compute Sparse regression

x0 = 2*rand(2,numvalidation);
for ii = 1:length(x0)
    [t,x]=ode45(@(t,x)rhs(x),tspan,x0(:,ii),options);  % integrate
    xA{ii} = x;
    xA{ii} = x + eps*randn(size(x));
%     for i=1:length(x)
%     dxt(:,i) = A*(x(i,:).^3)';
%     end
%     dxA{ii}= dxt' + eps*randn(size(dxt'));   % add noise
end
%%
Thetalib.Theta = Theta;
Thetalib.normTheta = 0;
Thetalib.dx = dx;
Thetalib.polyorder = polyorder;
Thetalib.usesine = usesine;

crossval.x0 = x0;
crossval.tvec = t;
crossval.xA = xA;
crossval.options= options;

lassovals.numlasso = 20;
lassovals.lassostart = -4;
lassovals.lassoend = 1;

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

