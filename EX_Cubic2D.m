% Copyright 2016, All Rights Reserved
% Code by Niall Mangan
% For Paper, "Model selection for dynamical systems via sparse regression
% and information criteria"
% by N. M. Mangan, J. N. Kutz, S. L. Brunton, and J. L. Proctor


addpath('utils')
addpath('models')
changeplot
clear all, close all, clc

filename ='Cubic2D';
plottag = 2; % for plotting measurements and cross validation plottag = 1
% for output plots only plotag =2

numvalidation = 100;
eps = 0.001; % noise

%% generate Data
n = 2;          % 2D system
A = [-.1 2; -2 -.1];
rhs = @(y)A*y.^3;   % ODE right hand side
tspan=[0:.01:25];   % time span
x0 = [2; 0];        % initial conditions
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));
[t,x]=ode45(@(t,x)rhs(x),tspan,x0,options);  % integrate

% compute Derivative 
for i=1:length(x)
    dx(i,:) = A*(x(i,:).^3)';
end
% add noise to the state variables
rng(10);
x = x + eps*randn(size(x));   

%% pool Data  (i.e., build library of nonlinear time series)
polyorder = 5;  % search space up to fifth order polynomials
usesine = 0;    % no trig functions
laurent = 0;
Theta = poolData(x,n,polyorder,usesine, laurent);
m = size(Theta,2);

Thetalib.Theta = Theta;
Thetalib.normTheta = 0;
Thetalib.dx = dx;
Thetalib.polyorder = polyorder;
Thetalib.usesine = usesine;


%% compute Sparse regression

lambdavals.numlambda = 20;
lambdavals.lambdastart = -4;
lambdavals.lambdaend = 1;
% find coefficient vectors
[Xicomb, numcoeff, lambdavec] = multiD_Xilib(Thetalib, lambdavals);

% % display coefficients
% for ii = 1:length(Xicomb)
%     Xicomb{ii}
% end

%% Simulate time series for validation

x0 = 10.^(-3 + (1+3)*rand(2,numvalidation));
% the smaller initial conditions (10^-3) are necessary to produce
% discrepencies in the data for models with lots of small coefficents 
% that are fit to noise.

for ii = 1:length(x0)
    [t,x]=ode45(@(t,x)rhs(x),tspan,x0(:,ii),options);  % integrate
    tA{ii} = t;
    xA{ii} = x + eps*randn(size(x));
end

val.x0 = x0;
val.tA = tA;
val.xA = xA;
val.options= options;

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
