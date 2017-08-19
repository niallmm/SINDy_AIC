% Copyright 2016, All Rights Reserved
% Code by Niall Mangan
% For Paper, "Model selection for dynamical systems via sparse regression
% and information criteria"
% by N. M. Mangan, J. N. Kutz, S. L. Brunton, and J. L. Proctor

addpath('utils')
addpath('models')
changeplot
clear all, close all, clc



filename = 'Lorenz';
plottag = 2; % for plotting measurements and cross validation plottag = 1
% for output plots only plotag =2

numvalidation = 100; % number of crossvalidation experiments
eps = 0.001; % noise
%% generate Data
n = 3;      % 3D system
sigma = 10;  % Lorenz's parameters (chaotic)
beta = 8/3;
rho = 28;
x0=[-8; 8; 27];  % Initial condition

% Integrate
tspan=[.001:.001:100];
N = length(tspan);
options = odeset('RelTol',1e-12,'AbsTol',1e-5*ones(1,n),'Events', @events_diverge);
[t,x]=ode45(@(t,x) lorenz(t,x,sigma,beta,rho),tspan,x0,options);

% compute Derivative

for i=1:length(x)
    dx(i,:) = lorenz(0,x(i,:),sigma,beta,rho);
end
% add noise to the state variables
rng(10);
%NOTE: **** This is not the same noise instance as the simulation in the paper.
% changing the noise initialization to other values sometimes
%causes the AIC ranking to fail to score the correct model with the lowest
%score (ie rng(8), rng(9)). Sometimes there are many other models that recieve supported
% scores (as is true for rng(6) intialization) We believe this is because the
% validation time series extend beyond the Lyapunov time of the system (t = 5>1).

x = x+eps*randn(size(x));

%% pool Data  (i.e., build library of nonlinear time series)
polyorder = 2;  % search space up to 3rd order polynomials
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
lambdavals.numlambda = 30;
lambdavals.lambdastart = -6;
lambdavals.lambdaend = 2;

% find coefficient vectors
[Xicomb, numcoeff, lambdavec] = multiD_Xilib(Thetalib, lambdavals);

% % display coefficients
% for ii = 1:length(Xicomb)
%     Xicomb{ii}
% end

%% calculate validation data for new intial conditions.

tspanx=[.001:.001:5];
x0 = 10.^(-1 + (3+1)*rand(n,numvalidation));
for ii = 1:numvalidation
    [t2,x2]=ode45(@(t,x) lorenz(t,x,sigma,beta,rho),tspanx,x0(:,ii),options);
    tA{ii} = t2;
    xA{ii}= x2 + eps*randn(size(x2));
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