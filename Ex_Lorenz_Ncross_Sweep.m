% Copyright 2016, All Rights Reserved
% Code by Niall Mangan
% For Paper, "Model selection for dynamical systems via sparse regression
% and information criteria"
% by N. M. Mangan, J. N. Kutz, S. L. Brunton, and J. L. Proctor

% This file will run a sweep over varying noise levels and save the results
% in .mat files. To plot, load each .mat file and run "AnalyzeOutput.m"
% Note: this file will take quite a while to run. 
addpath('utils')
addpath('models')
changeplot
clear all, close all, clc
savedir1=  [datestr(now, 'mmddyyyy') 'Nval_sweep/'];
mkdir(savedir1)
savetag = 1; % save everything (including plots if turned on)
plottag = 0;

%NOTE: **** the noise instance for this code is not the same as that in the
%paper, however the trends are the same. 

nvalsweep = [2 4 60 80 100];

% define system
n = 3;      % 3D system
sigma = 10;  % Lorenz's parameters (chaotic)
beta = 8/3;
rho = 28;
x0val=[-8; 8; 27];  % Initial condition

% generate training time-series without noise
tspan=[.001:.001:100];
N = length(tspan);
options = odeset('RelTol',1e-12,'AbsTol',1e-5*ones(1,n),'Events', @events_diverge);
[t,x_nonoise]=ode45(@(t,x) lorenz(t,x,sigma,beta,rho),tspan,x0val,options);

% compute Derivative
for i=1:length(x_nonoise)
    dx(i,:) = lorenz(0,x_nonoise(i,:),sigma,beta,rho);
end

rng(10); %initialize random number generator to the same place 
% measurement error value
eps=0.5; 
x = x_nonoise+eps*randn(size(x_nonoise));
  
% calculate validation time-series set without noise
tspanx=[.001:.001:1];


% Theta library params:
polyorder = 2;  % search space up to 3rd order polynomials
usesine = 0;    % no trig functions
laurent = 0;

% lambda search parameters
lambdavals.numlambda = 30;
lambdavals.lambdastart = -6;
lambdavals.lambdaend = 2;

% simulate validation time series
x0val_all = 10.^(-1 + (2+2)*rand(n,max(nvalsweep)));
for ii = 1:max(nvalsweep)
        [t2,x2]=ode45(@(t,x) lorenz(t,x,sigma,beta,rho),tspanx,x0val_all(:,ii),options);
        tA{ii} = t2;
        % add noise to validation time series
        xA{ii}= x2+ eps*randn(size(x2));
end

for i1 = 1:length(nvalsweep)
    
    nval = nvalsweep(i1); % number of validation instances
    filename = sprintf('Lorenz_numcross%d_eps%f', nval,eps);
    x0val = x0val_all(:,1:nval);
    % take subset of time series

    val.x0 = x0val;
    val.tA = tA(1:nval);
    val.options= options;
    val.xA = xA(1:nval);
    
    % build library
    Theta = poolData(x,n,polyorder,usesine, laurent);
    m = size(Theta,2);
    Thetalib.Theta = Theta;
    Thetalib.normTheta = 0;
    Thetalib.dx = dx;
    Thetalib.polyorder = polyorder;
    Thetalib.usesine = usesine;
    
    % compute Sparse regression and find coefficient vectors
    clear Xicomb numcoeff lambdavec
    [Xicomb, numcoeff, lambdavec] = multiD_Xilib(Thetalib, lambdavals);
    
    % Generate model library
    clear abserror RMSE tB xB IC AIC_rel
    for nn = 1:length(Xicomb)
        Xi = Xicomb{nn}
        clear error RMSE1 savetB savexB
        [error, RMSE1, savetB, savexB] = validateXi(Xi, Thetalib, val, plottag);
        ICtemp = ICcalculations(error', numcoeff(nn), nval);
        abserror(:,nn) = error';
        RMSE(:,nn) = RMSE1;
        tB{nn} =savetB;
        xB{nn} = savexB;
        IC(nn) = ICtemp;
    end
    % store AICc_rel in a vector for plotting purposes
    AIC_rel =cell2mat({IC.aic_c})-min(cell2mat({IC.aic_c}));
    
    if savetag>=1
        save([savedir1 filename '.mat'])
    end
end
rmpath('utils')
rmpath('models')