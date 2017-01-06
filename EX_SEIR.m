% Copyright 2016, All Rights Reserved
% Code by Niall Mangan
% For Paper, "Model selection for dynamical systems via sparse regression
% and information criteria"
% by N. M. Mangan, J. N. Kutz, S. L. Brunton, and J. L. Proctor


changeplot
clear all, close all, clc

filename = 'SEIR';
saveplotstag = 1; % save yes =1
savetag = 1;
plottag = 2; % for plotting measurements and cross validation plottag = 1
% for output plots only plotag =2



%% 
n = 3; % number of parameters
% generate Data
N  = 250; % number of time steps
eps=  1; % noise

% all others are zero
% Transfer Parameters
B_SE = 0.3;
B_EI = 0.4;
B_IR = 0.04;


Ntot = 1e4; % total population

% Initial Conditions
S(1) = 0.99*Ntot; % number of suceptibles in population
E(1) = 0.01*Ntot;
I(1) = 0;


% disease tranfer model
for ii =2:N
    S(ii) = S(ii-1) - B_SE*S(ii-1)*I(ii-1)/Ntot;
    E(ii) = E(ii-1) + B_SE*S(ii-1)*I(ii-1)/Ntot - B_EI*E(ii-1);
    I(ii) = I(ii-1) + B_EI*E(ii-1) - B_IR*I(ii-1);
end

% create x and dx matrices with all variables:
x = [S(1:end-1)' E(1:end-1)' I(1:end-1)'];% R(1:end-1)'];
% x = I(1:N-1);
dx = [S(2:end)' E(2:end)' I(2:end)'];% R(2:end)'];

x = x+eps*randn(size(x));
figure(6)
plot(x/Ntot, 'o')
xlabel('time step')
ylabel(['% of population size: ' num2str(Ntot)])
legend('S', 'E', 'I', 'R')
legend('boxoff')

%% pool Data
polyorder = 2;  % search space up to 3rd order polynomials
usesine = 0;    % no trig functions
laurent = 0;
Theta = poolData(x,n,polyorder,usesine, laurent);
m = size(Theta,2);

% calculate cros validation data for new intial conditions.
ncross = 100; % number of crossvalidation experiments

x0cross = 0.5*rand(n,ncross)*Ntot;
for jj = 1:ncross
    % initialize
    S = zeros(N,1); % susceptibles
    E = zeros(N,1); % Latent/exposed
    I = zeros(N,1); % infected
  %  R = zeros(N,1); % recovered
    % Initial Conditions
    S(1) = x0cross(1, jj); % number of suceptibles in population
    E(1) = x0cross(2, jj);
    I(1) = x0cross(3, jj);
 %   R(1) = x0(4, jj)*Ntot;
    
    for ii =2:N
        S(ii) = S(ii-1) - B_SE*S(ii-1)*I(ii-1)/Ntot;
        E(ii) = E(ii-1) + B_SE*S(ii-1)*I(ii-1)/Ntot - B_EI*E(ii-1);
        I(ii) = I(ii-1) + B_EI*E(ii-1) - B_IR*I(ii-1);
   %     R(ii) = R(ii-1) + B_IR*I(ii-1);
    end
    % create x and dx matrices with all variables:
    x2= [S(1:end-1) E(1:end-1) I(1:end-1)];
    xA{jj} = x2 +eps*rand(size(x2));
    dxA{jj} = [S(2:end) E(2:end) I(2:end)]';
    
end



%%
% % normalize the columns of Theta
for k=1:m
    normTheta(k) = norm(Theta(:,k));
    Theta(:,k) = Theta(:,k)/normTheta(k);
end 
%%
% t = N; % for discrete time, feed in the number of steps.
Thetalib.Theta = Theta;
Thetalib.normTheta = normTheta;
Thetalib.dx = dx;
Thetalib.polyorder = polyorder;
Thetalib.usesine = usesine;

crossval.x0 = x0cross;
crossval.tvec = N;
crossval.xA = xA;
crossval.options= 0;

lassovals.numlasso = 20;
lassovals.lassostart = -4;
lassovals.lassoend = 4;

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

