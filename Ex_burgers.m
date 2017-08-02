% Copyright 2016, All Rights Reserved
% Code by J.N. Kutz and  N. M. Mangan
% For Paper, "Model selection for dynamical systems via sparse regression
% and information criteria"
% by N. M. Mangan, J. N. Kutz, S. L. Brunton, and J. L. Proctor

% This code runs simulations of Burgers' equations and the PDEs found using
% PDE-find on time series data from Burgers' equation for 100 initial 
% conditions. Then it calculates the relative AICc for each model.
% The PDEs definition functions are stored in the models/Burgers_Proc folder. 
% The PDEs are solved in fourier space.
addpath('utils')
addpath('models/Burgers_Proc')
changeplot
clear all, close all, clc

% Burgers' equation
%  u_t + u*u_x - eps*u_xx =0

plottag = 0; % turn plotting on and off. 
% varying parameters for initial conditions.
amp = linspace(0.2, 2, 10); % vary the amplitude of the inital pulse
mid = linspace(-3, 3, 10); % vary the location of the initial pulse

% setup
dt=0.1;
t=0:dt:10; eps=0.1;
L=16; n=256;
x2=linspace(-L/2,L/2,n+1); x=x2(1:n); k=(2*pi/L)*[0:n/2-1 -n/2:-1].';
k2=fftshift(k);
for kk = 1:10
    for ii = 1:10
        % initialize a pulse
        u = amp(kk)*exp(-(x+2).^2).'; 
        ut=fft(u); 
        
        % solve in fft space
        [t,utsol]=ode45('burgers_rhs',t,ut,[],k,eps);
        
        % transform solution back for each time step
        for jj=1:length(t)
            usol(:,jj)=ifft( utsol(jj,1:n).' );
        end
        
        % plot the solution in waterfall plot
        if plottag ==1
        figure(1),
        waterfall(x,t,real(usol.')); colormap([0 0 0]);
        view(15,35), set(gca,'Fontsize',[12])
        set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2],'Ylim',[0 10],'Ytick',[0 10],'Zlim',[0 1],'Ztick',[0 0.5 1])
        end
        
        % add noise (same as the noise size in the original PDE-find paper)
        % it is based on the magnitude of the data
        noisesize =  0.001814019805;
        usol = usol+noisesize*randn(size(usol));
        
        %% OTHER MODELS
        % solve all 7 models in the same way.
        [t,utsol1]=ode45('burgers_rhs1',t,ut,[],k,eps);
        [t,utsol2]=ode45('burgers_rhs2',t,ut,[],k,eps);
        [t,utsol3]=ode45('burgers_rhs3',t,ut,[],k,eps);
        [t,utsol4]=ode45('burgers_rhs4',t,ut,[],k,eps);
        [t,utsol5]= ode45('burgers_rhs5',t,ut,[],k,eps);
        [t,utsol6]= ode45('burgers_rhs6',t,ut,[],k,eps);
        [t,utsol7]= ode45('burgers_rhs7',t,ut,[],k,eps);
        % transform all models back to time-space
        for jj=1:length(t)
            usol1(:,jj)=ifft( utsol1(jj,1:n).' );
            usol2(:,jj)=ifft( utsol2(jj,1:n).' );
            usol3(:,jj)=ifft( utsol3(jj,1:n).' );
            usol4(:,jj)=ifft( utsol4(jj,1:n).' );
            usol5(:,jj)=ifft( utsol5(jj,1:n).' );
            usol6(:,jj)=ifft( utsol6(jj,1:n).' );
            usol7(:,jj)=ifft( utsol7(jj,1:n).' );
        end
        
        % calculate the error between ground truth (burgers solution) and
        % each model
        abserror1m(ii,kk) = sum(sum(abs(usol-usol1)))/numel(usol);
        abserror2m(ii,kk) = sum(sum(abs(usol-usol2)))/numel(usol);
        abserror3m(ii,kk) = sum(sum(abs(usol-usol3)))/numel(usol);
        abserror4m(ii,kk) = sum(sum(abs(usol-usol4)))/numel(usol);
        abserror5m(ii,kk) = sum(sum(abs(usol-usol5)))/numel(usol);
        abserror6m(ii,kk) = sum(sum(abs(usol-usol6)))/numel(usol);
        abserror7m(ii,kk) = sum(sum(abs(usol-usol7)))/numel(usol);
    end
end
%% Calculate AICc for each discovered model

abserror1 = reshape(abserror1m, [],1);
IC1 = ICcalculations(abserror1, 15, 100);
abserror2 = reshape(abserror2m, [],1);
IC2 = ICcalculations(abserror2, 14, 100);
abserror3 = reshape(abserror3m, [],1);
IC3 = ICcalculations(abserror3, 13, 100);
abserror4 = reshape(abserror4m, [],1);
IC4 = ICcalculations(abserror4, 9, 100);
abserror5 = reshape(abserror5m, [],1);
IC5 = ICcalculations(abserror5, 4, 100);
abserror6 = reshape(abserror6m, [],1);
IC6 = ICcalculations(abserror6, 3, 100);
abserror7 = reshape(abserror7m, [],1);
IC7 = ICcalculations(abserror7, 2, 100);

numterms = [15, 14, 13, 9, 4, 3, 2];
ICvec = [IC1.aic_c IC2.aic_c IC3.aic_c IC4.aic_c IC5.aic_c IC6.aic_c IC7.aic_c]-IC7.aic_c;
% these are the relative AICc valuse listed in Table 1 of the paper.
disp('relative AICc for models 1-7:')
ICvec'

rmpath('models/Burgers_Proc')
%%
if plottag ==1
figure(2)
plot(numterms, ICvec, 'o')

figure(2)

subplot(2,3,1)
waterfall(x,t,real(usol1.')); colormap([0 0 0]);
view(15,35), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2],'Ylim',[0 10],'Ytick',[0 10],'Zlim',[0 1],'Ztick',[0 0.5 1])

subplot(2,3,2)
waterfall(x,t,real(usol2.')); colormap([0 0 0]);
view(15,35), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2],'Ylim',[0 10],'Ytick',[0 10],'Zlim',[0 1],'Ztick',[0 0.5 1])

subplot(2,3,3)
waterfall(x,t,real(usol3.')); colormap([0 0 0]);
view(15,35), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2],'Ylim',[0 10],'Ytick',[0 10],'Zlim',[0 1],'Ztick',[0 0.5 1])

subplot(2,3,4)
waterfall(x,t,real(usol4.')); colormap([0 0 0]);
view(15,35), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2],'Ylim',[0 10],'Ytick',[0 10],'Zlim',[0 1],'Ztick',[0 0.5 1])

subplot(2,3,5)
waterfall(x,t,real(usol5.')); colormap([0 0 0]);
view(15,35), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2],'Ylim',[0 10],'Ytick',[0 10],'Zlim',[0 1],'Ztick',[0 0.5 1])

subplot(2,3,6)
waterfall(x,t,real(usol6.')); colormap([0 0 0]);
view(15,35), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2],'Ylim',[0 10],'Ytick',[0 10],'Zlim',[0 1],'Ztick',[0 0.5 1])

%
figure(3)
plot(x,real(usol(:,end))','o')
hold on
plot(x,real(usol1(:,end)), '-')
plot(x,real(usol2(:,end)),'-')
plot(x,real(usol3(:,end)),'-')
plot(x,real(usol4(:,end)),'-')
plot(x,real(usol5(:,end)),'-')
plot(x, real(usol6(:,end)),'-')
plot(x, real(usol7(:,end)),'-')
xlabel('x')
ylabel('u')
end

