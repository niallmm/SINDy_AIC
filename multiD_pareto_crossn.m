function [Lasso_CV] = multiD_pareto_crossn(Thetalib,crossval,lassovals, plottag)
% Run a multidimensional pareto front w/ cross validation.
% Theta is data matrix
% dx is data for derivatives in a Matrix of
% size (time series X n state variables)

Theta = Thetalib.Theta;
normTheta = Thetalib.normTheta;
dx = Thetalib.dx;
polyorder = Thetalib.polyorder;
usesine = Thetalib.usesine;

x0 = crossval.x0;
tvec = crossval.tvec;
xA = crossval.xA;
options = crossval.options;

numlasso = lassovals.numlasso;
lassostart = lassovals.lassostart;
lassoend = lassovals.lassoend;

LambdaLasso = logspace(lassostart,lassoend, numlasso);
% LambdaLasso = [logspace(-6, -2, numlasso) LambdaLasso1(2:end)];


[ntimeseries, nstates] = size(dx);
[ntimeseries, nfunc] = size(Theta);
[nstates,ninits] = size(x0);
% create lasso index


% declare stuff to be correct size
XiLasso_L2 = zeros(nfunc,numlasso,nstates);


% options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,nstates), 'Events', @events_lorenz);

% n is the size of the system

% Run LASSO independently on each state variable and do a Least-squares
% regression on the output for each lambda.
for nn = 1:nstates
    % run Lasso
%     [Xilasso(:,:,nn), FitInfo{nn}] = lasso(Theta, dx(:,nn), 'Lambda',LambdaLasso, 'Alpha', 1);
%    % perform Least squares on each output LASSO model.
%     for mm = 1:numlasso
%         Thetascreen = Theta*diag(abs(Xilasso(:,mm,nn))>1e-12);
%         XiLasso_L2(:,mm,nn) = Thetascreen\dx(:,nn)
%        XiLasso_L2(abs(XiLasso_L2(:,mm,nn))<1e-14,mm,nn) = 0
%     end
    for jj = 1:numlasso
        XiLasso_L2(:,jj,nn) = sparsifyDynamics(Theta,dx(:,nn), LambdaLasso(jj), 1);
    end
end


lasso_index = permn(1:numlasso,nstates);
Xiold= 0;
nn = 1;

Xilib= [];
% calcualate crossvalidated error
for ll = 1:length(lasso_index)

    lassonum_state = lasso_index(ll,:);
    % build the Xi matrix for this test (columns are coefficents for each
    % variable)
    Xibuild = [];
    for ii = 1:length(lassonum_state)
        if length(normTheta)>1
            Xibuild = [Xibuild XiLasso_L2(:,lassonum_state(ii),ii)./normTheta'];
        else
      Xibuild = [Xibuild XiLasso_L2(:,lassonum_state(ii),ii)];
        end
        Xibuild(abs(Xibuild)<1e-10) =0;
    end
    % put all coefficients into a vector to check if we have tried this Xi
    % before
    Xitest = reshape(Xibuild, nstates*nfunc,1);
    [junk, nlib] = size(Xilib);
    
    if (ll>1)&& any(all(abs((repmat(Xitest,1,nlib)-Xilib)<1e-10),1)) % only run for unique Xi
        continue
    else
        Xilib = [Xilib Xitest];
    end


    numberofcoefficients = nnz(Xibuild);
    for kk = 1:ninits
        
        % check if system is discrete or continuous.
        if length(tvec) == 1
%             disp('system is discrete')
            xtemp(:,1) = x0(:,kk);
            
            for jj = 2:tvec
                dx = sparseGalerkin(0,xtemp(:,jj-1),Xibuild,polyorder,usesine);
                xtemp(:,jj) = dx;
            end
            
            savetB{kk, nn} = 1:tvec-1;
            savexB{kk, nn} = xtemp(:,1:end-1)';
            
            tlength = min([length(savexB{kk,nn}) length(xA{kk})]);
            abserror(kk,nn) = sum(sum(abs(xA{kk}(1:tlength, :)-savexB{kk,nn}(1:tlength,:)))/tlength);
            RMSE(kk,nn) = sqrt(sum(sum(abs(xA{kk}(1:tlength, :)-savexB{kk,nn}(1:tlength,:)).^2)/tlength));
            
        else %if continuous
            sol=ode45(@(t,x)sparseGalerkin(t,x,Xibuild,polyorder,usesine),tvec,x0(:,kk),options);  % approximate
            
            savetB{kk,nn} = sol.x;
            savexB{kk,nn} = sol.y';
            
            tlength = min(length(sol.x), length(xA{kk}(:,1)));
            abserror(kk,nn) = sum(sum(abs(xA{kk}(1:tlength,:)-savexB{kk,nn}(1:tlength,:)))/tlength);
            RMSE(kk,nn) = sqrt(sum(sum(abs(xA{kk}(1:tlength,:)-savexB{kk,nn}(1:tlength,:)).^2)/tlength));
        end
        Xicomb{nn} = Xibuild;
        
%         if numberofcoefficients == 6
%             disp('pause for correct num coeff')
%             RMSE(kk,nn)
%         end

        if (plottag ==1)%% && ((RMSE(kk,nn)/abs(mean(mean(savexB{kk,nn}))))<0.5)
            Xibuild
            vlength = min(length(savetB{kk,nn}), length(xA{kk}));
            figure(111)
            plot(savetB{kk,nn}(1:vlength), xA{kk}(1:vlength,:), 'o')
            hold on
            vlength = min(length(savetB{kk,nn}), length(savexB{kk,nn}));
            plot(savetB{kk,nn}(1:vlength), savexB{kk,nn}(1:vlength,:))
            xlabel('time')
            ylabel('concentration')
            title('crossvalidation')
            drawnow
            hold off
        end
    end
    
    % calculate IC criteria for these errors
    
    numberoftrials = ninits
    
    [IC(nn)] = ICcalculations(abserror(:,nn),  numberofcoefficients , ninits);
    %abserror(:,mm,nn)
    IC(nn).aic_c;
    numcoeff(nn) = numberofcoefficients;
    nn = nn+1
    
    
end
Lasso_CV.tB = savetB;
Lasso_CV.xB = savexB;
Lasso_CV.abserror = abserror;
Lasso_CV.RSME = RMSE;
Lasso_CV.IC = IC;
Lasso_CV.numcoeff = numcoeff;
Lasso_CV.Xi = Xicomb;
Lasso_CV.lambda = LambdaLasso;
end


