function [abserror, RMSE, savetB, savexB] = ...
    validateXi(Xi, Thetalib,validation, plottag)
% perform a cross validation given Xi
% returns absolute error and root mean squared error for each test case
% also returns generated test time series data in data structures
%

polyorder = Thetalib.polyorder;
usesine = Thetalib.usesine;

x0 = validation.x0;
tA = validation.tA;
xA = validation.xA;
options = validation.options;

[nterms, ninits] = size(x0);
for kk = 1:ninits % for each initial condition
    % check if system is discrete or continuous.
    tvec = tA{kk};
    if length(tvec) == 1 % if discrete:
        
        % set initial condition
        xtemp(:,1) = x0(:,kk);
        
        % directly calculate each next time step
        for jj = 2:tvec
            dx = sparseGalerkin(0,xtemp(:,jj-1),Xi,polyorder,usesine);
            xtemp(:,jj) = dx;
        end
        
        % save results in structure
        savetB{kk} = 1:tvec-1;
        savexB{kk} = xtemp(:,1:end-1)';
        
        % save error
        tlength = min([length(savexB{kk}) length(xA{kk})]);
        abserror(kk) = sum(sum(abs(xA{kk}(1:tlength,:)-savexB{kk,nn}(1:tlength,:)))/tlength);
        RMSE(kk) = sqrt(sum(sum(abs(xA{kk}(1:tlength,:)-savexB{kk,nn}(1:tlength,:)).^2)/tlength));
        
    else %if continuous integrate using the recovered functional form
        sol=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),...
            tvec,x0(:,kk),options);  % approximate
        % if the solution diverges quickly:
        dydx = abs(sol.y(end)-sol.y(end-1)/(sol.x(end)-sol.x(end-1)));

            
        % check if the test model produced a time series longer than
        % the validation time-series or vice-versa
        ind_tlast = find(tvec<=(max(sol.x)), 1, 'last');
        ind_tlast2 = find(sol.x<=(max(tvec)), 1, 'last');
        
        if dydx > 1e10
            savetB{kk} = tvec;
            savexB{kk} = [sol.y(1); sol.y(end)*ones(length(tvec)-1,1)];
            ind_tlast = length(tvec);
        elseif ind_tlast<ind_tlast2 % if the validation time-series
            % is longer than the test model, evaluate the test model at
            % the validation time-series points, until the last time
            % point that is calculated by both
            savetB{kk} = tvec(1:ind_tlast);
            savexB{kk} = deval(sol, tvec(1:ind_tlast))';
            
        elseif ind_tlast>=ind_tlast2 % if the test model results are
            % longer than or equal to the validation time-series,
            % evaluate at validation points only
            savetB{kk} = tvec;
            savexB{kk} = deval(sol, tvec)';
        else
            disp('some weird case!!!')
        end
        % define vectors to compare, making sure they are the same
        % length, and calculate length
        xAcomp =xA{kk}(1:ind_tlast, :);
        xBcomp = savexB{kk};
        tlength = length(xBcomp);
        
        % calculate errors
        if length(xAcomp) == tlength
            abserror(kk) = sum(sum(abs(xAcomp-xBcomp))/tlength)
            RMSE(kk) = sqrt(sum(sum(abs(xAcomp-xBcomp).^2)/tlength));
        else
            error('model and validation time-series are not the same length')
        end
    end
 
    % plot the validation data and current model
    if (plottag ==1)
        Xi
        figure(111)
        plot(tvec(1:ind_tlast), xAcomp, 'o')
        hold on
        plot(tvec(1:ind_tlast), xBcomp)
        xlabel('time')
        ylabel('concentration')
        title('crossvalidation')
        drawnow
        hold off
    end
end



