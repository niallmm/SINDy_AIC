function [abserror, RMSE, savetB, savexB] = ...
crossvalidateXi(Xi, fun, tspan, x0, polyorder, usesine, options, plottag, ...
savetA, savexA)
% perform a cross validation given Xi
% returns absolute error and root mean squared error for each test case
% also returns generated test time series data in data structures
% 
[nterms, nmeasure] = size(x0);
for kk = 1:nmeasure
%     x0(:,kk)
    % integrate approximate solution from Xi
    [tB,xB]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan,x0(:,kk),options);  % approximate
    savetB{kk} = tB;
    savexB{kk} = xB;
    
    abserror(kk,:) = sum(abs(savexA{kk}(1:length(xB),1:nterms)-xB))/length(xB);
    RMSE(kk,:) = sqrt(sum(abs(savexA{kk}(1:length(xB),1:nterms)-xB).^2)/length(xB));
  if plottag ==1
figure(111)
plot(savetA{kk}, savexA{kk}, 'o')
hold on
plot(tB, xB)
xlabel('time')
ylabel('concentration')
title('crossvalidation')
drawnow 
hold off
end
 end

