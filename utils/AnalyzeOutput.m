
%% plots
% plot everything on 2-D plot with the total number of terms on the x axis,
% and the aic for all models of that size on the y-axis.
if plottag >=1
    if exist('t')
        figure(111)
        plot(t,x, 'o')
        xlabel('t')
        ylabel('x')
        title('time series measurements')
        
    else
        figure(111)
        plot(x/Ntot, 'o')
        xlabel('t')
        ylabel('x')
        title('time series measurements')

    end


    %%
    figure(2)
    plot(numcoeff, cell2mat({IC.aic_c}), 'o')
    xlabel('number of terms')
    ylabel('AICc')
    
    figure(3);
    plot(numcoeff, AIC_rel, 'o')
     axis([min(numcoeff) max(numcoeff) 0 max(AIC_rel)])
%    axis([min(numcoeff) 20 0 max(AIC_rel)])
    ylabel('relative AICc')
    patch([0 max(numcoeff)+1 max(numcoeff)+1 0], ...
        [10 10  max(AIC_rel) max(AIC_rel)], ...
        [0.2 0.2 0.2], 'FaceAlpha', 0.6)
    patch([0 max(numcoeff)+1 max(numcoeff)+1 0], [4 4 7 7], [0.2 0.2 0.2], 'FaceAlpha', 0.3)
    patch([0 max(numcoeff)+1 max(numcoeff)+1 0], [0 0 2 2], [0.2 0.2 0.2], 'FaceAlpha', 0.1)
    

end

%% Analyze results
minind = find(min(cell2mat({IC.aic_c})) == cell2mat({IC.aic_c}), 1, 'first');
x1coeff = Xicomb{minind}
mincoeff = numcoeff(minind);

% find models with index below or equal to relAICc = 7;
indless7 = find(AIC_rel<=7);
[AIC_rel_uni, ia, ic] = unique(AIC_rel);
for kk = 1:max(numcoeff)
    kk
    % find indices of terms with kk coefficients
    ind = find(numcoeff==kk);
    % find unique indices in AIC_rel
    
    % find indices of unique AIC_rel valuse that have kk terms and
    % relAICc<=7
    indtog = intersect(intersect(indless7, ind), ia);
    
    if length(indtog)>0
        %     AIC_rel(indtog)
        %  figure(5)
        %     plot(kk, AIC_rel(indtog), 'og')
        Xiless7= Xicomb{indtog}
    end
end