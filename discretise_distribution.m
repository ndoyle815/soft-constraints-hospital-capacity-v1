function f = discretise_distribution(x,disttype,mean,sd)
% function generates a discretised distribution analogous to its continuous
% distribution counterpart using the method of Cori et al, 2013

% INPUTS:
% x: range of values over which to compute the PMF ("support")
% disttype: 'Gamma' or 'Normal'
% mean: mean of continuous distribution
% sd: standard deviation of continuous distribution

% initalise
Nx = length(x);
f = zeros(1,Nx);

% obtain shape and scale parameters (relevant for Gamma only)
gam_scale = sd^2/mean;
gam_shape = mean/gam_scale;

% main iteration - iterating for all x
for k = 1:Nx

    % define linspace over which to integrate
    deltax = x(2) - x(1);
    if k > 1 && k < Nx
        intValsStart = x(k-1);
        intValsEnd = x(k+1);
    elseif k == 1
        intValsStart = x(1) - deltax;
        intValsEnd = x(2);
    else
        intValsStart = x(k-1);
        intValsEnd = x(k) + deltax;
    end
    intVals = linspace(intValsStart, intValsEnd, 10001);  % assume x's evenly spaced

    % define function over which to integrate
    if isequal(disttype,'Gamma')
        funcVals = (deltax - abs(intVals - x(k))).*gampdf(intVals, gam_shape, gam_scale);

    elseif isequal(disttype,'Normal')
        funcVals = (deltax - abs(intVals - x(k))).*normpdf(intVals, mean, sd);

    else
        % invalid input
        disp('Distribution must be specified as Gamma or Normal')
        break

    end

    f(k) = trapz(intVals, funcVals);  % integrate
end

% normalise if necessary
f = f./sum(f);
