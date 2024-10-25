function C = compute_cost(sim,para,whichconstr,v,inc_or_prev,indirect)

% function to compute the cost of a model simulation using either constraint 
% INPUTS:
% sim:         results from model simulation
% alpha:       upper and lower bounds for cost per ICU induvidual
% para:        model distribution, capacity and cost parameters
% whichconstr: 'hard' or 'soft' to compute either constraint
% v:           scaling of logistic constraint
% inc_or_prev: flag to use different methods for calculation
% indirect:    whether to apply indirect costs

% simulation length
T = length(sim.new_ICU);

% survival distribution is the "normalised" time spent in ICU distribution
survival_dist = para.Time_C./sum(para.Time_C);
nts = length(survival_dist);

% background capacity (num patients in care without disease)
capacity = para.Ibar;
background = para.eta*para.Ibar;

% apply the hard or soft constraint
if isequal(whichconstr,'hard')
    
    % find where capacity is exceeded
    idxs_over_capacity = find(sim.in_ICU + background >= capacity);

    if inc_or_prev == 1

        % first term - smaller cost per hospitalisation
        C = para.alpha(1).*sim.new_ICU;

        % second term - higher cost per hospitalisation
        C(idxs_over_capacity) = para.alpha(2).*sim.new_ICU(idxs_over_capacity);

    else
    
        C = zeros(1,T+nts);
    
        % at each time point, the current incidence will be penalised (delta*1)
        % for current and future times below capacity and will be penalised
        % (delta*alpha) for current and future times above capacity, such that
        % this sums to the incidence case above
        for deltat = 1:T
            dt_lengthofstay = deltat:deltat+nts-1;
            
            times_over_capacity = ismember(dt_lengthofstay,idxs_over_capacity);
            dC_vec = sim.new_ICU(deltat).*(survival_dist.*(para.alpha(1).*ones(1,nts) + (para.alpha(2)-para.alpha(1)).*times_over_capacity));
            C(dt_lengthofstay) = C(dt_lengthofstay) + dC_vec;
        end
        C = C(1:T);
    
    end

    % append indirect costs
    if indirect
        Cindirect = max(0, min( sim.in_ICU+background-capacity, background ) );
        C = C + Cindirect;
    end

elseif isequal(whichconstr,'soft')

    %NB: v needs to be an NV x 1 vector for this to work - transpose if
    %need be
    if size(v,1) < size(v,2)
        v = v';
    end

    if inc_or_prev == 1
    
        % scaling of cost per new infection in time - logistic function
        alpha_t = para.alpha(1) + (para.alpha(2) - para.alpha(1))./(1 + exp(-v.*((sim.in_ICU + background)./capacity - 1)));
        C = alpha_t.*sim.new_ICU;
    
    elseif inc_or_prev == 0

        % append prevalence vector (ICU occupancy)
        Prev = [sim.in_ICU zeros(1,nts+1)];
    
        C = zeros(length(v),T+nts);
    
        % at each time point, the current incidence will be penalised (delta*1)
        % for current and future times below capacity and will be penalised
        % (delta*alpha) for current and future times above capacity, such that
        % this sums to the incidence case above
        for deltat = 1:length(sim.new_ICU)
            dt_lengthofstay = deltat+1:deltat+1+nts-1;

            alpha_t = para.alpha(1).*ones(length(v),nts) + (para.alpha(2)-para.alpha(1))./(1 + exp(-v.*((Prev(dt_lengthofstay)+background)./capacity - 1)));
            dC_vec = survival_dist.*alpha_t;
            C(:,dt_lengthofstay) = C(:,dt_lengthofstay) + dC_vec;
        end

        C = sim.new_ICU.*C(:,1:T);
    
    else

        % the faster way to compute this, but different (correct?) interpretation for
        % alpha

        % append incidence vector (ICU admissions)
        Prev = sim.in_ICU;
        Inc = [zeros(1,nts) sim.new_ICU];

        % time-dependent cost as a function of occupancy alpha_t
        alpha_t = para.alpha(1).*ones(length(v),length(Prev)) + (para.alpha(2)-para.alpha(1))./(1 + exp(-v.*((Prev+background)./capacity - 1)));

        % generate incidence matrix
        IncMat = zeros(T,nts); for j=1:T; IncMat(j,:) = Inc(j:j+nts-1); end

        % generate incidence times survival dist row vec
        DistMat = IncMat*flip(survival_dist)';

        C = alpha_t.*DistMat';

    end

    % append indirect costs
    if indirect
        logfunc = 1./(1 + exp(-v.*((sim.in_ICU+background)./capacity - 1)));
        %Cindirect = max( 0, min( (sim.in_ICU+background-capacity).*logfunc, background ) );
        Cindirect = max( 0, min( (sim.in_ICU+background-capacity), background ) );
        C = C + Cindirect;
    end

else

    disp('whichconstr flag must be "hard" or "soft"')
    C = 0;

end

