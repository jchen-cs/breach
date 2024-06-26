function [val__, time_values__] = STL_EvalThom_Gen(Sys, phi, P, trajs, partition, relabs, t)
%STL_EVALTHOM_GEN computes the satisfaction function of a property for one
% or many trajectory(ies) with respect to a partition of signals. 
% This function uses a variable time step robust
% monitoring algorithm.
%
% Synopsis: [val__, time_values__] = STL_EvalThom_GEN(Sys, phi, P, trajs, partition, relabs, t)
%
% Input:
%  - Sys   : is the system
%  - phi   : is a STL property
%  - P     : is a parameter set which contains one parameter vector only
%            used for properties parameters.
%  - trajs : is a structure with fields X and time. It may contains many
%            trajectories. In this case, all will be checked considering
%            the property parameters described in P.
%  - partition  : is a set of signals given as an array of strings:
%                 the robustness computation is computed in these signals.
%  - relabs : is a string indicating how to treat variables that are 
%             not in the partition: 'rel' for -inf/+inf or 'abs' for 0.
%  - t     : (optional, default=traj.time) is the time point(s), so
%            possibly an array, when to eval the satisfaction of the
%            property. All time points not belonging to traj.time will be
%            linearly interpolated.
%
% Output:
%  - val__         : is a cell array of dimension 1 x numel(trajs). Each
%                    cell contains a line array describing the evaluation
%                    of phi at each time of time_values__. If numel(trajs)
%                    is 1, val__ is the content of its only cell to avoid a
%                    useless cell array (so, it is a line array).
%  - time_values__ : is a cell array of dimension 1 x numel(trajs). Each
%                    cell contains an array indicating time point(s) at
%                    which the formula is evaluated or interpolated. If the
%                    parameter t is provided, each cells of time_values__
%                    is equal to t. If numel(trajs) is 1, time_values__
%                    contains the content of the only cell (this avoid a
%                    useless cell array), and thus becomes a line array.
%
%See also STL_Eval SEvalProp
%

%
%  This works with piecewise affine signals.
%

%% defines the parameter as global variables so that they are available for
% all subsequent computations

global BreachGlobOpt;

if ~isempty(P.ParamList)
    BreachGlobOpt.GlobVarsDeclare = ['global ', sprintf('%s ',P.ParamList{:})]; % contains parameters and IC values (can remove IC if phi is optimized)
    eval(BreachGlobOpt.GlobVarsDeclare); % These values may be used in generic_predicate and GetValues
else
    BreachGlobOpt.GlobVarsDeclare = ''; % contains parameters and IC values (can remove IC if phi is optimized)
end

ii=1;
num_dim = size(P.pts,1); 
eval_str = [P.ParamList(1:num_dim);num2cell(1:num_dim)];
eval_str = sprintf('%s=P.pts(%d,ii);',eval_str{:});
eval(eval_str);

%% for each trajectory, compute values and times

numTrajs = numel(trajs);
val__ = cell(1, numTrajs);
time_values__ = cell(1, numTrajs);

if isstruct(trajs)||isa(trajs, 'matlab.io.MatFile')
    trajs = {trajs};
end

for ii=1:numTrajs % we loop on every traj in case we check more than one
    if (Psize_pts(P)==1)
        Pii = P;
    else
        Pii = Sselect(P, ii);
        eval(eval_str); % needed, as parameters can change from one Pii to another
    end
    
    % Ensures that traj.X and traj.time are double precision
    traj.time = double(trajs{ii}.time);
    traj.X = double(trajs{ii}.X);
    
    % Robusthom doesn't like singular intervals - should be optimized one
    % of these days ...
    if exist('t','var')
        time_values__{ii} = t;
        if(numel(t)==1)
            tn = find(traj.time>t,1);
            if isempty(tn)
                interval = [t t+1]; % Maybe not the best choice, but we have to make one !
            else
                interval = [t traj.time(1,tn)];
            end
        else
            interval = [t(1) t(end)];
        end
        [val, ~, ~, time_values] = GetValues(Sys, phi, Pii, traj, partition, relabs, interval);
        try
 %           if(numel(t)==1) % we handle singular times
 %               val__{ii} = val(1);
 %           else
                if isfield(BreachGlobOpt, 'disable_robust_linear_interpolation')&&BreachGlobOpt.disable_robust_linear_interpolation
                    val__{ii} = interp1(time_values, val, t, 'previous');
                else
                    val__{ii} = interp1(time_values, val, t);
                end
%            end
        catch % if val is empty
            val__{ii} = NaN(1,numel(t));
        end
    else
        interval = [0 traj.time(1,end)];
        [val__ii, ~, ~, time_values__ii] = GetValues(Sys, phi, Pii, traj, partition, relabs, interval);
        val__{ii} = val__ii(time_values__ii<=traj.time(1,end));
        time_values__{ii} = time_values__ii(time_values__ii<=traj.time(1,end));
        
        if (time_values__{ii}(end) < traj.time(end))
            vend__ = interp1(time_values__ii, val__ii,traj.time(1,end));
            time_values__{ii}(end+1) = traj.time(1,end);
            val__{ii}(end+1)= vend__;
        end
        
    end
end

if(numTrajs==1)
    val__ = val__{1};
    time_values__ = time_values__{1};
else
    if numel(t)==1
        val__ = cell2mat(val__);
        time_values__ = cell2mat(time_values__);
    end
end
end
%%

function [valarray, valarray_P, valarray_N, time_values] = GetValues(Sys, phi, P, traj, partition, relabs, interval)
global BreachGlobOpt;
eval(BreachGlobOpt.GlobVarsDeclare);

% Constants for TeLEx semantics
TeLEx_beta = 1;
TeLEx_gamma = 0.01;

switch(phi.type)
    case 'predicate'
        time_values = GetTimeValues(traj, interval); % Gets time steps of trace
        params = phi.params;
        params.Sys = Sys;
        params.P = P;
        evalfn = @(t) phi.evalfn(0, traj, t, params); % will call generic_predicate
        
        try   % works if predicate can handle multiple time values
            valarray = evalfn(time_values);
        catch %#ok<CTCH>
            valarray = arrayfun(evalfn, time_values);
        end
        % valarray = robustness score at each time step
        sigs = STL_ExtractSignals(phi); % sigs = signals involved in phi
        switch relabs
            case 'rel'
                % Check whether the predicate talks only about the signals
                % that are not in the partition. If yes, treat the predicate
                % qualitatively
                if ~isempty(sigs) % parameters only are always treated quantitatively
                    all_outside_partition = isempty(intersect(sigs, partition));
                    if(all_outside_partition)
                        valarray = Inf*(2*(valarray>0)-1);
                    end
                end
                
            case 'abs'
                % Checks whether the predicate has at least one signal that
                % is outside the partition. If yes, treat the predicate
                % qualitatively
                if ~isempty(sigs) % parameters only are always treated quantitatively
                    one_outside_partition = ~isempty(setdiff(sigs, partition));
                    if (one_outside_partition)
                        valarray = zeros(size(valarray));
                    end
                end
        end
        valarray_P = s_nu(time_values, max(0, valarray), phi.semantics);
        valarray_N = s_mu(time_values, min(0, valarray), phi.semantics);
        
    
        
    case 'not'
        [~, P, N, time_values] = GetValues(Sys, phi.phi, P, traj, partition, relabs, interval);
        valarray_P = -N;
        valarray_N = -P;
        
    case 'or'
        [~, valarray_P1, valarray_N1, time_values1] = GetValues(Sys, phi.phi1, P, traj, partition, relabs, interval);
        [~, valarray_P2, valarray_N2, time_values2] = GetValues(Sys, phi.phi2, P, traj, partition, relabs, interval);
        [t, P1, P2] = matchTraces(time_values1, time_values2, valarray_P1, valarray_P2);
        [~, N1, N2] = matchTraces(time_values1, time_values2, valarray_N1, valarray_N2);
        valarray_P = s_beta(P1, P2, phi.semantics);
        valarray_N = -s_alpha(-N1, -N2, phi.semantics);
        time_values = t;
        
        
        
    case 'and'
        [~, valarray_P1, valarray_N1, time_values1] = GetValues(Sys, phi.phi1, P, traj, partition, relabs, interval);
        [~, valarray_P2, valarray_N2, time_values2] = GetValues(Sys, phi.phi2, P, traj, partition, relabs, interval);
        [t, P1, P2] = matchTraces(time_values1, time_values2, valarray_P1, valarray_P2);
        [~, N1, N2] = matchTraces(time_values1, time_values2, valarray_N1, valarray_N2);
        valarray_P = s_alpha(P1, P2, phi.semantics);
        valarray_N = -s_beta(-N1, -N2, phi.semantics);
        time_values = t;
        
        
%     case 'andn'
%         n_phi = numel(phi.phin);
%         valarray_P = cell(1,n_phi);
%         valarray_N = cell(1,n_phi);
%         valarray = cell(1,n_phi);
%         time_values = cell(1,n_phi);
%         for ii=1:n_phi
%             [valarray_P{ii}, valarray_N{ii},time_values{ii}] = GetValues(Sys, phi.phin(ii), P, traj, partition, relabs, interval);
%             valarray{ii} = valarray_P{ii} + valarray_N{ii};
%         end
%         [time_values, valarray] = RobustAndn(time_values,valarray);
        
    case '=>'
        [~, valarray_P1, valarray_N1, time_values1] = GetValues(Sys, phi.phi1, P, traj, partition, relabs, interval);
        [~, valarray_P2, valarray_N2, time_values2] = GetValues(Sys, phi.phi2, P, traj, partition, relabs, interval);
        [t, P1, P2] = matchTraces(time_values1, time_values2, valarray_P1, valarray_P2);
        [~, N1, N2] = matchTraces(time_values1, time_values2, valarray_N1, valarray_N2);
        valarray_P = s_beta(-N1, P2, phi.semantics);
        valarray_N = -s_alpha(P1, -N2, phi.semantics);
        time_values = t;
        
        
    case 'always'
        % Using the equivalence G \phi = ~F(~\phi) = ~(T U ~\phi)
        I___ = eval(phi.interval);
        I___ = max([I___; 0 0]);
        I___(1) = min(I___(1), I___(2));
        % Interval of first predicate is from beginning of robustness
        % interval to the end, plus the Until operator's time horizon
        interval1 = [interval(1), I___(2)+interval(2)];
        interval2 = I___+interval;

        % True is infinite positive robustness, zero negative robustness
        valarray_P1 = inf(size(interval1));
        valarray_N1 = zeros(size(interval1));
        time_values1 = interval1;
        
        [~, valarray_P2, valarray_N2, time_values2] = GetValues(Sys, phi.phi, P, traj, partition, relabs, interval2); 
        [time_values, valarray_P_until, valarray_N_until] = PlusMinusUntil(time_values1, valarray_P1, valarray_N1, time_values2, -valarray_N2, -valarray_P2, interval, I___, phi.semantics);
        valarray_P = -valarray_N_until;
        valarray_N = -valarray_P_until;
       
%     case 'av_eventually'
%         I___ = eval(phi.interval);
%         I___ = max([I___; 0 0]);
%         I___(1) = min(I___(1), I___(2));
%         next_interval = I___+interval;
%         [valarray_P1, valarray_N1, time_values1] = GetValues(Sys, phi.phi, P, traj, partition, relabs, next_interval);
%         valarray1 = valarray_P1 + valarray_N1;
%         if(I___(end)~=inf)
%             time_values1 = [time_values1 time_values1(end)+I___(end)];
%             valarray1 = [valarray1 valarray1(end)];
%         end
%         [time_values, valarray] = RobustAvEvRight(time_values1, valarray1, I___);
        
    case 'eventually' 
        % Using the property that F \phi = T U \phi
        % Get time horizon for the Until operator
        I___ = eval(phi.interval);
        I___ = max([I___; 0 0]);
        I___(1) = min(I___(1), I___(2));
        % Interval of first predicate is from beginning of robustness
        % interval to the end, plus the Until operator's time horizon
        interval1 = [interval(1), I___(2)+interval(2)];
        interval2 = I___+interval;

        % True is infinite positive robustness, zero negative robustness
        valarray_P1 = inf(size(interval1));
        valarray_N1 = zeros(size(interval1));
        time_values1 = interval1;
        
        [~, valarray_P2, valarray_N2, time_values2] = GetValues(Sys, phi.phi, P, traj, partition, relabs, interval2); 
        [time_values, valarray_P, valarray_N] = PlusMinusUntil(time_values1, valarray_P1, valarray_N1, time_values2, valarray_P2, valarray_N2, interval, I___, phi.semantics);
        
        
%     case 'once'
%         I___ = eval(phi.interval);
%         I___ = max([I___; 0 0]);
%         I___(1) = min(I___(1), I___(2));
%                         
%         next_interval = interval-[I___(1)+I___(2), I___(1)];
%         next_interval(1) = max(0, next_interval(1));        
%         
%         [valarray_P1, valarray_N1, time_values1] = GetValues(Sys, phi.phi, P, traj, partition, relabs, next_interval);
%         valarray1 = valarray_P1 + valarray_N1;
%         % Flipping time, taking into account constant interpolation with
%         % previous 
%         Tend__ =  time_values1(end)+1; 
%         past_time_values1 = fliplr(Tend__-[time_values1 Tend__]);
%         past_valarray1 =    fliplr([valarray1(1) valarray1]);
%         
%         switch phi.semantics
%             case 'max'
%                 [past_time_values, past_valarray] = RobustEv(past_time_values1, past_valarray1, I___);  
%             case 'add'
%                 [past_time_values, past_valarray] = RobustAlways(past_time_values1, -past_valarray1, I___);  
%                 past_valarray = -past_valarray;
%             case 'vbool_v1'
%                 [past_time_values, past_valarray] = RobustAlways_v1(past_time_values1, -past_valarray1, I___);  
%                 past_valarray = -past_valarray;
%             case 'MARV'
%                 % On this level, MARV is just standard robustness, since
%                 % MARV only applies to top-level "always"-operator.
%                 [past_time_values, past_valarray] = RobustEv(past_time_values1, past_valarray1, I___);  
%             case 'constant'
%                 [past_time_values, past_valarray] = RobustEv(past_time_values1, past_valarray1, I___);  
%             otherwise
%                 error('Unknown objective function!');
%         end
%         
%         % Flipping back 
%         time_values = fliplr(Tend__-[past_time_values Tend__]);
%         valarray = fliplr([past_valarray(1) past_valarray]);   
%         
%     case 'historically'
%         I___ = eval(phi.interval);
%         I___ = max([I___; 0 0]);
%         I___(1) = min(I___(1), I___(2));
%         
%         next_interval = interval-[I___(1)+I___(2), I___(1)];
%         next_interval(1) = max(0, next_interval(1));        
%                 
%         [valarray_P1, valarray_N1, time_values1] = GetValues(Sys, phi.phi, P, traj, partition, relabs, next_interval);   
%         valarray1 = valarray_P1 + valarray_N1;
%         % Flipping time, taking into account constant interpolation with
%         % previous 
%         Tend__ =  time_values1(end)+1; 
%         past_time_values1 = fliplr(Tend__-[time_values1 Tend__]);
%         past_valarray1 =    fliplr([valarray1(1) valarray1]);  
%         
%         switch phi.semantics
%             case 'max'
%                 [past_time_values, past_valarray] = RobustEv(past_time_values1, -past_valarray1, I___);  
%             case 'add'
%                 [past_time_values, past_valarray] = RobustAlways(past_time_values1, past_valarray1, I___);  
%                 past_valarray = -past_valarray;
%             case 'vbool_v1'
%                 [past_time_values, past_valarray] = RobustAlways_v1(past_time_values1, past_valarray1, I___);  
%                 past_valarray = -past_valarray;
%             case 'MARV'
%                 % On this level, MARV is just standard robustness, since
%                 % MARV only applies to top-level "always"-operator.
%                 [past_time_values, past_valarray] = RobustEv(past_time_values1, -past_valarray1, I___);  
%             case 'constant'
%                 [past_time_values, past_valarray] = RobustEv(past_time_values1, -past_valarray1, I___);  
%             otherwise
%                 error('Unknown objective function!');
%         end
%         
%         % Flipping back
%         time_values = fliplr(Tend__-[past_time_values Tend__]);
%         valarray = -fliplr([past_valarray(1) past_valarray]);        

    case 'until'
        % Get time horizon for the Until operator
        I___ = eval(phi.interval);
        I___ = max([I___; 0 0]);
        I___(1) = min(I___(1), I___(2));
        % Interval of first predicate is from beginning of robustness
        % interval to the end, plus the Until operator's time horizon
        interval1 = [interval(1), I___(2)+interval(2)];
        interval2 = I___+interval;

        [~, valarray_P1, valarray_N1, time_values1] = GetValues(Sys, phi.phi1, P, traj, partition, relabs, interval1);
        [~, valarray_P2, valarray_N2, time_values2] = GetValues(Sys, phi.phi2, P, traj, partition, relabs, interval2);
        [time_values, valarray_P, valarray_N] = PlusMinusUntil(time_values1, valarray_P1, valarray_N1, time_values2, valarray_P2, valarray_N2, interval, I___, phi.semantics);
end

valarray = valarray_P + valarray_N;

%%  Sanity checks

% time progress
isblocked = (diff(time_values)<=0);
if find(isblocked)
    %  warning('Some time values not strictly increasing for property %s', disp(phi));
    keep = [1 find(~isblocked)+1];
    valarray = valarray(keep);
    time_values = time_values(keep);
end

ibof = isnan(time_values)|isinf(time_values);
if ~isempty(find(ibof, 1))
    val_ok = valarray(~ibof);
    time_ok = time_values(~ibof);
    valarray = interp1(time_ok, val_ok, time_values, 'nearest');
end

end

function time_values = GetTimeValues(traj,interval)
%GETTIMEVALUES provides time points belonging to traj.time strictly
% included in interval plus the bounds of interval.
%
% Note: for now it does not deal correctly with operations on time, e.g.
% x[t-a] of x[ g(t) ] where g is some function
%
% Also if we wanted to deal correctly with open or closed intervals that
% would be the place to look into. For now, the interpretation is rather
% that of closed intervals.
%

if(interval(1)==interval(end))
    time_values = interval(1);
    return ;
end

% first time instant
ind_ti = find(traj.time>=interval(1),1);
if isempty(ind_ti)
    if ~isempty(traj.time)
        time_values = [traj.time(1,end) traj.time(1,end)+1];
    else
        time_values = [];
    end
    return
end

if(traj.time(1,ind_ti)==interval(1))
    time_values = traj.time(1,ind_ti);
    ind_ti = ind_ti+1;
else
    time_values = interval(1);
end

% Last time instant
if(interval(end)==inf)
        time_values = [time_values traj.time(1,ind_ti:end)];
else
    ind_tf = find(traj.time >= interval(end),1);
    if isempty(ind_tf)
        time_values = [time_values traj.time(1,ind_ti:end) interval(end)];
    elseif(traj.time(1,ind_tf)==interval(end))
        time_values = [time_values traj.time(1,ind_ti:ind_tf)];
    else
        time_values = [time_values traj.time(1,ind_ti:ind_tf-1) interval(end)];
    end
end

end
