classdef alw_monitor < stl_monitor
    properties
        interval
        subphi
    end
    
    methods
        function this = alw_monitor(formula)
            this = this@stl_monitor(formula);
            this.subphi = get_children(this.formula);
            this.subphi = this.subphi{1};
            this.signals = {[this.formula_id '_violation'],... 
                            [this.formula_id '_quant_sat']};
                
            this.interval = get_interval(formula);
            this.init_P();
        end
        
        function [time, Xout] = computeSignals(this, time, X, p)
            this.init_tXp(time,X,p);

            % compute robustnes of top formula
            
            idx  = this.get_time_idx_interval(time,p);
            
            [~ , rob] = this.get_standard_rob(this.subphi, time(idx));            
            Xout = nan(2, numel(time));
            Xout(1, idx) = rob<0;
            Xout(2,idx) = rob;            
   
        end
        
        function [v, t, Xout] = eval(this, t, X,p)
            TeLEx_gamma = 0.01;
            [t, Xout] = this.computeSignals(t, X,p);
            idx  = this.get_time_idx_interval(t,p);
            Xout(end-1,:) = Xout(end,:)<0;         % violation flags
            Xout(end-1:end, ~idx) = NaN;
            
            % Xout(end, idx) array of obj function values inside alw()
            % t(idx) array of corresponding times
            
            % We want to use the same implementation as in
            % @STL_Formula/private/RobustAlways.m. However, since it's a
            % private function, we cannot reach it. Therefore, we have
            % copies called
            % alw_monitor_RobustAlways.m
            % alw_monitor_RobustAlways_v1.m
            % These function are in the Core/OutputGen folder. 
            time_values = t(idx);
            valarray = Xout(end,idx);
            valarray_P2 = max(valarray, 0);
            valarray_N2 = min(valarray, 0);
            I___ = [time_values(1) time_values(end)];

            % G phi = ~(T U ~phi)
            % rho_+(G phi) = -rho_-(T U ~phi)
            % rho_-(G phi) = -rho_+(T U ~phi)

            % Trace 1 will be just True
            valarray_P1 = inf(size(I___));
            %valarray_P1 = 500 * ones(size(I___));
            valarray_N1 = zeros(size(I___));
            time_values1 = I___;
            semantics = get_semantics(this.formula);

            [~, valarray_N_until, valarray_P_until] = PlusMinusUntil(time_values1, valarray_P1, valarray_N1, time_values, -valarray_N2, -valarray_P2, I___, I___, semantics);
            valarray_N_until = -valarray_N_until;
            valarray_P_until = -valarray_P_until;
            v = valarray_P_until + valarray_N_until;
            v = v(1);
            
           

        end
        
        function varargout = disp(this)
            phi = this.subphi;
            st = sprintf(['%s := alw_%s (%s)\n'], this.formula_id,this.interval,  disp(phi,1));
            
            if ~strcmp(get_type(phi),'predicate')
                st_pred = [];
                preds = STL_ExtractPredicates(phi);
                for ip = 1:numel(preds)
                    id = get_id(preds(ip));
                    status(ip)= STL_CheckID(id);
                end
                
                if any(status==1)
                    st_pred = '  where \n';
                    for ip = 1:numel(preds)
                        if status(ip)==1
                            st_pred =   sprintf([ st_pred '%s := %s \n' ],id,disp(preds(ip)));
                        end
                    end
                end
                st = [st st_pred];
            end
     
            if nargout == 0
                varargout = {};
                fprintf(st);
            else
                varargout{1} = st;
            end            
     
        end
        
        function plot_diagnostics(this, F)
            % Assumes F has data about this formula 
            
            F.AddSignals(this.signals_in);
            sig= this.signals{end};
            F.HighlightFalse(sig);
        
        end
    
    end
    
    methods  (Access=protected)
        function idx__ = get_time_idx_interval(this, t__, p__)
            this.assign_params(p__);
            interval = eval(this.interval);
            idx__ = (t__>= interval(1))&(t__<=interval(end));
            if ~any(idx__)
                idx__(end) = 1;
            end
        end
        
        
    end
end