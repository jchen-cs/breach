clear; close all;

filedatestring = datetime("now", 'format', 'yyyy-MM-dd_HHmmss-SSS');
diary autotrans.diary
diary on
fprintf("\n\n\n%s =====================================================\n", filedatestring)
Sys = BreachSimulinkSystem('Autotrans_shift');
Sys.PrintParams();
Sys.PrintSignals();
SysFalsify = Sys.copy();

% %throttle_gen = random_signal_gen({'Throttle'});
% throttle_gen = fixed_cp_signal_gen({'Throttle'});
% InputGen = BreachSignalGen({throttle_gen});
% 
% %InputGen.SetParam({'Throttle_min', 'Throttle_max', 'Throttle_dt_min', 'Throttle_dt_max', 'Throttle_seed'}, [0; 100; 30/7; 30/7; 42]);
% SysFalsify.SetInputGen(InputGen);
% SysFalsify.SetParamRanges({'Throttle_u0'}, [0 100]);


%sig_gen = var_step_signal_gen({'throttle', 'brake'}, 5);
%SysFalsify.SetInputGen(sig_gen);

% SysFalsify.SetParamRanges({'dt_u0', 'dt_u1', 'dt_u2', 'dt_u3'}, ...
%                   [.1 10  ;  .1 10;    0.1 10;    0.1 10]);
% SysFalsify.SetParamRanges({'throttle_u0','brake_u1', 'throttle_u2', 'brake_u3'}, ... 
%                   [0 100;        0 325;      0 100;         0 325]);


% This code from ARCH2020 RE
%STaLiRo transmission benchmark seems to use 7 control points
sig_gen = var_cp_signal_gen({'throttle', 'brake'}, [7 2], 'previous'); 
SysFalsify.SetInputGen(sig_gen);

for i = 0:6
    SysFalsify.SetParamRanges(['throttle_u' num2str(i)], [0 100]);
end

SysFalsify.SetParamRanges({'throttle_dt0','throttle_dt1', 'brake_dt0' }, [0 25]);
SysFalsify.SetParamRanges('brake_u0', [0 325]);
SysFalsify.SetParamRanges('brake_u1', [0 325]);


%phi1_20 = STL_Formula('phi', 'alw_[0, 30] ((RPM[t] < 4500) & (Speed[t] < 120))');


phi1_20 = STL_Formula('phi', 'ev_[0, 20] (RPM[t] >= 2000)');
phi1_30 = STL_Formula('phi', 'ev_[0, 30] (RPM[t] >= 2000)');
phi1_40 = STL_Formula('phi', 'ev_[0, 40] (RPM[t] >= 2000)');

phis = [phi1_20 phi1_30 phi1_40];
%phi = STL_Formula('phi', 'alw(not (gear[t] == 3 and speed[t] < 30))')

% req = BreachRequirement(phi);
% 
% falsify = FalsificationProblem(SysFalsify, req);
% 
% falsify.solve();
% 
% 
% cex = falsify.GetFalse();
% cex.PlotSignals({'Throttle', 'Speed', 'RPM'});

semantics = ["max-breach"];

num_repeat = 10;
results_iterations = nan(1, numel(semantics)*num_repeat*numel(phis));
results_time = nan(1, numel(semantics)*num_repeat*numel(phis));
results_phi = strings(1, numel(semantics)*num_repeat*numel(phis));
results_semantics = strings(1, numel(semantics)*num_repeat*numel(phis));
exp_count = 0;
for i=1:numel(semantics)
    semantics(i)
    for k=1:numel(phis)
        phi = phis(k);
        for j=1:num_repeat
            exp_count = exp_count + 1;
            fprintf("iteration %d of formula %s of semantics %s\n", j, display(phi), semantics(i));
            phi_test = set_semantics(phi, semantics(i));
            req = BreachRequirement(phi_test);
            falsify = FalsificationProblem(SysFalsify, req);
            
            % uncomment for Nelder-Mead
            falsify.solver_options.num_corners = 0;
            falsify.solver_options.num_quasi_rand_samples = 1;
        
            % % uncomment for Simulated Annealing
            % falsify.solver = "simulated_annealing";
            % falsify.solver_options.start_sample = [];
        
        
            falsify.max_obj_eval = 500;
            falsify.solve();
            results_time(exp_count) = falsify.time_spent;
            results_iterations(exp_count) = falsify.nb_obj_eval;
            results_phi(exp_count) = string(display(phi));
            results_semantics(exp_count) = string(semantics(i));
        end
    end
end
diary off
writematrix([["semantics" "formula" "time" "iterations"];[results_semantics' results_phi' results_time' results_iterations']], sprintf('autotrans_%s.csv', filedatestring));
