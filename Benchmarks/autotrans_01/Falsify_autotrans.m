clear; close all;

filedatestring = datetime("now", 'format', 'yyyy-MM-dd_HHmmss-SSS');
diary autotrans.diary
diary on
fprintf("\n\n\n%s =====================================================\n", filedatestring)
Sys = BreachSimulinkSystem('sldemo_autotrans_mod01');
Sys.PrintParams();
Sys.PrintSignals();
SysFalsify = Sys.copy();

%throttle_gen = random_signal_gen({'Throttle'});
throttle_gen = fixed_cp_signal_gen({'Throttle'});
InputGen = BreachSignalGen({throttle_gen});

%InputGen.SetParam({'Throttle_min', 'Throttle_max', 'Throttle_dt_min', 'Throttle_dt_max', 'Throttle_seed'}, [0; 100; 30/7; 30/7; 42]);
SysFalsify.SetInputGen(InputGen);
SysFalsify.SetParamRanges({'Throttle_u0'}, [0 70]);

phi = STL_Formula('phi', 'alw_[0, 30] ((RPM[t] < 4500) & (Speed[t] < 120))');

% req = BreachRequirement(phi);
% 
% falsify = FalsificationProblem(SysFalsify, req);
% 
% falsify.solve();
% 
% 
% cex = falsify.GetFalse();
% cex.PlotSignals({'Throttle', 'Speed', 'RPM'});

%"max-breach", "const-breach", "plus-breach", "telex", "belta", "agm-product", "sum-product", "sum-min", "max-product", "minonly","smoothrect"
%semantics = ["max-breach", "const-breach", "plus-breach", "belta", "sum-min", "smoothrect", "smooth1", "smooth2"];
semantics = ["max-breach", "const-breach", "plus-breach", "telex", "belta", "sum-product", "sum-min", "max-product", "minonly", "smoothrect", "smooth1"];
results_iterations = nan(1, numel(semantics));
results_time = nan(1, numel(semantics));
for i=1:numel(semantics)
    semantics(i)
    phi_test = set_semantics(phi, semantics(i));
    req = BreachRequirement(phi_test);
    falsify = FalsificationProblem(SysFalsify, req);
    falsify.max_obj_eval = 500;
    falsify.solve();
    results_time(i) = falsify.time_spent;
    results_iterations(i) = falsify.nb_obj_eval;
end
diary off
writematrix([["semantics" "time" "iterations"];[semantics' results_time' results_iterations']], sprintf('autotrans_%s.csv', filedatestring));
