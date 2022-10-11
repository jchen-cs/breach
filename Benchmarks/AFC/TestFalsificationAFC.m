%% Initialization
% The following script creates a default interface with the
% AbstractFuelControl model.
clear; close all;

filedatestring = datetime("now", 'format', 'yyyy-MM-dd_HHmmss-SSS');
diary AFC.diary
diary on
fprintf("\n\n\n%s =====================================================\n", filedatestring)

BrDemo.InitAFC();
BrAFC

% % First we create the parameter search domain and load specifications. 
AFC_Falsify = BrAFC.copy();
AFC_Falsify.SetParamRanges({'Pedal_Angle_pulse_period', 'Pedal_Angle_pulse_amp'}, [10 20; 10 60]);
STL_ReadFile('AFC_simple_spec.stl');

%% Benchmark
phi = AF_alw_ok;
SysFalsify = AFC_Falsify;

%semantics = ["max", "add", "MARV", "constant", "TeLEx"];
semantics = ["max-breach", "const-breach", "plus-breach", "telex", "belta", "sum-product", "sum-min", "max-product", "minonly", "smoothrect", "smooth1"];
results_iterations = nan(1, numel(semantics));
results_time = nan(1, numel(semantics));
for i=1:numel(semantics)
    semantics(i)
    phi_test = set_semantics(phi, semantics(i));
    req = BreachRequirement(phi_test);
    %falsify = MaxSatProblem(SysFalsify, req);
    falsify = FalsificationProblem(SysFalsify, req);
    falsify.solve();
    results_time(i) = falsify.time_spent;
    results_iterations(i) = falsify.nb_obj_eval;
end
diary off
writematrix([["semantics" "time" "iterations"];[semantics' results_time' results_iterations']], sprintf('AFC_%s.csv', filedatestring));
