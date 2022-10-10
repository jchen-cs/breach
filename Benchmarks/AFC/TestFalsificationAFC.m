%% Initialization
% The following script creates a default interface with the
% AbstractFuelControl model.
clear; close all;
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
semantics = ["smooth1"];
results_iterations = nan(1, numel(semantics));
results_time = nan(1, numel(semantics));
for i=1:numel(semantics)
    phi_test = set_semantics(phi, semantics(i));
    req = BreachRequirement(phi_test);
    falsify = MaxSatProblem(SysFalsify, req);
    %falsify = FalsificationProblem(SysFalsify, req);
    falsify.solve();
    results_time(i) = falsify.time_spent;
    results_iterations(i) = falsify.nb_obj_eval;
end
