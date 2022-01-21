clear; close all;

Sys = BreachSimulinkSystem('cannonball');
Sys.PrintParams();
Sys.PrintSignals();
SysFalsify = Sys.copy();

angle_gen = fixed_cp_signal_gen({'LaunchAngle', 'LaunchVelocity'});
InputGen = BreachSignalGen({angle_gen});

SysFalsify.SetInputGen(InputGen);
SysFalsify.SetParamRanges({'LaunchAngle_u0', 'LaunchVelocity_u0'}, [0 90; 0 50]);
%SysFalsify.SetParamRanges({'LaunchAngle_u0'}, [0 90]);
%SysFalsify.SetParam({'LaunchVelocity_u0'}, [20]);
%phi = STL_Formula('phi', 'alw(((Y[t] < 10) and (Y[t] > 0)) => ((X[t] > 45) or (X[t] < 35)))');

phi = STL_Formula('phi', 'alw(((Y[t] < 0.5) and (Y[t] > 0)) => ((X[t] > 40.5) or (X[t] < 39.5)))');


% %semantics = ["max", "add", "MARV", "constant", "TeLEx"];
% semantics = ["TeLEx"];
% results_iterations = nan(1, numel(semantics));
% results_time = nan(1, numel(semantics));
% for i=1:numel(semantics)
%     phi_test = set_semantics(phi, semantics(i));
%     req = BreachRequirement(phi_test);
%     falsify = FalsificationProblem(SysFalsify, req);
%     falsify.max_obj_eval = 1000;
%     falsify.solve();
%     results_time(i) = falsify.time_spent;
%     results_iterations(i) = falsify.nb_obj_eval;
% end


%% Plot robustness
xs = 0:5:90;
vals = zeros(1, numel(xs));
phi = STL_Formula('phi', 'alw(((Y[t] < 10) and (Y[t] > 0)) => ((X[t] > 40.15) or (X[t] < 39.85)))');


phi_test = set_semantics(phi, "TeLEx");
for i = 1:numel(xs)
    SysFalsify.SetParam({'LaunchAngle_u0', 'LaunchVelocity_u0'}, [xs(i); 20]);
    SysFalsify.Sim();
    vals(i) = SysFalsify.CheckSpec(phi_test);
end

figure
plot(xs, vals, xs, zeros(1, numel(xs)))
title('Robustness vs Angle');
xlabel('Angle (°)');
ylabel('Robustness, max semantics');
%cex = falsify.GetFalse();
%cex.PlotSignals({'LaunchAngle', 'LaunchVelocity', 'X', 'Y'});