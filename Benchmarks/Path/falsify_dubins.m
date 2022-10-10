clear; close all;

Sys = BreachSimulinkSystem('dubins');
Sys.PrintSignals();
SysFalsify = Sys.copy();

velocity_gen = fixed_cp_signal_gen({'Velocity'});
turn_gen = cp_signal_gen({'TurnRate'}, 2, 'previous');
InputGen = BreachSignalGen({velocity_gen, turn_gen});
Sys.PrintParams();


SysFalsify.SetInputGen(InputGen);
SysFalsify.SetParam({'TurnRate_u0'}, [0]);
SysFalsify.SetParamRanges({'Velocity_u0', 'TurnRate_t1', 'TurnRate_u1'}, [0.1 3; 0.01 2; -2 2]);


phi = STL_Formula('phi', 'alw(Dist2Target[t] > 0.01)');

%semantics = ["max", "add", "MARV", "constant", "TeLEx"];
%semantics = ["max-breach", "const-breach", "plus-breach", "belta", "sum-min", "smoothrect"];
semantics = ["sum-product"];
results_iterations = nan(1, numel(semantics));
results_time = nan(1, numel(semantics));
for i=1:numel(semantics)
    phi_test = set_semantics(phi, semantics(i));
    req = BreachRequirement(phi_test);
    falsify = FalsificationProblem(SysFalsify, req);
    %falsify.max_obj_eval = 1000;
    falsify.solve();
    results_time(i) = falsify.time_spent;
    results_iterations(i) = falsify.nb_obj_eval;
end
writematrix([["semantics" "time" "iterations"];[semantics' results_time' results_iterations']], 'results.csv');

% 
% cex = falsify.GetFalse();
% cex.PlotSignals({'Velocity', 'TurnRate', 'Dist2Target'});
% fig = figure;
% plot(cex.GetSignalValues({'X'})', cex.GetSignalValues({'Y'})')
% axis equal
% circle = annotation('ellipse');
% circle.Parent = fig.CurrentAxes;
% circle.Position = [0.99 0.99 0.02 0.02];
% %% Plot robustness
% xs = 0.1:0.1:3;
% vals = zeros(1, numel(xs));
% 
% 
% phi_test = set_semantics(phi, "max");
% for i = 1:numel(xs)
%     SysFalsify.SetParam({'Velocity_u0', 'TurnRate_u0'}, [1; xs(i)]);
%     SysFalsify.Sim();
%     vals(i) = SysFalsify.CheckSpec(phi_test);
% end
% 
% figure
% plot(xs, vals, xs, zeros(1, numel(xs)))
% title('Robustness vs Turn Velocity');
% xlabel('Turn Velocity');
% ylabel('Robustness, max semantics');
% %cex = falsify.GetFalse();
% %cex.PlotSignals({'LaunchAngle', 'LaunchVelocity', 'X', 'Y'});