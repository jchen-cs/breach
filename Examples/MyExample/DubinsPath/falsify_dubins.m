clear; close all;

Sys = BreachSimulinkSystem('dubins');
Sys.PrintSignals();
SysFalsify = Sys.copy();

traj_gen = fixed_cp_signal_gen({'Velocity', 'TurnRate'});
InputGen = BreachSignalGen({traj_gen});
Sys.PrintParams();


SysFalsify.SetInputGen(InputGen);
SysFalsify.SetParamRanges({'Velocity_u0', 'TurnRate_u0'}, [0.1 3; 0.1 2]);


phi = STL_Formula('phi', 'alw(Dist2Target[t] > 0.01)');

% semantics = ["max", "add", "MARV", "constant"];
% results_iterations = nan(1, numel(semantics));
% results_time = nan(1, numel(semantics));
% for i=1:numel(semantics)
%     phi_test = set_semantics(phi, semantics(i));
%     req = BreachRequirement(phi_test);
%     falsify = FalsificationProblem(SysFalsify, req);
%     falsify.solve();
%     results_time(i) = falsify.time_spent;
%     results_iterations(i) = falsify.nb_obj_eval;
% end

req = BreachRequirement(phi);
falsify = FalsificationProblem(SysFalsify, req);
falsify.solve();
cex = falsify.GetFalse();
cex.PlotSignals({'Velocity', 'TurnRate', 'Dist2Target'});
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