clear; close all;

filedatestring = datetime("now", 'format', 'yyyy-MM-dd_HHmmss-SSS');
diary dubins2.diary
diary on
fprintf("\n\n\n%s =====================================================\n", filedatestring)

Sys = BreachSimulinkSystem('path2');
STL_ReadFile('path2_spec.stl');
Sys.PrintSignals();
SysFalsify = Sys.copy();

signal_gen = cp_signal_gen({'v', 'omega'}, [5; 5], 'previous');
InputGen = BreachSignalGen({signal_gen});
Sys.PrintParams();


SysFalsify.SetInputGen(InputGen);
SysFalsify.SetParam({'v_t1', 'v_t2', 'v_t3', 'v_t4'}, [1 2 3 4]);
SysFalsify.SetParam({'omega_t1', 'omega_t2', 'omega_t3', 'omega_t4'}, [1 2 3 4]);
SysFalsify.SetParamRanges({'v_u0', 'v_u1', 'v_u2', 'v_u3', 'v_u4'}, [0 1]);
SysFalsify.SetParamRanges({'omega_u0', 'omega_u1', 'omega_u2', 'omega_u3', 'omega_u4'}, [-1 1]);


phi = no_collide;
global usegradient
usegradient = false;

%semantics = ["max", "add", "MARV", "constant", "TeLEx"];
%semantics = ["max-breach", "const-breach", "plus-breach", "belta", "sum-min", "smoothrect"];
%semantics = ["max-breach", "const-breach", "plus-breach", "telex", "belta", "sum-product", "sum-min", "max-product", "minonly", "smoothrect", "smooth1"];
semantics = ["max-breach"];
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
writematrix([["semantics" "time" "iterations"];[semantics' results_time' results_iterations']], sprintf('dubins2_%s.csv', filedatestring));

% 
cex = falsify.GetFalse();
cex.PlotSignals({'v', 'omega'});
position = cex.GetSignalValues({'x', 'y'})';
plot(position(:,1), position(:,2));
rectangle('Position', [3 3 3 3]);
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