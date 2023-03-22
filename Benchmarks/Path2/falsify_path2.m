clear; close all;

filedatestring = datetime("now", 'format', 'yyyy-MM-dd_HHmmss-SSS');
diary dubins2.diary
diary on
fprintf("\n\n\n%s =====================================================\n", filedatestring)

Sys = BreachSimulinkSystem('path2');
STL_ReadFile('path2_spec.stl');
Sys.PrintSignals();
SysFalsify = Sys.copy();



num_pts = 12;
dt = 0.1;
v_max = 1;
v_min = 0;
omega_max = 1;
omega_min = -1;

signal_gen = cp_signal_gen({'v', 'omega'}, [num_pts; num_pts], 'previous');
InputGen = BreachSignalGen({signal_gen});
Sys.PrintParams();

SysFalsify.SetInputGen(InputGen);

for i=1:(num_pts-1)
    SysFalsify.SetParam(['v_t' num2str(i)], dt*i);
    SysFalsify.SetParam(['omega_t' num2str(i)], dt*i);
end

for i=0:(num_pts-1)
    SysFalsify.SetParamRanges(['v_u' num2str(i)], [v_min v_max]);
    SysFalsify.SetParamRanges(['omega_u' num2str(i)], [omega_min omega_max]);
end

SysFalsify.PrintParams();


phi = no_collide;
%global usegradient
%usegradient = false;

%semantics = ["max", "add", "MARV", "constant", "TeLEx"];
%semantics = ["max-breach", "const-breach", "plus-breach", "belta", "sum-min", "smoothrect"];
%semantics = ["max-breach", "const-breach", "plus-breach", "telex", "belta", "sum-product", "sum-min", "max-product", "minonly", "smoothrect", "smooth1"];
semantics = ["smooth2"];
results_iterations = nan(1, numel(semantics));
results_time = nan(1, numel(semantics));
for i=1:numel(semantics)
    semantics(i)
    phi_test = set_semantics(phi, semantics(i));
    req = BreachRequirement(phi_test);
    falsify = FalsificationProblem(SysFalsify, req);
    falsify.max_obj_eval = 1000;
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