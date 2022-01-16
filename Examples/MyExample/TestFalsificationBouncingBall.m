% Simple test program that just checks a spec
clear;
mdl = 'BouncingBall'; % A system with a linear function

Sys = BreachSimulinkSystem(mdl); 
Sys.PrintSignals(); % Should get us an output named "ModelOut"

SysFalsify = Sys.copy();

input_gen      = fixed_cp_signal_gen({'InitialPosition'}); % Just a constant signal
        
InputGen = BreachSignalGen({input_gen});

InputGen.SetParam({'InitialPosition_u0'}, [0]);

SysFalsify.SetInputGen(InputGen);

SpeedLimit = STL_Formula('VelocityBelow20', 'alw((Velocity[t] < 20) & (Velocity[t] > -20))');

Time = 0:0.05:10;
global iterations;
iterations = 0;

SysFalsify.SetParamRanges({'InitialPosition_u0'}, [0 10]);
req = BreachRequirement(SpeedLimit);

falsify = FalsificationProblem(SysFalsify, req);
%falsify.solver = 'simulated_annealing'
%falsify.startSample = []
falsify.solve();
fprintf('%d iterations\n', round(iterations))

cex = falsify.GetFalse();
cex.PlotSignals({'Position', 'Velocity'});