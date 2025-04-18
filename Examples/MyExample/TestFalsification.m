% Simple test program that just checks a spec
clear;
mdl = 'LinearInputSystem'; % A system with a linear function

Sys = BreachSimulinkSystem(mdl); 
Sys.PrintSignals(); % Should get us an output named "ModelOut"

SysFalsify = Sys.copy();

input_gen      = fixed_cp_signal_gen({'ModelIn'}); % Just a constant signal
        
InputGen = BreachSignalGen({input_gen});

InputGen.SetParam({'ModelIn_u0'}, [-0.01]);

SysFalsify.SetInputGen(InputGen);

AboveZero = STL_Formula('AboveZero', 'alw(ModelOut[t] > 0)');

Time = 0:0.05:10;

AboveZero = set_semantics(AboveZero, 'max-breach')
SysFalsify.SetParamRanges({'ModelIn_u0'}, [-0.2 1]);
req = BreachRequirement(AboveZero);

falsify = FalsificationProblem(SysFalsify, req);
falsify.solver = 'simulated_annealing'
falsify.solve();

cex = falsify.GetFalse();
cex.PlotSignals({'ModelIn', 'ModelOut'});