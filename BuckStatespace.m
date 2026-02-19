%%
function [] = xxx_buck_6()
%
% Buck switch mode power supply simulation in Matlab.
%
%{
This is a tidy version of the previous iterations to show the core of
the simulation loop. It includes peak current limiting and discontinious
mode detection.
State space is awseome for Matlab and time domain circuits, let alone the
usefulness of it for SMPS analysis. As a reminder for the future, state
space writes the differential equations in matrix form where each "state
varible" is inductor current or capacitor voltage, and the general
form is
	xdot = A.x + B.u
	y    = C.x + D.u
where u are the inputs, and x is a vector of the state variables.
For SMPS analysis, one gets the state space equation for ON, OFF, DCM etc. 
Then in Matlab we might as well use the built in ode solvers to solve this,
with a bit of fettling to switch between the operating states. This is
clearly better than discretizing the state space and solving at fixed time
steps.
[discritzitation uses bilinear approximation to exp(At)]
Papers often use state space averaging to combine the different modes into
a single state space matrix, which they then use to derive small signal
models etc.
Previous versions investigated different ways of dealing with discontinious
mode operation, and different ode solvers.
Now with more playing the simple mode switch is the right way, and use the
ode solver "events" function to detect DCM and peak current control.
%}
close all;
clc;
%
% Choose which ode solver to use
%
% ODESOLVER = @(t, y, varargin) ode45(t, y, varargin{:}); 
ODESOLVER = @(t, y, varargin) ode23(t, y, varargin{:}); 
%
% Input data for a simple fixed duty cycle simulation - i.e. this code is
% not a SMPS just simulation of the switches and L-C tank
%
Vin = 8;
L = 5e-6;
C = 100e-6;
Fs = 100e3;
% R = 1.0; iLimit = -Inf; % CCM, no peak current limiting
% R = 2.0; iLimit = 5; % CCM, peak current limiting
R = 10.0; iLimit = -Inf; % DCM, no peak current limit
R = 10.0; iLimit = 3; % DCM, peak current limit
D1 = [linspace(1e-3, 0.75, 50), linspace(0.75, 0.75, 100)];
cycles = numel(D1);
%%
ax = [];
figure(11);
ax.iL = subplot(2,1,1);
set(ax.iL, 'nextplot', 'add');
grid(ax.iL, 'on');
ax.vC = subplot(2,1,2);
set(ax.vC, 'nextplot', 'add');
grid(ax.vC, 'on');
%%
% State space matricies for ON and OFF and DCM operation
% x = [iL; vC]'
ON = [];
ON.A = [
	0		-1/L
	1/C		-1/(R*C)
	];
ON.B = [
	1/L
	0
	];
ON.C = [1 1];
ON.D = [
	0
	0
	];
OFF = [];
OFF.A = [
	0		-1/L
	1/C		-1/(R*C)
	];
OFF.B = [
	0
	0
	];
OFF.C = [1 1];
OFF.D = [
	0
	0
	];
DCM = [];
DCM.A = [
	0		0
	0		-1/(R*C)
	];
DCM.B = [
	0
	0
	];
DCM.C = [1 1];
DCM.D = [
	0
	0
	];
u = Vin;
ON.odefn = @(t,x) ON.A*x + ON.B*u;
OFF.odefn = @(t,x) OFF.A*x + OFF.B*u;
DCM.odefn = @(t,x) DCM.A*x + DCM.B*u;
%%
% Simplest version that does not bother with modelling a diode - i.e. a
% synchronous converter. In this case discontinious mode is not relevant
% and current can flow back from the inductor to the input source.
%
% For each PWM cycle, we use ODESOLVER to solve the state space eqns for the ON
% portion and OFF portion on the cycle.  All we need to do, is trim out
% duplicate timesteps at the switching break points.
%
% set up initial conditions and vectors for the results - as we don't know
% how many timesteps we need, we cannot preallocate easily
x = [0; 0];
t = [0];
n = [];
for cc = 1:cycles
	% ON portion
	[tode, xode] = ODESOLVER(ON.odefn, [0 D1(cc)/Fs], x(:,end));
	t = horzcat(t, tode(2:end)'+(cc-1)/Fs);
	x = horzcat(x, xode(2:end,:)');
	n(cc) = numel(tode);
	
	% OFF portion
	[tode, xode] = ODESOLVER(OFF.odefn, [D1(cc)/Fs 1/Fs], x(:,end));
	ii = x(:,1) < 0;
	t = horzcat(t, tode(2:end)'+(cc-1)/Fs);
	x = horzcat(x, xode(2:end,:)');
	n(cc) = n(cc) + numel(tode);
	
end;
iL = x(1,:)*ON.C(1);
vC = x(2,:)*ON.C(2);
plot(ax.iL, t, iL, 'o-r');
plot(ax.vC, t, vC, 'r');
linkaxes([ax.iL, ax.vC], 'x');
fprintf('mean iterations per pwm cycle = %.0f\n', mean(n));
%%
% Now using the 'event' capability of the ode solver to detect peak current and DCM. It
% really is very simple to add these into the model. Although it did take a
% few hours of work to whittle it down to this and I learnt more about the
% ode solver capabilities.
%
x = [0; 0];
t = [0];
n = [];
for cc = 1:cycles
	
	% ON portion
	% Will be terminated early if the event function triggers
	[tode, xode, tOFF] = ODESOLVER(ON.odefn, [0 D1(cc)/Fs], x(:,end), odeset('Events', @(t, y) iPeakEvent(t, y, iLimit)));
	t = horzcat(t, tode(2:end)'+(cc-1)/Fs);
	x = horzcat(x, xode(2:end,:)');
	n(cc) = numel(tode);
	
	% Ensure that we start the OFF phase at the end of the ON phase
	if (isempty(tOFF))
		tOFF = D1(cc)/Fs;
	end;
	
	% OFF portion
	% Will be terminated early if the event function triggers
	[tode, xode, tDCM] = ODESOLVER(OFF.odefn, [tOFF 1/Fs], x(:,end), odeset('Events', @dcmEvent));
	t = horzcat(t, tode(2:end)'+(cc-1)/Fs);
	x = horzcat(x, xode(2:end,:)');
	n(cc) = n(cc) + numel(tode);
		
	% DCM portion, if any
	if (~isempty(tDCM))
		[tode, xode] = ODESOLVER(DCM.odefn, [tDCM 1/Fs], x(:,end));
		n(cc) = n(cc) + numel(tode);
		t = horzcat(t, tode(2:end)'+(cc-1)/Fs);
		x = horzcat(x, xode(2:end,:)');		
	end;
	
end;
iL = x(1,:)*ON.C(1);
vC = x(2,:)*ON.C(2);
plot(ax.iL, t, iL, 'o-m', 'linewidth', 2.0);
plot(ax.vC, t, vC, 'm', 'linewidth', 2.0);
fprintf('mean iterations per pwm cycle = %.0f\n', mean(n));
end
%%
function [value,isterminal,direction] = iPeakEvent(t, y, iLimit)
% Locate the time when iL passes through zero in a decreasing direction and stop solver
value = y(1) - iLimit;     % Detect iL above a thresold
isterminal = 1;   % Stop the integration
direction = 1;   % Positive direction only 
end
%%
function [value,isterminal,direction] = dcmEvent(t, y)
% Locate the time when iL passes through zero in a decreasing direction and stop solver
value = y(1);     % Detect iL
isterminal = 1;   % Stop the integration
direction = -1;   % Negative direction only 
end
