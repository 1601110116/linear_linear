function h = viscose(g)
% This function calculates the sum of the diffusion terms of g
%  the resulting term h excludes boundaries
global viscosityxdt difzdt 
h = viscosityxdt .* (g(3:end, 2:end-1, 2:end-1) ...
	- 2*g(2:end-1, 2:end-1, 2:end-1) + g(1:end-2, 2:end-1, 2:end-1)) + ...
	viscosityxdt .* (g(2:end-1, 3:end, 2:end-1) ...
	- 2*g(2:end-1, 2:end-1, 2:end-1) + g(2:end-1, 1:end-2, 2:end-1)) + ...
	difzdt .* (g(2:end-1, 2:end-1, 3:end) - 2*g(2:end-1, 2:end-1, 2:end-1) ...
	+ g(2:end-1, 2:end-1, 1:end-2));
