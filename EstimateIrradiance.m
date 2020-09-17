fiberRadius = 0.1;
fiberNA = 0.22;
refractiveIndex = 1.36; %(grey matter)
rho =  fiberRadius * sqrt((refractiveIndex/fiberNA)^2 -1);
scatterCoeff = 11.2; % 11.2mm-1 for mouse and 10.3 mm-1 for rats
distance = 0.5; %500um

P0 = 3.96; % power output from fiber tip, here measured at 10mW stim
I0 =  round(P0/(pi*0.1^2),0);

% Irradiance at target distance (mW/mm^2)
Itarget = I0 * (rho^2 / ((scatterCoeff*distance) * (distance + rho)^2));