function f = myfun(x)
% global CONSTANTS
% f = -(CONSTANTS.Theta(1)+CONSTANTS.Theta(2)*x+CONSTANTS.Theta(3)*(x^2));

global PARAMETERS
f = -(PARAMETERS.VFCS(1)+PARAMETERS.VFCS(2)*x+PARAMETERS.VFCS(3)*(x^2));
