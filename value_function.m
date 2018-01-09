function [value]=value_function(x, theta, option)

%Given an inventory level, x, calculate the value for that level.
value = option * (theta( 1) + theta( 2)*x + theta( 3)*( x^2));

end


