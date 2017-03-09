function [bool] = approx(A,B,tolerance)
% Roughly equal function for tests

bool = ( A < B + tolerance && A > B - tolerance );

end

