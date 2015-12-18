function y = huberPenaltyAdaptive( x, mode, delta, weights )
% compute huber penalty function cost value, derivative, and curvature
%   input:
%       x       - image to be penalized 2D or 3D in ([nx ny], [nx ny nz])
%       mode    - ( 0 cost function )
%                 ( 1 derivative )
%                 ( 2 curvature )
%       delta   - parameter for huber function (default 0.01)
%  output:
%       y       - results
%
% Meng Wu at Stanford University
% 2012 - 2013

y = huberPenalty( x, mode, delta ) .* weights;