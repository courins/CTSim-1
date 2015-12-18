function [x, phis, rmsds] = reconPwlsADMMOs( y, w, geom, beta, pfun, itnlim, delta, numos, img0, imgOpt   )
% Penalized-Weight-Least-Squares recosntruction for single energy sinograms
% using Nonuniform SQS algorithm
% input:
%       y       - log sinogram
%       w       - weights for pwls
%       geom    - system geometry
%       beta    - weight on penalty function
%       pfun    - roughness penalty function (default 'huber')
%       itnlim  - maximum number of iterations
%       delta   - parameter for penalty function
%       numos   - number of ordered subsets
%       img0    - initial image
%       imgOpt  - converged image
% output:
%       x       - reconstruction result
%       phis    - cost function values
%       rmsds   - root means squared differences with previous iteration or converged solution 
%
% Based Nien, Hung, and Jeffrey A. Fessler. �Fast X-Ray CT Image
% Reconstruction Using the Linearized Augmented Lagrangian Method with
% Ordered Subsets.� Optimization and Control; Learning; Machine Learning.
% arXiv Preprint arXiv:1402.4381 (February 18, 2014): 21.
%
% 2012-13 Meng Wu at Stanford University

fprintf('Reconstructing with single energy data using PWLS ADMM ... \n');
tRecon = tic;


if nargin < 4
    beta = 1e3;
end

if nargin < 5
    pfun = 'huber';
end

if nargin < 6
    itnlim = 200;
end

if nargin < 7
    delta           = 0.001;
end

if nargin < 8
    numos           = 24;
end

% use FBP to compute initial image
if nargin < 9
    img0 = reconFBP( y, geom );
end

stopCrt         = 1e-5;
k               = coneTruncatedSlices( geom );

rho = 0.1;

% load operators for projection and regularization
[A, At ] = loadPojectors( geom, numos );

rw = At(w);
rw = extendVoi( rw, k );
[R, S, T ]  = loadPenaltyOperator( pfun, delta, rw );

a = A( ones( size(x), 'single' ) );
precom = At( w .* a );
precom = extendVoi( precom, k );

fprintf('\nbeta  = %11.2e, rho  = %11.2e\n', beta, rho );
fprintf('\titnlim = %10g', itnlim);
hdg1 = '   itn       x(0)          PHI(x)        b*R(x)        RMSD';
fprintf('\n%s'      , hdg1  );

phis = zeros(itnlim, 1);
rmsds = zeros(itnlim, 1);

d = zeros( size(y), 'single');
u = A( x );
t = 1 ./ max( precom(:) );

for itn = 1 : itnlim
    
    xold = x;
    
    % step 2: u^{(k+1)} = argmin{ L( u )  + rho / 2 ||  A x^{k+1} - u - d^{k} ||^2  }
    p = A( x );
    u = u - ( w .* ( u - y ) + rho * ( u + d - p ) ) ./ ( w + rho );
    u( isnan(u) ) = 0;
    u( isinf(u) ) = 0;
    
    % step 3: d^{(k+1)} = d^{(k)} - A x^{k+1} + u^{(k+1)}
    d = d - p + u;
    
    
    % step 1: x^{(k+1 )} = argmin{ R(x) + rho / 2 ||  Ax - u^{k} - d^{k} ||^2  }
    
    s =  t / rho * At( w.* ( p - d - u )) ;
    for j = 1:10
        x = x - ( beta * S(x) + rho / t * ( x - xold +  s )  ) ./ ( beta * T(x) + rho / t  );
        
        x( x < 0 ) = 0;
        x( isnan(x) ) = 0;
        x( isinf(x) ) = 0;
        
    end
    
    
     phi = sum(  w(:).* (p(:) - y(:)).^2 ) / 2 + beta * R(x);
    if nargin < 10
        rmsd = sqrt( mean( (x(:) - xold(:)).^2 ) );
    else
        if ndims( x ) == 3
            a = x(:,:,end/2);
            b = imgOpt(:,:,end/2);
        else
            a = x;
            b = imgOpt;
        end
        rmsd = sqrt( mean( (a(:) - b(:)).^2 ) );
    end
    
    phis(itn) = phi;
    rmsds(itn) = rmsd;
    
    prnt = 0;
    if itn   <= 10       , prnt = 1; end
    if itn   >= itnlim-10, prnt = 1; end
    if rem(itn,10) == 0  , prnt = 1; end
    
    if prnt
        fprintf('\n%6g %13.3e'  , itn  , x(end/2,end/2,ceil(end/2)));
        fprintf(' %13.3e %13.3e %13.3e', phi  , beta * R(x), rmsd)
    end
    
    if itn > 5 && rmsd < stopCrt
        break;
    end
    
end

tRecon = toc(tRecon);
fprintf('\nDone in %dmin %0ds.\n\n\n', floor(tRecon/60), round(mod(tRecon, 60)));




end

