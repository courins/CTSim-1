function [x, phis, rmsds] = reconPwlsLALMOs14( y, w, geom, beta, pfun, itnlim, delta, numos, img0, imgOpt  )
% Penalized-Weight-Least-Squares recosntruction for single energy sinograms
% using linearized augmented largragian + ordered subset
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
% Eqn. (33)
%
% 2012-13 Meng Wu at Stanford University

fprintf('Reconstructing with single energy data using PWLS LALM OS 14... ');
tRecon = tic;

if nargin < 9
    img0 = reconFBP( y, geom );
end

stopCrt         = 1e-5;
k               = coneTruncatedSlices( geom );

% use FBP to compute initial image
x = img0;
x = extendVoi( x, k );

% load operators for projection and regularization
[A, At, Aos, Atos, Os ] = loadPojectors( geom, numos );

rw = At(w);
rw = extendVoi( rw, k );
[R, S, T ]  = loadPenaltyOperator( pfun, delta, rw );

a = A( ones( size(x), 'single' ) );
precom = At( w .* a );
precom = extendVoi( precom, k );

fprintf('\n\tbeta  = %11.2e, rho  = %11.2e\n', beta, 1 );
fprintf('\titnlim = %10g', itnlim);
hdg1 = '   itn       x(0)          PHI(x)        b*R(x)        RMSD';
fprintf('\n%s'      , hdg1  );

phis = zeros(itnlim, 1);
rmsds = zeros(itnlim, 1);

g = zeros( size(x), 'single');
l = 0;
for itn = 1 : itnlim
    
    for isub = 1 : numos
        
        if l >= 1
            dLold = dL;
        end
        
        xold = x;
        
        % step 3: d^{(k+1)} = d^{(k)} - A x^{k+1} + u^{(k+1)}
        r = Aos( x, isub ) - Os( y, isub );
        dL = Atos( Os( w, isub) .* r, isub  ) ;
        
        if l >= 1
            e = sum( (g(:) - dL(:) ) .* ( dL(:) - dLold(:)) );
            if e > 0
                l = 0;
            end
        end
        
        if l == 0
            rho = 1;
        else
            rho = max( pi/(l+1) * sqrt( 1- (pi / 2 / (l+1)) ), 0.001);
        end
        l = l + 1;
        
        g = rho / (rho + 1) * dL + 1/( rho + 1 ) * g;
        
        % step 1: s^{(k+1 )} = rho * dL( x^{(k)} ) + ( 1 - rho) g^{(k)}
        s = rho * dL + ( 1 - rho ) * g ;
        s( isnan(s) ) = 0;
        s( isinf(s) ) = 0;
        
        % step 2: x^{(k+1)} = prox_{ rho^-1 t }( x^{k} - ( rho^-1 t ) s^{k+1})
        % TODO: implemented by FISTA
        
        for j = 1:5
            x = x - ( beta * S(x) + rho * precom .* ( x - xold ) + s  ) ./ ( beta * T(x) + rho * precom );
            x( x < 0) = 0;
            x( isnan(x) ) = 0;
            x( isinf(x) ) = 0;
        end
        
        if isub == 1
            phi = 0;
        end
        
        
        temp = Os(w, isub).* r.^2;
        phi = phi + sum( temp(:) ) / 2;
        
    end
        
    phi = phi + beta * R(x);
    if nargin < 10
        rmsd = sqrt( mean( (x(:) - xold(:)).^2 ) );
    else
        if ndims( x ) == 3
            b = x(:,:,end/2) - imgOpt(:,:,end/2);
        else
            b = x - imgOpt;
        end
        rmsd = sqrt( mean( b(:).^2 ) );
    end
    
    if itn > 2 && phi > phis(itn-1)
        l = 0;
        numos = round( numos / 2 );
        fprintf('\treduce ordered subsets to %i', numos );
        [~, ~, Aos, Atos, Os ] = loadPojectors( geom, numos );
    end
    
    
    phis(itn) = phi;
    rmsds(itn) = rmsd;
    
    prnt = 0;
    if itn   <= 10      , prnt = 1; end
    if itn   >= itnlim-3, prnt = 1; end
    if rem(itn,10) == 0  , prnt = 1; end
    
    if prnt
        fprintf('\n%6g %13.3e'  , itn  , x( round(end/2), round(end/2),ceil(end/2)));
        fprintf(' %13.3e %13.3e %13.3e', phi  , beta * R(x), rmsd)
    end
    
    if itn > 5 && rmsd < stopCrt
        break;
    end
    
end



tRecon = toc(tRecon);
fprintf('\nDone in %dmin %0ds.\n', floor(tRecon/60), round(mod(tRecon, 60)));



end

