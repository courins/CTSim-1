function imgOut = bilateralFilterImageDomain( imgIn, w, sd, ss )

fprintf('Bilateral filter in image domain with sd = %2.2f, ss = %2.2f ... \n', sd, ss);

imgOut = zeros( size(imgIn), class(imgIn) );

for iz = 1 : size( imgIn, 3 )
    
    fprintf('Bilateral filter in image slice %2.2f \n', iz);
    
    view = squeeze( imgIn(:,:,iz) );
    %imgOut(:,:,iz) = BilateralFilter( view, w, sd, ss );
    h = fspecial('gaussian', [8 8], 2);
    imgOut(:,:,iz) = imfilter(view,h,'replicate');
    
    
    
end

fprintf('\t done.\n');

end

