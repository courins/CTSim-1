function final = segmentPorousMediaBurner( imgAir, imgKr, kernelSize )

%kernelSize = 5;
%imgSub = imgKr - imgAir;

% segmentation
solid = imgAir > 0.3;
solid = solid | imgKr > 0.3;

% for slices that have porous media below system resolution
solid( :, :, 380 : end ) = false;
solid = solid |imgAir > 0.5;
solid = solid | imgKr > 0.5;

% mophological blurring in all dimension
final = solid;
solid_blurred = solid; 

% in x direction
for i = 1 : size( imgAir, 1 )
    solid_blurred( i, :, : ) = imdilate( solid(i,:,:) , ones(kernelSize) ); 
end
final = final | solid_blurred;

% in y direction
for i = 1 : size( imgAir, 1 )
    solid_blurred( :, i, : ) = imdilate( solid(:,i,:) , ones(kernelSize) ); 
end
final = final | solid_blurred;

% in z direction
for i = 1 : size( imgAir, 3 )
    solid_blurred( :, :, i ) = imdilate( solid(:,:,i) , ones(kernelSize) ); 
end
final = final | solid_blurred;


end