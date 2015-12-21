clear all;
close all;

data = load( 'burn_04.mat');
imgKr = data.imgBHC;

data = load( 'air_05.mat');
imgAir = data.imgBHC;

clear data;



%% Get image pixel that are not gas

kernelSize = 5;
%imgSub = imgKr - imgAir;

% segmentation
solid = imgAir > 0.3;
solid = solid | imgKr > 0.3;

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

final = final | abs( imgSub ) > 0.02;

% for slices that have porous media below system resolution
 final( :, :, 380 : end ) = false;


%% Now let's compute average density for each slice
close all;

% bounding box with the
x = [120 260];
y = [120 260];

%final( 120:140,120:140,:) = true;

att_curve = zeros( 1, size( imgAir, 3));
att_curve_bhc = zeros( 1, size( imgAir, 3));

for i = 1 : length( att_curve )
    
    slice = imgSub(x(1):x(2), y(1):y(2),i);
    slice_bhc = imgSubBHC(x(1):x(2), y(1):y(2),i);
    
    valid =  ~final( x(1):x(2), y(1):y(2),i);
    
    att_curve(i) = mean( slice( valid(:) ) );
    att_curve_bhc(i) = mean( slice_bhc( valid(:) ) );
end

figure; plot( att_curve ); hold on;
plot( att_curve_bhc, 'r' );
xlabel 'slice #', ylabel 'attenuation';


figure;
slice = imgSub(:,end/2,:);
slice( final(:,end/2,:) ) = 0;
imagesc( squeeze( slice )' , [0 0.02] ); axis image, colorbar, colormap jet;


figure;
slice = imgSubBHC(:,end/2,:);
slice( final(:,end/2,:) ) = 0;
imagesc( squeeze( slice )' , [0 0.02] ); axis image, colorbar, colormap jet;
