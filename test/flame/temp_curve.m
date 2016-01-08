%% Get image pixel that are not gas


final = segmentPorousMediaBurner( imgAir, imgKr, 5 );

%final( final1 ) = false;

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

figure; plot( att_curve, 'r' ); hold on;
plot( att_curve_bhc, 'b' );
xlabel 'slice #', ylabel 'attenuation';
axis([1 length( att_curve ) 0 0.015])
legend('Uncorrected', 'corrected')

figure;
slice = imgSub(:,end/2,:);
slice( final(:,end/2,:) ) = 0;
imagesc( squeeze( slice )' , [0 0.015] ); axis image, colorbar, colormap jet;
title 'Uncorrected'

figure;
slice = imgSubBHC(:,end/2,:);
slice( final(:,end/2,:) ) = 0;
imagesc( squeeze( slice )' , [0 0.015] ); axis image, colorbar, colormap jet;
title 'Beam hardening corrected'

