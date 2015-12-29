% load reconstruction parameters
load 'temp.mat';

geom = loadProjectionGeometryCT( p );
geom.detOffset(2) = geom.detOffset(2) - 64;

spectrum = loadSpectraCT(p, geom, 2e6);

Dir = 'D:\TabletopScannerData\Nov_10_2015_Study\3ppi_interface_60kV_50mA\';

%% load air scan data

dataPathAir = [Dir 'air_cold_04' '\'];
process_seq_file( Dir, 'air_cold_04' );
[sinoAttAir, sinoPCAir] = loadTableTopData( dataPathAir, geom, 0, [701 800], [101 500], 0 , -64 );

dataPathAir = [Dir 'air_cold_05' '\'];
process_seq_file( Dir, 'air_cold_05' );
sinoAttAir = sinoAttAir + loadTableTopData( dataPathAir, geom, 0, [701 800], [101 500], 0 , -64 );
sinoAttAir = sinoAttAir / 2;

% load burn scan data

dataPath = [Dir 'background_01' '\'];
process_seq_file( Dir, 'background_01' );
[sinoAtt, sinoPC] = loadTableTopData( dataPath, geom, 0, [701 800], [101 500], 0, -64 );

dataPath = [Dir 'background_02' '\'];
process_seq_file( Dir, 'background_02' );
sinoAtt = sinoAtt + loadTableTopData( dataPath, geom, 0, [701 800], [101 500], 0, -64 );
sinoAtt = sinoAtt / 2;

%% first pass reconstruction

sinoAttAirPoly = beamHardeningMaterialCorrection(sinoAttAir, spectrum, 'Quartz', 3 );
imgAir = reconFBP( sinoAttAirPoly, geom, 'hamming' );

sinoAttPoly = beamHardeningMaterialCorrection(sinoAtt, spectrum, 'Quartz', 3 );
imgKr = reconFBP( sinoAttPoly, geom, 'hamming' );

imgSub = imgKr - imgAir;

clear sinoAttAirPoly sinoAttPoly;

%% second pass beam hardening correction

mapTube = single( imgAir > 0.6 );
sinoTube = forwardProjectMex( mapTube, geom ) ;
sinoTube = imfilter3( sinoTube, fspecial('gaussian', [5 5], 1 ) );

sinoSubBHC = beamHardeningMaterialCorrectionBurner(sinoAtt - sinoAttAir, sinoAttAir, sinoTube, spectrum);

imgSubBHC = reconFBP( sinoSubBHC, geom, 'hamming' );

%% PWLS reconstruction

weights = computeWeightsPwls( ( sinoPC + sinoPCAir ) / 2 , 0, spectrum.electronicNoise );

inner = segmentPorousMediaBurner( imgAir, imgKr, 1 );
burner = segmentPorousMediaBurner( imgAir, imgKr, 5 );
burner( inner ) = false; 

img_pwls_1 = reconPMBPwlsLALMOs14( sinoSubBHC, weights, geom, 10, 'isotv', 30, 0, 12, imgSubBHC, burner );

img_pwls_2 = reconPMBPwlsLALMOs14( sinoSubBHC, weights, geom, 1, 'isotv', 30, 0, 12, imgSubBHC, burner );

img_pwls_3 = reconPMBPwlsLALMOs14( sinoSubBHC, weights, geom, 0.1, 'isotv', 30, 0, 12, imgSubBHC, burner );

%% Display results

figure;
slice = imgSubBHC(:,end/2,:);
imagesc( squeeze( slice )' , [0 0.015] ); axis image, colorbar, colormap jet;
title 'FBP reconstruction'

figure;
slice = img_pwls_1(:,end/2,:);
imagesc( squeeze( slice )' , [0 0.015] ); axis image, colorbar, colormap jet;
title 'beta = 10'

figure;
slice = img_pwls_2(:,end/2,:);
imagesc( squeeze( slice )' , [0 0.015] ); axis image, colorbar, colormap jet;
title 'beta = 1'

figure;
slice = img_pwls_3(:,end/2,:);
imagesc( squeeze( slice )' , [0 0.015] ); axis image, colorbar, colormap jet;
title 'beta = 0.1'


