% load reconstruction parameters
load 'temp.mat';

geom = loadProjectionGeometryCT( p );
geom.detOffset(2) = -  geom.detOffset(2);

spectrum = loadSpectraCT(p, geom, 2e6);

Dir = 'E:\Data\NasaFlame\Nov_10_2015_Study\3ppi_interface_60kV_50mA\';

%% load air scan data

dataPathAir = [Dir 'air_cold_05' '\'];
process_seq_file( Dir, 'air_cold_05' );
sinoAttAir = loadTableTopData( dataPathAir, geom, 0, [701 800], [101 500] );

% load burn scan data

dataPath = [Dir 'background_01' '\'];
process_seq_file( Dir, 'background_01' );
sinoAtt = loadTableTopData( dataPath, geom, 0, [701 800], [101 500] );


%% first pass reconstruction

sinoAttAirPoly = beamHardeningMaterialCorrection(sinoAttAir, spectrum, 'Quartz', 3 );
imgAir = reconFBP( sinoAttAirPoly, geom, 'hamming' );

sinoAttPoly = beamHardeningMaterialCorrection(sinoAtt, spectrum, 'Quartz', 3 );
imgKr = reconFBP( sinoAttPoly, geom, 'hamming' );

imgSub = imgKr - imgAir;

clear sinoAttAirPoly sinoAttPoly;

%% second pass beam hardening correction

mapTube = single( imgAir > 0.7 );
sinoTube = forwardProjectMex( mapTube, geom ) ;
sinoTube = imfilter3( sinoTube, fspecial('gaussian', [5 5], 1 ) );

%%

sinoSubBHC = beamHardeningMaterialCorrectionBurner(sinoAtt - sinoAttAir, (sinoAtt + sinoAttAir)/2, sinoTube, spectrum);

imgSubBHC = reconFBP( sinoSubBHC, geom, 'hamming' );

