% load reconstruction parameters
load 'temp.mat';

geom = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 2e6);

dir = 'D:\MATLAB\CTData\Dec_01_2015_Study\3ppi_sic_60KV_50mA_lowflow\';

file_list_air = { 'air_01';  'air_02'; 'air_03'; 'air_04'; 'air_05';  'air_06'; 'air_07';  };

file_list_burn = { 'burn_00'; 'burn_01'; 'burn_02'; 'burn_03'; 'burn_04'; 'burn_05'; 'burn_06';  };

for i = 1 : min(  length( file_list_air ), length( file_list_burn ) )
    
    
    %% load air scan data
    dataPathAir = [dir file_list_air{i} '\'];
    process_seq_file( dir, file_list_air{i} );
    sinoAttAir = loadTableTopData( dataPathAir, geom, 0, [701 800], [101 500] );
    
    % load burn scan data
    dataPath = [dir file_list_burn{i} '\'];
    process_seq_file( dir, file_list_burn{i} );
    sinoAtt = loadTableTopData( dataPath, geom, 0, [701 800], [101 500] );
    
    
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
    
    %%
    
    sinoSubBHC = beamHardeningMaterialCorrectionBurner(sinoAtt - sinoAttAir, sinoAttAir, sinoTube, spectrum);
    
    imgSubBHC = reconFBP( sinoSubBHC, geom, 'hamming' );
    
    save([seq_filename '.mat'], 'imgAir', 'imgKr', 'imgSub', 'imgSubBHC');
    
end



