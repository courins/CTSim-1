% load reconstruction parameters
load 'temp.mat';

geom = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 2e6);

dir = 'D:\MATLAB\CTData\Dec_01_2015_Study\3ppi_sic_60KV_50mA_lowflow\';

file_list = { 'air_00'; 'air_01';  'air_02'; 'air_03'; 'air_04'; 'air_05';  'air_06'; 'air_07'; ...
    'burn_00'; 'burn_01'; 'burn_02'; 'burn_03'; 'burn_04'; 'burn_05'; 'burn_06';  };

for i = 1 : length( file_list )
    
    seq_filename = file_list{i};
    
    dataPath = [dir seq_filename '\'];
    
    process_seq_file( dir, seq_filename )
    
    sinoAtt = loadTableTopData( dataPath, geom );
    
    %% first path reconstruction
    
    sinoAttPoly = beamHardeningMaterialCorrection(sinoAtt, spectrum, 'Quartz', 3 );

    img = reconFBP( sinoAttPoly, geom, 'hamming' );
        
    %% second pass beam hardening correction
    
    mapTube = single( img > 0.65 );
    
    sinoTube = forwardProjectMex( mapTube, geom ) ;
    
    sinoAttBHC = beamHardeningMaterialCorrectionBurner(sinoAtt, sinoTube, spectrum);
    
    % final reconstruction
    
    imgBHC = reconFBP( sinoAttBHC, geom, 'hamming' );
    
    save([seq_filename '.mat'], 'img', 'imgBHC' );
    
end



