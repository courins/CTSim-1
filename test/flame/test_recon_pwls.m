% load reconstruction parameters
load 'temp.mat';

geom = loadProjectionGeometryCT( p );

spectrum = loadSpectraCT(p, geom, 2e6);

dir = 'D:\MATLAB\CTData\Dec_01_2015_Study\3ppi_sic_60KV_50mA_lowflow\';

file_list = { 'air_05';  'burn_05';  };

for i = 1 : length( file_list )
    
    seq_filename = file_list{i};
    
    dataPath = [dir seq_filename '\'];
    
    process_seq_file( dir, seq_filename )
    
    sinoAtt = loadTableTopData( dataPath, geom );
    
    sinoAtt = beamHardeningMaterialCorrection(sinoAtt, spectrum, 'Quartz', 10 );
    
    img = reconFBP( sinoAtt, geom, 'hamming' );
    
    save([seq_filename '.mat'], 'img' );
    
end



