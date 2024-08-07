function strategies = InitCmprx(idx)
       

        %for idx =1: in,
            
            switch idx
             case 1  
             strategies.points_per_antenna = 16;
             
             case 2  
             strategies.points_per_antenna = 24;             
                 bits_Ampl_first_peak=8;
                 bits_phase_first_peak=8;
                 bits_Ampl_residual_peak=8;
                 bits_phase_residual_peak=8;

             case 3  
             strategies.points_per_antenna = 4;
                 bits_Ampl_first_peak=4;
                 bits_phase_first_peak=6;
                 bits_Ampl_residual_peak=4;
                 bits_phase_residual_peak=6;             
             case 4  
             strategies.points_per_antenna = 2;
             
            end
            
            vector_amplitude=[bits_Ampl_first_peak, bits_Ampl_residual_peak*ones(1,strategies.points_per_antenna-1)];
            vector_phase=[bits_phase_first_peak, bits_phase_residual_peak*ones(1,strategies.points_per_antenna-1)];
            strategies.amplitudeFormula = vector_amplitude;
            strategies.phaseFormula = vector_phase;
            strategies.isDifferentialAmplitude = 0;
            strategies.compression_scenario = 3; % Use only IFFT for frequency correlation (get taps in time domain) (D-MIMO)
            
        %end
            
       
end