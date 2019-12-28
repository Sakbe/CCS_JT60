close all
clear all

%%%%%%%%%%%%%%%%%%%%% CCS - FBC - XCS - CREATE configuratio file
%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%


%%%%%%% Select # of control points (19,8 or 6)
    nOfFluxCntrl =19;
    
    %%%%% Select mode: XSC controller (mode =1) , shape will change during transition (mode =2), FBC
    %%%%% controller (mode=3)
    mode=3;
    
      %%% select disturbance 
    %%! 1->URano  2->comp ELM  3->ELM   4->minor disrup  5-> ITER-like
    %%minor disrup
    
    dist=5;
    

if (nOfFluxCntrl== 8 ||nOfFluxCntrl== 6 || nOfFluxCntrl == 19)

    
    if(mode ==1)
      
    
    if(dist>=1 && dist<= 5)
             JT_60SA_SOF18d66_japs_sensors_isoflux_simulator
             Smlnkmodel=('JT_60SA_scheme_japanese_sensors_isoflux_test_XSC');
            open_system(Smlnkmodel);

    else
            disp('Not valid disturbance')
    end
    
    elseif(mode == 3)
    if(dist>=1 && dist<= 5)
            JT60_CD_Invessel_simulator
            Smlnkmodel=('JT_60SA_scheme_isoflux_CCS_FBC_CDinvessel');
            open_system(Smlnkmodel)
    else
            disp('Not valid disturbance')
    end 
    
    elseif(mode == 2)
        
            if(dist>=1 && dist<= 5)
            JT_60SA_simulator_VS_SOF18d66_japs_sensors_isoflux_sqzd
            Smlnkmodel=('JT_60SA_scheme_japanese_sensors_isoflux_XSC_sqzd');
            open_system(Smlnkmodel);
    else
            disp('Not valid disturbance')
    end 

    else
        disp('Not valid selection of mode ')
    end
else
    disp('Not valid number of control points')
end