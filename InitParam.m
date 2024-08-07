function Param = InitParam(Scenario)
addpath 'func'
%  Params    - structure of Parameters  
%% Initialization parameters 
 Speed_idx =              2;% 1  - speed 0kmh; 2 - speed 3kmh; 3 - speed 5kmh
 SubBandNum = 24; 
 sig_term   = -167;  % termal noise level dB/Hz for UE
 
switch Scenario    
    case '4x4x2'    
         load ('E:\Channel_mobility\D-MIMO_Channel\TDD\Scenario-4x4x2\Channel_DMIMO.mat');    
         tmp = permute(H(Speed_idx,1).H_channel,[1 4 2 3]);  % RB\ant\Ues\TTIs
         Param.H_channel = tmp(:,:,:,6:end);  % shift on the 5 TTI in time 
    case '4x4x4'
         load ('E:\Channel_mobility\D-MIMO_Channel\TDD\Scenario-4x4x4\Channel_DMIMO_4TX.mat');             
         Param.H_channel = (H(Speed_idx,1).H_channel);     % RB\ant\Ues\TTIs
         disp ('necessary correct dims!!!') 
         return;
end
   file_pathloss='E:\Channel_mobility\D-MIMO_Channel\TDD\PathLoss.mat';  
     
  Param.PL =PathLos(file_pathloss); % calculation Path loss for each Ue in dB. 
  
%  Param.PathLoss = PathLoss; % in dB!
 Param.SRS_period = 12;  % SRS period of measurements (12 in TDD d-MIMO according 20ms)
 %% ------------- Powers --------------- 
 Param.Power_DMIMO          = 52;   % in dBm Power on all DMIMO cluster!    
 Param.Power_Termal_noise   = -174; % in dBm/Hz
 Param.Noise_ue             = -100; % dBm/5MHz (LTE)                                   [TS 36.101 Table 7.3.1-1] Band e-UTRA 38 (TDD 2.6GHz)
%  Param.snr_DL               = 20;   % SNR of DL channel in dB
 Param.P_srs_ue             = 20;  % Ue srs power in dBm/BW [TS 36.104 ]
       Noise_BS             = -101.5; % Reference noise level BS in dBm/BW[ 25Rb LTE] (TS 36.104 Table 7.2.1-1)!         
 %% --------------- Initialization parameters --- Antenna selection ------------------------------------------------
 Param.TxNum = size(H(Speed_idx,1).H_channel,4) / 4;% - number of antennas per DMIMO subcluster
 Param.TxGroup=           4;% - number of DMIMO subclusters       
 Param.RRUNum =           4;% - number of RRU in DMIMO subclusters
 Param.Scenario =        Scenario;
 Param.H_speed = H(Speed_idx,1).H_info;
 [Param.RbNum, Param.FullTxNum, Param.UeNum, Param.TTInum ] = size(Param.H_channel); 
%  Param.SNR_SRS=25;
 Param.Nbits = 8; % Number bits per component (I/Q)  
  
   
%% ----------------- Channel renormalization -------------------------------------------------
% Multiplication on sum pathloss for each ue 
    PL=sum(10.^(Param.PL/10));  % summary pathloss per UE
    PL_mas = permute(repmat(PL',[1 Param.RbNum Param.FullTxNum Param.TTInum]),[2 3 1 4]);
    Param.H_channel =  Param.H_channel.*sqrt(PL_mas);  % Initial channel with PL   h*PL(UE) 
%% ---------------------------------------------------------------------------------------------    
   Param.RbNum=24;   % definition number subBands per BW   
   Param.H_channel = squeeze(mean(reshape(Param.H_channel(1:end-1,:,:,:), ...
                              [(25-1)/Param.RbNum  Param.RbNum Param.FullTxNum, Param.UeNum, Param.TTInum ]),1));        
% search SNR_srs    

   N_srs = 24/Param.RbNum*6; % number SRS symbols in 1 RB;   
%    Noise_BS =  Param.Power_Termal_noise+10*log10(15000*12); % Termal noise level in dBm/RB
   
   SNR_SRS=10*log10(squeeze(sum(Param.H_channel.*conj(Param.H_channel),1)))-Noise_BS+10*log10(sqrt(N_srs)); % SNR_srs=||h||*P_srs^UE*sqrt(N_srs)/Noise_BS
    
%%
 Param.P_srs_ue=Param.P_srs_ue+10*log10(N_srs);  % Ue srs power in dBm/BW
 
 %% ----------------------- generation noise SRS according  ----------------------
 noise = 10^(Noise_BS/20)* 1/sqrt(2)*(randn(size(Param.H_channel))+1i*randn(size(Param.H_channel)));
 
 
 %% --------------------------------------------------------------------------

 for idxTTI=1:Param.TTInum 
     for idxUe = 1:Param.UeNum                 
     %% --------------------add noise in frequency ------------------    
       Param.H_channel_eq(:,:,idxUe,idxTTI) = Param.H_channel(:,:,idxUe,idxTTI)*10^(Param.P_srs_ue/20)+1/Param.RbNum*fft(noise(:,:,idxUe,idxTTI),[],1);  
       
    %% ----------------------- Quantization  ----------------------------------------------------------------------------
             maxx_I= max(max(real(Param.H_channel_eq(:,:,idxUe,idxTTI))));
             maxx_Q = max(max(imag(Param.H_channel_eq(:,:,idxUe,idxTTI))));
             mmaxx=max(maxx_I, maxx_Q);
            Param.H_channel_eq(:,:,idxUe,idxTTI) = fix(2^(Param.Nbits-1)/mmaxx*Param.H_channel_eq(:,:,idxUe,idxTTI));   %       quantization in fotmat I - Nbits / Q - Nbits 
    %% --------------------------------------------------------           
% % % % %          [~, amp] = normVectorsInMatrix(Param.H_channel(:,:,idxUe,idxTTI),1:25,1:32);
% % % % %          amp = sqrt(mean(Param.H_channel(:,:,idxUe,idxTTI).*conj(Param.H_channel(:,:,idxUe,idxTTI))));
% % % % %          Param.H_channel_eq(:,:,idxUe,idxTTI) = fix(2^7/mmaxx*(repmat(amp/10^(Param.SNS_SRS/20),[25 1]).*noise(:,:,idxUe,idxTTI)+(real(Param.H_channel(:,:,idxUe,idxTTI))+ ... 
% % % % %                          1i*imag(Param.H_channel(:,:,idxUe,idxTTI)))));
    
%          Param.H_channel_eq(:,:,idxUe,idxTTI) = fix(2^7/mmaxx*(repmat(amp*10^(-Param.SNS_SRS/20),[25 1]).*noise(:,:,idxUe,idxTTI)+Param.H_channel(:,:,idxUe,idxTTI)));
                     
%             nr=sum(mean(Param.H_channel_eq(:,:,idxUe,idxTTI).*conj(Param.H_channel_eq(:,:,idxUe,idxTTI)),1),2);         

            
 % --------------------- Normalization of quantized channel and calculation eigenvector ----------------------------           
            Param.H_channel_eq(:,:,idxUe,idxTTI)=normVectorsInMatrix(Param.H_channel_eq(:,:,idxUe,idxTTI).',1:32,1:Param.RbNum).';
%             Param.H_channel(:,:,idxUe,idxTTI)=normVectorsInMatrix(Param.H_channel(:,:,idxUe,idxTTI).',1:32,1:25).';
     end
 end  


%% ----Antenna selection ---------
 Param.TTIidx_Ant_selection = 1;    %  index of TTI for antenna selection  (1 time per second!!)
 
 for Ueidx = 1: Param.UeNum 
     
    H_tmp=squeeze(Param.H_channel_eq(:,:,Ueidx,Param.TTIidx_Ant_selection));
%% --------- 1 Stage RRU selection -----------------------------------------------
    Param.Alg_Selection.type= 0;% -  0/1 ('0' - 1 Stage RRU selection / '1' - 2 Stage RRU selection)
    Param.Alg_Selection.Thrd= 10;%-  threshold of antenna selection in dB (1 stage = 10dB, 2 stage = 13dB);    

    Param.Ant_selection_active(Ueidx) = RRU_Selection(H_tmp,Param);    
    
    Param.Alg_Selection.type= 0;% -  0/1 ('0' - 1 Stage RRU selection / '1' - 2 Stage RRU selection)
    Param.Alg_Selection.Thrd= 13;%-  threshold of antenna selection in dB (1 stage = 10dB, 2 stage = 13dB);    

    Param.Ant_selection_measure(Ueidx) = RRU_Selection(H_tmp,Param);        
    
    
 end       
  
 
%  UeSet =  [23 17 24 29 7 16 28 5];  % [2 4 6 11 14 22 24 29]; %1 : Param.UeNum;
% 
%  figure
%  id=1;
%  for idxUe=UeSet
%     
%     subplot(4,2,id)
%     idxAnt_mes= Param.Ant_selection_measure(idxUe).AntIdx_selected;
%     idxAnt_act= Param.Ant_selection_active(idxUe).AntIdx_selected;
%     P_srs_ue=[0:23];
%     plot((repmat(SNR_SRS(idxAnt_mes,idxUe,1),[1 24])+repmat(P_srs_ue,[length(idxAnt_mes) 1]))','-','LineWidth',1)
%     hold on 
%     plot((repmat(SNR_SRS(idxAnt_act,idxUe,1),[1 24])+repmat(P_srs_ue,[length(idxAnt_act) 1]))','*','LineWidth',2)
%     hold on 
%    
% %     plot ([zeros(1,thr-1) ones(1,24-thr+1)]*30,'k--')
%     title(['idxUE - ' num2str(idxUe)])
%     box off
%     grid on 
%     ylim([0 30])
%      id=id+1;
%  end
 
 

end