% Taskl #1 Channel compression 
%% ------------- Load Channel-------------------------------
% Structure of channel H: 
%  'H' - channel data of following demensions:
%        [Ue Antennas/ BS Antennas / Frequency subcarriers / Realizations(T) ] - complex double format       
       load('Channel_test2.mat')  
%% --- Initialization of parameters and Channel redundancy -------------

   Param = Scenarios.Scena1(Hh2(:,:,:,1:100));       
  %     Param.H_channel_eigV  - array include BaseLine of Eigenvectors of
  %     input channel 'H' for compression !
%% ---------------------------------------------------------------------  
 name_of_cmprx = dir('+Cmprx\*m'); 
 [NumTX, NumRb, NumTTIs]=size(Param.H_channel_eigV);
 Results =[];
 lg1=[];lg2=[];
 
  figure (101)
 LineCl='brgkcmy'; 
 LineTp=' o*xsd^v><ph';
%% -------------------- Compare different compression / decompression algorithms ------------ 
 for idxAlg = 1:2 %size(name_of_cmprx)     
     BaseLine_bytes = NumTX*NumRb*2*Param.IQbits/8; %NumTX*NumRb*2*IQbits/8 [bytes] 
     for idxTTI =1:NumTTIs                                
             alg_name=name_of_cmprx(idxAlg).name; 
             Htmp  = Param.H_channel_eq(:,:,:,idxTTI); % Channel 
             HeigV = Param.H_channel_eigV(:,:,idxTTI); % 1st eiginvector of channel
       %% ----------- Compression/Decompression Algorithms -----------------              
          % ----------- BaseLine Moscow RTT compression ------------
          [CmprxH,  cmprxData_bytes]   = Cmprx.(alg_name(1:end-2))(HeigV,1);    
          [RestoreH, ~] = Cmprx.(alg_name(1:end-2))('tmp.bin',0);                                     
       
       %% --------  Error calculation --------------------------------------          
        nr1 = sqrt(sum(RestoreH.*conj(RestoreH),1))';
        nr2 = sqrt(sum(Param.H_channel_eigV(:,:,idxTTI).*conj(Param.H_channel_eigV(:,:,idxTTI)),1))';
          Results(idxAlg).Err(idxTTI)      = ...              
              mean(1-abs(diag(RestoreH'*Param.H_channel_eigV(:,:,idxTTI)))./nr1./nr2)*100;   %Error in [%]                
%           Results(idxAlg).CmprxH(idxTTI)   = CmprxH;              
%           Results(idxAlg).RestoreH(idxTTI) = RestoreH;              
          Results(idxAlg).name           = alg_name(1:end-2);              
          Results(idxAlg).Rate(idxTTI)   = BaseLine_bytes/cmprxData_bytes;              
     end
     %% ----------- Error plotting -----------------------------
      er=Results(idxAlg).Err;
      er_bt=mean(Results(idxAlg).Rate);
      alg_name(alg_name=='_')=' ';
      lg1{idxAlg}=['[\bf \mu(Err)=' num2str(fix(mean(er*100))/100) '%]-' Results(idxAlg).name]; 
      lg2{idxAlg}=['[\bf \mu(Rate)=' num2str(fix(mean(er_bt*100))/100) ']-' Results(idxAlg).name]; 
      
    subplot (1,2,1),       
       hold on, plot(mean(er),mean(er_bt),[LineCl(mod(idxAlg,6)) LineTp(fix(idxAlg/5)+2)],'LineWidth',2)
%      legend(lg2,'Location','southoutside','FontSize',12);
       legend(lg2);
       title('\bf Results of compare Compression Algorithms','FontSize',16)
       ylabel('\bf Compression Rate of channel reconstruction','FontSize',14)
       xlabel('MSE of channel reconstruction, %','FontSize',14)
%      set(gca, 'yscale','log');
%     ylim([0 100]);
    grid on 
    box off
   
   subplot (1,2,2),   
     hold on, plot(er,[LineCl(mod(idxAlg,6)) '-' LineTp(fix(idxAlg/5)+1)],'LineWidth',2)
%    legend(lg1,'Location','southoutside','FontSize',12);
     legend(lg1);
     title('Results of compare Reconstruction Error','FontSize',16)
     ylabel('MSE of channel reconstruction, %','FontSize',14)
     xlabel('TTIs','FontSize',14)
     set(gca, 'yscale','log'); 
     ylim([0 100]);
     
    grid on 
    box off   
 end
 

