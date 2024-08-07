function out = Scena1(Hin) 
% Defenition senario for channel redundancy in Frequency/Space domain 
% -----------------------------------------------------------------------------------
% Inputs:
%  Hin - channel data of following demensions:
%        [Ue Antennas/ BS Antennas / Frequency subcarriers / Realizations(T) ]
% -----------------------------------------------------------------------------------
% Outputs:
%  out - structure of variables 
%     out.ScNum          - frequency subcarriers number    
%     out.ScSB           - number subcarriers in 1 subband 
%     out.SbNum          - frequency subband numbers
%     out.IQbits         - bits number per Real/Imagenary part of complex number of channel H 
%     out.H_channel_eq   - quantized channel by predefined bits number forman (IQbits)
%              [Ue Antennas/ BS Antennas / Frequency subbands / Realizations(T)]
%     out.H_channel_eigV - quantized 1st eigenvectors of channel 'out.H_channel_eq' 
%                          by predefined bits number format (IQbits)
%              [BS Antennas / Frequency subbands / Realizations(T)]
%-------------------------------------------------------------------------------------
%% Initialization parameters  
[out.RxNum, out.TxNum, nSC , out.Snapshots] = size(Hin); 
 out.ScNum      = nSC; 
 out.ScSB       = 1; 
 out.SbNum      = out.ScNum/out.ScSB;  
 out.IQbits     = 8;  % bits number per Real/Imagenary part of complex number of channel H  
 flg_H_Aver     = 'mean'; % 'mean' or 'svd' 'decim' - different type of channel averaging
 
%% ----------------- Channel averaging ------------------------------------------------------
if out.ScNum>nSC
   disp('Error! Cutted frequency subband is out of subcarriers set'); 
   return
end
    
Hin0 = Hin(:,:,nSC/2-out.ScNum/2+1:nSC/2+out.ScNum/2,:);

Hin_ext = permute(reshape(Hin0,[out.RxNum, out.TxNum, out.ScSB, out.SbNum, out.Snapshots]),[2 3 1 4 5]);  
switch flg_H_Aver
    case 'svd'
        for idSnap = 1:out.Snapshots            
           for idRx = 1:out.RxNum 
            for idRb = 1:out.SbNum
              tmp = Hin_ext(:,:,idRx,idRb,idSnap);     
              [u,~,~] = svd(tmp'*tmp);
              Hin_aver(idRx,:,idRb,idSnap)=tmp*u(:,1);   
            end
           end
        end
    case 'mean'
       Hin_aver = permute(mean(Hin_ext,2),[3 1 4 5 2]);      
    case 'decim'   
       Hin_aver = permute(Hin_ext(:,1,:,:,:),[3 1 4 5 2]);      
end 


%% ----------------- Channel renormalization -------------------------------------------------
     tmp = sum(sum(sum(Hin_aver.*conj(Hin_aver),3),2),1); 
     Hin_norm = Hin_aver./sqrt(repmat(tmp,[out.RxNum out.TxNum out.SbNum 1]));    
%      Hin_norm = Hin_aver;
 
% ----------------------- generation noise   (option)  ----------------------
%  noise =  1/sqrt(2)*(randn(size(Hin_norm))+1i*randn(size(Hin_norm))); 
% --------------------------------------------------------------------------

 for idxTTI=1:out.Snapshots     
% --------------------add noise in frequency (optional) ------------------    
%    out.H_channel_eq(:,:,idxUe,idxTTI) = Hin_norm(:,:,idxUe,idxTTI)*10^(Param.P_srs_ue/20)+1/Param.RbNum*fft(noise(:,:,idxUe,idxTTI),[],1);         

    %% ----------------------- Quantization  ---------------------------------------------------------
         tmp = Hin_norm(:,:,:,idxTTI);      
         mmaxx=max([real(tmp(:))' imag(tmp(:))']);      
         out.H_channel_eq(:,:,:,idxTTI) = fix(2^(out.IQbits-1)/mmaxx*Hin_norm(:,:,:,idxTTI));   %       quantization in fotmat I - out.IQbits / Q - out.IQbits                 
         for idRb = 1:out.SbNum
           tmp = out.H_channel_eq(:,:,idRb,idxTTI);   
           [u,~,~] = svd(tmp'*tmp);
           mmaxx = max([real(u(:,1))' imag(u(:,1))']);
           out.H_channel_eigV(:,idRb,idxTTI) = fix(2^(out.IQbits-1)/mmaxx*u(:,1));                 
         end 
 end  

end