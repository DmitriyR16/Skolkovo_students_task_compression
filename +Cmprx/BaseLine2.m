function [out, cmprxData_bytes] = BaseLine(in,flg_direction)  
% ---------------------------- Description ------------------------------------------
% Inputs:
% in            – input data of size [N x Nf] (for compression: flg_direction=1)  
% in            – full path to file ‘tmp.bin’     (for decompression: flg_direction=0)
% flg_direction – flag of compression / decompression => respectevly (1/0) 
% --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- 
% Outputs:
% out             - compressed version of input data   ‘in’ (for compression: flg_direction=1)  
% out             - decompressed version of input data ‘in’ (for decompression: flg_direction=0)  
% cmprxData_bytes - size of compressed channel in bytes 
% -----------------------------------------------------------------------------------
prm.bits    = 3;
prm.nTX     = 64;
prm.nSb     = 128;
% ------------- Initialization compression/decompression algorithm ------------------        
cmp = Initialization(); % if necessary to initialize parameters or define class  
% ----------------------------------------------------------------------------------- 
  if   flg_direction
    % --------------- Compressing procedures (Example) ------------------------    
     mmaxx = max([real(in(:))' imag(in(:))']);
     out = fix(2^(prm.bits-1)/mmaxx*in); % compressing     
     cmprxData_bytes = Byte_calcIn(out,prm); 
  else
    % --------------- DeCompressing  procedures (Example)-------------------       
     data = load_file(in,prm);
     out =  data; % decompressing 
     cmprxData_bytes = 0;
  end
  
end

function out = load_file(in,prm)
% ---------- Loading data from file ---------------------
if ~isnumeric(in)  
      fid = fopen(in,'r')  ; 
      tmp = fread(fid,[prm.nTX prm.nSb*2],['ubit' num2str(prm.bits)]);    
      out = (tmp(:,1:end/2)+1i*tmp(:,end/2+1:end))-(1+1i)*2^(prm.bits-1);
      fclose(fid) ;  
else
    out = in;
end 
    
    
end

function out = Initialization()
% ---------------------------- Description ------------------------------------------
% Here can be defined structure, variables, which can help 
% to evaluation compressed and decompressed version of channel measurements  
out = 0;
% ................
end

function bt = Byte_calcIn(in,prm)
% ---------------------------- Description ------------------------------------------
% Input: 
% in    – input compressed data for saving into file.
% prm   – input parameters for save assistance 
% --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- 
% Output:
% bt    - size of file ‘tmp.bin’ in bytes 
% -----------------------------------------------------------------------------------

 fid=fopen('tmp.bin','w');
 fwrite(fid,[real(in(:))+2^(prm.bits-1) imag(in(:))+2^(prm.bits-1)],['ubit' num2str(prm.bits)]);      
 fclose(fid);
 fl=dir('tmp.bin');
 bt=fl.bytes;
%  delete('tmp.bin'); 

end