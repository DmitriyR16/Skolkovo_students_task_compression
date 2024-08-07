function [out, cmprxData_bytes] = MainDCTvalues(in,flg_direction)  
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
prm.nTX     = 64;
prm.nSb     = 50; 
prm.accuracy = 0.96;  % desired compression accuracy in fractions of unit
prm.bit_size_coeff_inds = 12; % a number of bits used to describe a single DCT coefficient's position
prm.bit_size_num_coeffs = 10; % a number of bits used to describe an amount of valuable coefficients
prm.bit_size_num_bits = 4; % a length of bit word containing the number of bits used for precise 
% bit representation of a single component of a complex coefficient
% ----------------------------------------------------------------------------------- 
  if   flg_direction
    % --------------- Compressing procedures ------------------------    
    data_array = in(:);  % an array of input matrix elements
    transf_data = DCT(data_array, flg_direction);  % a sequence of discrete cosine transform terms
    % transf_data = dct(data_array); % computationally efficient transform from Image Processing Toolbox
    [~,ind] = sort(abs(transf_data),'descend'); % obtain an array of indices sorted by absolute value
    transf_data_norm = sqrt(sum(transf_data.*conj(transf_data))); % the energy of the sequence X
    % Initialize a loop to define coefficients which represent (prm.accuracy*100)% of the
    % energy in the sequence X
    i = 1; % an initial index
    while sqrt(sum(transf_data(ind(1:i)).*conj(transf_data(ind(1:i))))) / transf_data_norm < prm.accuracy
       i = i + 1; % increase the index until the condition is met
    end
    out.coeff_inds = ind(1:i); % an array of indices of the most valuable DCT coefficients
    out.coeffs = transf_data(out.coeff_inds); % an array of the values, which have to be compressed
    out.num_coeffs = i; % a number of selected coefficients, considered to be up to 1024
    mmaxx = [max(abs(real(out.coeffs))) max(abs(imag(out.coeffs)))];
    out.num_bits = nextpow2(mmaxx) + 1; % a number of bits that is used to get an accurate bit
    % representation of a single component of a complex DCT coefficient,
    % considered to be up to 16 
    cmprxData_bytes = Byte_calcIn(out,prm);
  else
    % --------------- DeCompressing  procedures (Example)-------------------       
    data = load_file(in,prm);
    rest_transf_data = zeros(prm.nTX * prm.nSb, 1); % initialize an array for data restoration step
    rest_transf_data(data.coeff_inds) = data.coeffs; % fill it in with loaded values in appropriate positions
    rest_data = DCT(rest_transf_data, flg_direction); % decompressing 
    % rest_data = idct(rest_transf_data); % computationally efficient transform from Image Processing Toolbox
    out = reshape(rest_data, [prm.nTX, prm.nSb]); 
    cmprxData_bytes = 0;
  end
  
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
    fwrite(fid, in.num_coeffs, ['ubit' num2str(prm.bit_size_num_coeffs)]);
    fwrite(fid, in.num_bits, ['ubit' num2str(prm.bit_size_num_bits)]);
    fwrite(fid, real(in.coeffs)+2^(in.num_bits(1)-1), ['ubit' num2str(in.num_bits(1))]);
    fwrite(fid, imag(in.coeffs)+2^(in.num_bits(2)-1), ['ubit' num2str(in.num_bits(2))]);   
    fwrite(fid, in.coeff_inds, ['ubit' num2str(prm.bit_size_coeff_inds)]);  
    fclose(fid);
    fl=dir('tmp.bin');
    bt=fl.bytes;

end

function out = load_file(in,prm)
% ---------- Loading data from file ---------------------
% Input: 
% in    – input filepath with compressed data for reading.
% prm   – input parameters for loading 
% --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- 
% Output:
% out   - a structure array containing compressed data
% -----------------------------------------------------------------------------------
if ~isnumeric(in)   
    fid = fopen('tmp.bin','r'); 
    num_coeffs = fread(fid, 1, ['ubit' num2str(prm.bit_size_num_coeffs)]);
    num_bits = fread(fid, 2, ['ubit' num2str(prm.bit_size_num_bits)]);
    Real_vals = fread(fid, num_coeffs, ['ubit' num2str(num_bits(1))]);
    Imag_vals = fread(fid, num_coeffs, ['ubit' num2str(num_bits(2))]);
    out.coeff_inds = fread(fid, num_coeffs, ['ubit' num2str(prm.bit_size_coeff_inds)]);
    out.coeffs = Real_vals + 1i*Imag_vals - 2^(num_bits(1)-1)*ones(num_coeffs, 1) - 1i*2^(num_bits(2)-1)*ones(num_coeffs, 1);
    fclose(fid);
else
    out = in;
end 
       
end

function b = DCT(a, flag)
% ---------- Naive discrete cosine transform implementation ---------------------
% Input: 
% a    – an input array
% flag – binary variable determining a direction of the transform
% (flag = 1 for the direct transform and flag = 0 for the inverse one)
% --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- 
% Output:
% b - an output array of transform coefficients
% -----------------------------------------------------------------------------------
    N = length(a);
    b = zeros(1, N); % Initialize the output array
    for k = 0:N-1
        sum = 0;        
        for n = 0:N-1            
            if flag
                if k == 0
                    % Apply the scaling factor
                    alpha = sqrt(1/N);
                else
                    alpha = sqrt(2/N);
                end
                sum = sum + alpha * a(n+1) * cos(pi/N * (n + 0.5) * k);
            else
                if n == 0
                    % Apply the scaling factor
                    alpha = sqrt(1/N);
                else
                    alpha = sqrt(2/N);
                end                
                sum = sum + alpha * a(n+1) * cos(pi/N * (k + 0.5) * n);
            end
        end
        b(k+1) = sum;
    end
end