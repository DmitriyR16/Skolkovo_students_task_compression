function bt = Byte_calc(in,prm)
fid=fopen('tmp.bin','w');

%   fwrite(fid,in.numAntTx,'ubit7');  
%   fwrite(fid,in.numSubbands,'ubit7')  ;
 fwrite(fid,in.PointValues(2:end,1,:),['ubit' num2str(prm.amplitudeFormula(1))]) ; 
 fwrite(fid,in.PointValues(:,2,:),['ubit' num2str(prm.amplitudeFormula(2))])  ;
 fwrite(fid,in.PointIndexes,['ubit' num2str(nextpow2(in.numAntTx*in.numSubbands))]);  
%  fwrite(fid,in.Peak_Ref_Amp,'float'); 

 fclose(fid);
 fl=dir('tmp.bin');
 bt=fl.bytes;
 delete('tmp.bin');
end