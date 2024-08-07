clear all
% ----- Random MISO channel generation ------------------
Ntx = 64;
Sc = 1200;
D=1e5;
alpha = 1;
Rank = 2;
Hphase = rand(Ntx,Sc)*pi/6; % iniform phase in [0, 2pi]
for idNtx = 1:Ntx
  Hamp(idNtx,:) = 1/sqrt(2*alpha*pi*D)*exp(-([-Sc/2:Sc/2-1]-randn*Sc/2).^2/2/D);
end



H=Hamp.*exp(1i*Hphase);
pw = sum(H.*conj(H),2);
H=H./sqrt(repmat(pw,[1 Sc]))/Ntx;
Hn=zeros(size(H));
Hn([1 4 7 32+[1 4 7]],:)=H([1 4 7 32+[1 4 7]],:);

figure, 
 subplot 121  
 mesh(abs(H))  
 view([40 10 120])
 
% Hn=reshape(fft(fft(reshape(Hn,[ 8 8 Sc]),8,1),8,2),[64 1200]);
Hnn=fft(Hn,64,1);
[u,s,v] = svd(H);
s(Rank+1:end,Rank+1:end)=0;
Hreg = u*s*v';

 
 subplot 122  
 mesh(abs(Hreg))  
 view([40 10 120])
 
 wishrnd