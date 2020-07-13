function [gci]=epochExtract(zfSig)


s=zfSig>=0;
k=s(2:end)-s(1:end-1);
gci=find(k>0);

p=zeros(length(zfSig),1);  
	
% gci=gci+3; 
p(gci)=1;
gci=p;