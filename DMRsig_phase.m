function [Y] = DMRsig(filename, RD,RP,phase, App, M, Fs, NS, NzDS);
% make dynamic moving ripple stimulus
% stores signal Y as floats in filename.bin and other variables in
%   
RDint   = interp1(NzDS*(0:Mnfft-1),RD,0:M-1,'cubic'); 
RPint   = interp1(NzDS*(0:Mnfft-1),RP,0:M-1,'cubic'); 

wb = waitbar(0,'Generating DMR signal');
for k=1:NS
    A = sin(2*pi*RDint*X(k)+ RPint);
    Y=Y + App*10.^(A/20).* sin( 2*pi*faxis(k)*(1:M)/Fs + phase(k) );
    waitbar(k/NS,wb);
end
close(wb)
% saving data and parameters
tofloat([filename '_dB.bin'],Y);
