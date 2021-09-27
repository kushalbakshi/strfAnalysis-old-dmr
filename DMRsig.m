function [Y phase] = DMRsig(filename, RD,RP,X,faxis, App, M, Fs, NS, NzDS);
% [Y phase] = DMRsig(filename, RD,RP,X,faxis, App, M, Fs, NS, NzDS);
% make dynamic moving ripple stimulus
% stores signal Y as floats in filename.bin and other variables in
%   


Y = zeros(1,M); % signal

%Initial Carrier Sinusoid Phases 
phase=2*pi*rand(1,NS);

% interpolate RD and RP
disp('Interpolating RD and RP.')
RDint   = interp1(NzDS*(0:length(RD)-1),RD,0:M-1,'pchip'); 
RPint   = interp1(NzDS*(0:length(RP)-1),RP,0:M-1,'pchip'); 
disp('Done.')

% open file
fullfilename = [filename '_dB.bin'];
fid = fopen(fullfilename,'r');

% generate signal one second at a time and write to data file carrier at a time
numsec = floor(M/Fs);
extrasamp = M-numsec*Fs;
wb = waitbar(0,'Generating DMR signal');
for i=1:numsec
    disp([num2str(i) '/' num2str(numsec)]);
    inds = (i-1)*Fs+(1:Fs);
    for k=1:NS
        A = (sin(2*pi*X(k)*RDint(inds)+ RPint(inds))-1)/2;
        Y(inds)=Y(inds) + (10.^(App*A/20)).* sin( 2*pi*faxis(k)*(inds-1)/Fs + phase(k) );
    end
    % saving data and parameters
    tofloat([filename '_dB.bin'],Y(inds));
    waitbar(i/numsec,wb);
end
if floor(M/Fs)<M/Fs
    inds = (numsec*Fs+1):M;
    % interpolate to get envelope at sample frequency
    RDint   = interp1(NzDS*(0:length(RD)-1),RD,inds-1,'pchip'); 
    RPint   = interp1(NzDS*(0:length(RP)-1),RP,inds-1,'pchip'); 
    for k=1:NS
        A = (sin(2*pi*X(k)*RDint+ RPint)-1)/2;
        Y(inds)=Y(inds) + (10.^(App*A/20)).* sin( 2*pi*faxis(k)*(inds-1)/Fs + phase(k) );
    end
    % saving data and parameters
    fwrite(fid,Y(inds),'float');
end    
    
close(wb)
close(fid)

