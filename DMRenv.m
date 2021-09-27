function [RD RP FM faxis X FsN MN MNfft] = DMRenv(flo, fhi, fRD, fFM, MaxRD, MaxFM, M, Fs, NS, NzDS,seed);
% [RD RP FM faxis X] = DMRenv(flo, fhi, fRD, fFM, MaxRD, MaxFM, M, Fs, NS, NzDS,seed);
% caculate envelope for dynamic moving ripple 
%    RD, FM, RP - Ripple density, temporal modulation, ripple phase
%    faxis,X - frequencies list in Hz/octaves above flo
%    FsN, MN - noise sampling frequency, 
%    MN, MNfft - number of samples for ripple noise, samples needed for fft
%   

%
FsN=Fs/NzDS;			%Noise Sampling Frequency
MN=floor(M/NzDS)+1;			%Noise Signal Length
MNfft=2^(ceil(log2(MN)));	%Next power of 2 - Used to calculate bandlimited uniform noise

XMax=log2(fhi/flo); % range of frequencies in octaves
X=XMax*(0:NS-1)/(NS-1); % frequencies list in octaves above flo
faxis=flo*2.^X;

if exist('seed')
	seedt=seed;
else
	seedt=1;
end

%Generating Ripple Density Signal (Uniformly Distributed Between [0 , MaxRD] )
RD=MaxRD*noiseunif(fRD,FsN,MNfft,seedt);

%Generating the Teporal Modulation Signal 
%Signal Varies Between [-MaxFM , MaxFM]
FM=2*MaxFM*(noiseunif(fFM,FsN,MNfft,seedt+2)-.5);

%Genreating Ripple Phase Components
%Note that: RP = 2*pi/FsN*intfft(FM);
% RP=2*pi/FsN*intfft(FM);
RP = (2*pi/FsN)*cumsum(FM);

disp('Done generating RD,FM and RP.');

%For Displaying Statistics
%figure(1),hist(RD,100),title('Rippe Density Distribution')
%figure(2),hist(FM,100),title('Modulation Frequency Distribution')
%figure(3),plot((1:length(FM))/FsN,FM),title('Modulation Frequency')
%figure(4),plot((1:length(RD))/FsN,RD),title('Ripple Density')


