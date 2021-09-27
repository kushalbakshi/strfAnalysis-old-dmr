%% Set parameters
bat_list = {'Tb90'};
data_location='S:\Smotherman_Lab\Auditory cortex\';
sr = 40000;
[b,a] = butter(4,[800/(sr/2), 3000/(sr/2)],'bandpass');
cmap=fireice(256);
channels_available=32;
STRF_save_path='U:\STRF Analysis\MUA STRFs\';
MTF_save_path = 'U:\STRF Analysis\TMR and SMR\';
shank_pos=xlsread('Shank Position.xlsx');

for bat_num = 1:numel(bat_list)
    d = dir(strjoin([data_location,bat_list(:,1),'\Matfile'],''));
    site_number = numel(d)-2;
    for site = 4 %1:site_number
        %% Import Data, preprocess, and extract spike times
        load(strjoin([data_location,string(bat_list(:,bat_num)),'\Matfile\',...
            d(site+2,1).name,'\event.mat'],''));
        cutoff = ts(end)*sr;
        marker = importdata(strjoin([data_location,string(bat_list(:,bat_num)),...
            '\Data\',string(bat_list(:,bat_num)),'_',num2str(site),'_marker_tc.mat'],''));
        
        for channel = 28%1:channels_available
            if isfile(strjoin([data_location,string(bat_list(:,bat_num)),'\Matfile\',...
                    d(site+2,1).name,'\Chn',num2str(channel),'.mat'],''))
                
                readfile = matfile(strjoin([data_location,string(bat_list(:,bat_num)),...
                    '\Matfile\',d(site+2,1).name,'\Chn',num2str(channel),'.mat'],''));
                data = readfile.data;
                data(1:cutoff) = [];
                data = filtfilt(b,a,data);
                thr = std(data)*4;
                [peaks,locs] = findpeaks(data*-1,'MinPeakHeight',thr);
                locs = (locs/sr)*1000;
                depth = marker.depth-shank_pos(channel,1);
                %% Find STRF
                spk = locs;
                spiketimes{1}=spk;
                processdata
                STA=DMRcells.STA;
%                 STA_avg=mean(STA(:));
%                 STA_stdev=std(STA(:))*2.5;
%                 minus_STA_stdev=STA_stdev*-1;
%                 STA(STA > minus_STA_stdev & STA < STA_stdev) = 0;
                [col,row]=find(STA==max(STA(:)));
                latency=taxis(row);
                peak_freq=faxis(col)/1000;
                %% Create STRF Figure
                figure('WindowState','maximized')
                set(gcf,'color','w')
                ss = surface(taxis,X,STA);
                set(ss,'edgecolor','none');
                xlabel('Time (ms)', 'FontSize',28,'FontName','Trebuchet MS')
                ylabel('Frequency (kHz)', 'FontSize',28,'FontName','Trebuchet MS')
                xlim([0 200])
                set(gca,'FontSize',24,'FontName','Trebuchet MS')
                set(gca,'ydir','normal');
                set(gca,'xdir','reverse');
                set(gca,'view',[0 90])
                trange = [min(taxis) max(taxis)];
                %             set(gca,'xlim',trange);
                %set(gca,'ylim',[min(faxis) max(faxis)]/1000)
                colormap(cmap)
                lim=caxis;
                caxis([-lim(1,2),lim(1,2)])
                c_bar=colorbar;
                ylabel(c_bar,'Firing rate')
                annotation('textbox',[.65 .9 .1 .1],'string',...
                    ['Depth: ',num2str(depth),' µm'],'Edgecolor','none',...
                    'FontSize',14,'FontName','Trebuchet MS')
                
                annotation('textbox',[.65 .87 .1 .1],'string',...
                    ['Latency of peak excitation: ',num2str(latency),' ms'],...
                    'Edgecolor','none','FontSize',14)
                
                annotation('textbox',[.75 .902 .1 .1],'string',...
                    ['Frequency at peak excitation: ',num2str(peak_freq),' kHz'],...
                    'Edgecolor','none','FontSize',14)
                
                title(['Site ',num2str(site), ' Chn ',num2str(channel),' MUA STRF'])
                print(strjoin([STRF_save_path,string(bat_list(:,bat_num)),...
                    '\Site ',num2str(site),' Chn ',num2str(channel),' STRF'],''),...
                    '-dtiff')
                save(strjoin([STRF_save_path,string(bat_list(:,bat_num)),...
                    '\Site ',num2str(site),' Chn ',num2str(channel),' STA'],''),...
                    'STA')
                save(strjoin([STRF_save_path,string(bat_list(:,bat_num)),...
                    '\Site ',num2str(site),' Chn ',num2str(channel),' faxis'],''),...
                    'faxis')
                save(strjoin([STRF_save_path,string(bat_list(:,bat_num)),...
                    '\Site ',num2str(site),' Chn ',num2str(channel),' taxis'],''),...
                    'taxis')
                %% Generate modulation functions
                MTF = fftshift(fft2(STA));
                % find 0 modulation row and col index
                spec0ind = ceil((size(MTF,1)+1)/2);
                temp0ind = ceil((size(MTF,2)+1)/2);
                % only use positive spectral modulation
                MTF = MTF(spec0ind:end,:);
                % find axis values
                Xrange = diff(X(1:2))*length(X);
                trange = diff(taxis(1:2))*length(taxis);
                specmod = (0:size(MTF,1)-1)/Xrange; % cycles per octave
                tempmod = (-(temp0ind-1):(temp0ind-1))/(trange/1000); % Hz
                %% Generate modulation rate figures
                figure;
                fig1 = surface(tempmod,specmod,abs(MTF));
                set(fig1, 'FaceColor','interp','EdgeColor','interp');
                xlim([-100 100])
                ylim([0 10])
                ylim([0 5])
                colorbar
                colormap jet
                title(['Site ',num2str(site), ' Chn ',num2str(channel),' MUA MTF'])
                print(strjoin([MTF_save_path,string(bat_list(:,bat_num)),...
                    '\Site ',num2str(site),' Chn ',num2str(channel),' MTF'],''),...
                    '-dtiff')
                save(strjoin([MTF_save_path,string(bat_list(:,bat_num)),...
                    '\Site ',num2str(site),' Chn ',num2str(channel),' MTF'],''),...
                    'MTF')
                save(strjoin([MTF_save_path,string(bat_list(:,bat_num)),...
                    '\Site ',num2str(site),' Chn ',num2str(channel),' specmod'],''),...
                    'specmod')
                save(strjoin([MTF_save_path,string(bat_list(:,bat_num)),...
                    '\Site ',num2str(site),' Chn ',num2str(channel),' tempmod'],''),...
                    'tempmod')
                close all
            end
            clearvars -except bat_list data_location sr b a cmap...
                channels_available channel STRF_save_path MTF_save_path...
                ts marker shank_pos d site bat_num cutoff
        end

    end
end