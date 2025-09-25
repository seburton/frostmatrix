clear
%create timeseries sections using hour-long data segments
%and apply matched filtering (waveform and analytic envelope)


%will run faster when set to 'false'
make_plots=true;

%save output for later analysis
save_results=false;

%pause between plots if desired
plot_pause=true;

%for lava lake focus; change later for any LL-distal events
range_thresh = 3;

%load templates
%eval('!/bin/ls Templates/template_*.mat > template_list.txt')
A=importdata('template_list.txt');
ntemplates=length(A);

for i=1:ntemplates
    %template_name1='Templates/template_345_040000_1380_1400.mat';
    template_name_multi{i}=char(A{i});
    eval(['load ',char(template_name_multi{i})]);
    template_multi{:,:,i}=template;
    template_sachdr_multi{i}=template_sachdr;
    template_range_multi{i}=template_range;
end
disp([num2str(ntemplates),' templates loaded'])
%
%genrate random phase templates
for i=1:ntemplates
    tmp_template=cell2mat(template_multi(:,:,i));
    tmp_template_rand=zeros(size(tmp_template));
    [m,n]=size(tmp_template);
    for j=1:n
        tmp_template_rand(:,j)=randomphase(tmp_template(:,j));
    end
    template_multi_rand{:,:,i}=tmp_template_rand;
end

fs=200;
%waveform bandpass
f1=1;
f2=20;
[b,a]=butter(4,[f1/(fs/2),f2/(fs/2)]);
%apply all three components
compstr={'ELZ','ELN','ELE'};
%convert counts to velocity within the instrument passband with this constant
sensitivity=6.02383e8;

if make_plots
    Hf=figure(1);
end

%DAY LOOP
%for day=345:362 %run all data days
for day=345 % single day test
    daystr=num2str(day);

    %master index for storing output for this day
    kindex=1;

    %HOUR LOOP
    for hour=0:23
    %for hour=4:4 %single hour test
        hourstr=num2str(hour);
        if hour < 10
            hourstr=['0',hourstr,'0000'];
        else
            hourstr=[hourstr,'0000'];
        end

        %station='ETS53';
        awkcmd='''{print $2}''';
        %create a list of complete data files for the desired day by doing an ls of
        %the relevant data directory for files

        [ax,bx]=system(['/bin/ls Data_Tomo_Erebus/Data_',daystr,'/Summit_Plateau/*ETS*',hourstr,'.SAC']);
        if  ~contains(bx,'No ')
            cmd=['!/bin/ls Data_Tomo_Erebus/Data_',daystr,'/Summit_Plateau/*ETS*',hourstr,'.SAC | awk -F. ',awkcmd,' | sort | uniq > Station_Lists/station_list_',daystr,'.txt'];


            eval(cmd);

            station_list=textread(['Station_Lists/station_list_',daystr,'.txt'],'%s');

            nstas=length(station_list);
            %three loops to read in each of the three seismogram components for each
            j=0;
            dlen=60*60*fs;

            %load data files (ensure path is correct here)
            %ELZ
            for i=1:nstas
                fname=['Data_Tomo_Erebus/Data_',daystr,'/Y4.',char(station_list(i)),'..',char(compstr(1)),'.M.2008.',daystr,'.',hourstr,'.SAC'];
                if isfile(fname)
                    j=j+1;
                    [sachdr(j),d]=load_sac(fname);
                    data(1:length(d),j)=d;
                    data(:,j)=filter(b,a,data(:,j));
                    data(:,j)=detrend(data(:,j))/sensitivity;

                    if length(data(:,j))<720001
                        if j==1
                            data=[data(:,j);zeros(720001-length(data(:,j)),1)];
                        else
                            data(:,j)=[data(:,j);zeros(720001-length(data(:,j)),1)];
                        end
                    end

                    sd(j)=mad(data(:,j));
                    %ndata(:,j)=data(:,j)/sd(j);
                    ndata(:,j)=data(:,j);
                end
            end

            %ELN
            for i=nstas+1:2*nstas
                fname=['Data_Tomo_Erebus/Data_',daystr,'/Y4.',char(station_list(i-nstas)),'..',char(compstr(2)),'.M.2008.',daystr,'.',hourstr,'.SAC'];
                if isfile(fname)
                    j=j+1;
                    [sachdr(j),d]=load_sac(fname);
                    data(1:length(d),j)=d;
                    data(:,j)=filter(b,a,data(:,j));
                    data(:,j)=detrend(data(:,j))/sensitivity;
                    sd(j)=mad(data(:,j));
                    %ndata(:,j)=data(:,j)/sd(j);
                    ndata(:,j)=data(:,j);
                end
            end

            %ELE
            for i=2*nstas+1:3*nstas
                fname=['Data_Tomo_Erebus/Data_',daystr,'/Y4.',char(station_list(i-2*nstas)),'..',char(compstr(3)),'.M.2008.',daystr,'.',hourstr,'.SAC'];
                if isfile(fname)
                    j=j+1;
                    [sachdr(j),d]=load_sac(fname);
                    data(1:length(d),j)=d;
                    data(:,j)=filter(b,a,data(:,j));
                    data(:,j)=detrend(data(:,j))/sensitivity;
                    sd(j)=mad(data(:,j));
                    %ndata(:,j)=data(:,j)/sd(j);
                    ndata(:,j)=data(:,j);
                end
            end
            %Nominal Lava Lake Location
            LL_lat=-77.5274;
            LL_lon=167.1645;
            LL_elev=3572;
            %get lava lake distances, angles for all stations
            nchans=length(sachdr);
            for i=1:nchans
                [azi(i,1),bazi(i,1),range(i,1),angle(i,1)] =  edist(sachdr(i).stla,sachdr(i).stlo,LL_lat,LL_lon);
            end

            t=(0:length(data)-1)/fs;


            nsamps=length(ndata);

            %counter for found, non-dead, channels, and optional range variable (km)


            %TEMPLATE LOOP
            for itemp=1:ntemplates
                template=cell2mat(template_multi(:,:,itemp));
                template_rand=cell2mat(template_multi_rand(:,:,itemp));
                template_name=char(template_name_multi{itemp});
                template_sachdr=template_sachdr_multi{itemp};
                template_range=template_range_multi{itemp};
                kstore=[0,0];
                k=1;
                %lnorm=sqrt(length(template)-1);

                %original matched filter flow (using the conv function with
                %time reversal to correlate; commented out here)
                for i=1:nchans
                    %must match station and component here betwen the working data set and
                    %the template to be sure that like channels are being correlated (and
                    %nonmatching ones are not used)
                    j = find(strcmp({template_sachdr.kstnm},sachdr(i).kstnm) & strcmp({template_sachdr.kcmpnm},sachdr(i).kcmpnm));
                    %if we have found a match in the data, and we are within template range
                    %(km), then correlate this channel
                    if ~isempty(j) && template_range(j) <= range_thresh

                        %implement moving-window correlation with template channel j and
                        %data channel i
                        %matched_output(:,k)=conv(ndata(:,i),flipud(template(:,j))/(std(template(:,j))*lnorm),'same');

                        %this (moving averaged rms) function provides the proper normalizaton so that the output is
                        %limited to [-1,1]
                        %dnorm(:,k)=sqrt(conv(ndata(:,i).^2,ones(size(template(:,j))),'same'));

                        %save the i and j match indices
                        kstore(k,:)=[i,j];
                        %disp([template_sachdr(j).kstnm,' ',template_sachdr(j).kcmpnm,':',sachdr(i).kstnm,' ',sachdr(i).kcmpnm])

                        k=k+1;
                    end
                end
                nk=length(kstore);

                %use SEC_C (fast correlation) function (Senobari et al.,
                %2019)
                clear D T TR
                D(:,1,:)=ndata(:,kstore(:,1));
                T(:,1,:,:)=template(:,kstore(:,2),1);
                TR(:,1,:,:)=template_rand(:,kstore(:,2),1);

                %moving normalization function-=
                ccc_sum=SEC_C(D,T,2^13,ones(nk,1),ones(nk,1),1);

                %randomized template phase correlation
                ccc_sum_rand=SEC_C(D,TR,2^13,ones(nk,1),ones(nk,1),1);

                %template correlation with data
                matched_output_med=[zeros(length(template)-1,1);ccc_sum]/nk;

                %random phase template correlation with data
                matched_output_med_rand=[zeros(length(template_rand)-1,1);ccc_sum_rand]/nk;

                %find isolated peaks
                threshr=10*std(matched_output_med_rand);

                %minimum peak distance (10 s)
                mpd=20*200/2;
                [PKS_c,LOCS_c]=findpeaks(matched_output_med,'MinPeakHeight',threshr,'MinPeakDistance',mpd,'MinPeakProminence',0.02);

                %show the results for each template correlation
                if make_plots
                set(0,'CurrentFigure',Hf);
                    clf
                    hold on
                    plotdec=10;
                    for i=1:nchans
                        %quick plot
                        plot(t(1:plotdec:end),ndata(1:plotdec:end,i)/max(ndata(1:plotdec:end,i))+range(i))
                    end

                    %templates were created with 1 s inital offset, shift detection by half
                    %window length as well
                    %d_shift = length(template)/2/fs-1;

                    d_shift=length(template-1)/fs;

                    plot(t-d_shift,matched_output_med-1,'k-','linewidth',3)
                    plot(t(LOCS_c)-d_shift,matched_output_med(LOCS_c)-1,'r*','Markersize',20)
                    ylabel('km, median channel matched filter output')
                    xlabel('Time (Seconds in Hour)')
                    xlim([0 3600])
                    title(['Detections ',daystr,' ',hourstr,' ',template_name],'interpreter','none')
                    bookfonts_TNR
                    hold off
                end

                %various waveform correlation results for inspection
                % figure(2)
                % clf
                % subplot(3,1,1)
                % plot(t,matched_output_mean)
                % ylabel('Mean')
                % xlim([0 60])
                % bookfonts_TNR
                %
                % subplot(3,1,2)
                % plot(t,matched_output_med)
                % hold on
                % plot(t(LOCS_c),matched_output_med(LOCS_c),'r*')
                % hold off
                % ylabel('Median')
                % xlim([0 60])
                % bookfonts_TNR
                %
                % subplot(3,1,3)
                % plot(t,matched_output_max)
                % xlabel('Time (s)')
                % ylabel('Max')
                % xlim([0 60])
                % bookfonts_TNR

                %event cutting window length (samples)
                cutw=20*fs;
                if make_plots
                    set(0,'CurrentFigure',Hf);
                    hold on
                    %event bracketed detections (green for event start and
                    %red for event end
                    if ~isempty(LOCS_c)
                        for i=1:length(LOCS_c)
                            plot([t(LOCS_c(i))-cutw/fs,t(LOCS_c(i))-cutw/fs],[3,0],'g-','linewidth',2);
                            plot([t(LOCS_c(i))+2*cutw/fs,t(LOCS_c(i))+2*cutw/fs],[3,0],'r-.','linewidth',2);
                        end
                    end

                    hold off
                    ylim([-1.1 3])
                    xlim([0 3600])
                end

                %save matched filter median outputs across all templates here
                template_matched_output.med_output{kindex}=single(matched_output_med);
                template_matched_output.day(kindex)=day;
                template_matched_output.hour(kindex)=hour;
                template_detects.PKS{kindex}=PKS_c;
                template_detects.LOCS{kindex}=LOCS_c;
                template_detects.mean_cc{kindex}=mean(matched_output_med);
                template_detects.std_cc{kindex}=std(matched_output_med);
                template_detects.med_cc{kindex}=median(matched_output_med);
                template_detects.mad_cc{kindex}=mad(matched_output_med);
                template_detects.nchan(kindex)=nk;
                template_detects.day(kindex)=day;
                template_detects.hour(kindex)=hour;
                disp([template_name ' ',num2str(day),' ',num2str(hour),...
                    ' (',num2str(length(template_detects.PKS{kindex})),' detections, max: ',num2str(max(template_detects.PKS{kindex})),')'])

                kindex=kindex+1;
                if plot_pause
                    disp('Pausing for figure')
                    pause(0.1)
                end
            end

            %save detection time series and metadata
            %eval(['save detects_',template_name(1:end-4),'_day_',daystr,'_hour_',hourstr(1:2),...
            %    ' matched_output_med', ...
            %    ' matched_output_med_h t day hour fs template_name sachdr range azi bazi ', ...
            %    ' PKS_c LOCS_c PKS_h LOCS_h']);
            hourstr;
            daystr;
            %pause
            %[~,ndata]=size(data);
            stnames={sachdr.kstnm};
            if save_results
                eval(['save Test_Results/erupt_detects_',daystr,'.mat template_detects template_name_multi ndata stnames'])
                eval(['save Test_Results/erupt_matched_output_',daystr,'.mat template_matched_output template_name_multi'])
            end
            clear sachdr data
        else
            disp('No data files found')
        end

    end %end day loop

end %end day loop