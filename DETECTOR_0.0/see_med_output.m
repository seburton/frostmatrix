clear
showfigs=true;
%width of detection window for gathering multi-template detections
%(samples)
detect_width=200*10;
%requirement for templates with consistent detections
%can be highgraded later using information in the saved .mat file at the end of this script
template_num_thresh=1/17;
%convolution kernel for estimating best detection times
detect_kernel=bartlett(detect_width);

%These are the intertemplate time shifts in seconds relative to template 17 estimated by
%compare_templates.m
Dshift=[2.1050
    -0.0250
    -0.6700
    0.1500
    -1.1200
    -0.6200
    -0.5550
    0.6750
    0.1350
    0.0150
    -0.9250
    2.1150
    -0.2650
    0.0200
    -0.7600
    -0.7150
    0];
%one hour time vector for plotting each segment (units of seconds)
tz=(0:60*60*200);

%days=345:361;
days=345; %single day test
%days=362;
for day=days
    %day index
    nday=day-min(days)+1;
    figure(nday)
    clf
    %load matched output files for this day as created by timeseries_run_fast_dev.m
    eval(['load Test_Results/erupt_matched_output_',num2str(day),'.mat']);
    %data structure
    % template_matched_output      1x1             1175097537  struct
    %  med_output: {1×408 cell}
    %  day: [345 345 345 345 345 345 345 345 345 345 345 345 345 345 345 345 345 345 345 345 345 … ] (1×408 double)
    %  hour: [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 … ] (1×408 double)
    %  template_name_multi          1x17                  3502  cell
    %load detection outputs from timeseries_run_fast_dev.m
    eval(['load Test_Results/erupt_detects_',num2str(day),'.mat']);
    %data structure
    % ndata                    1x1                   8  double
    % stnames                  1x159             21624  cell

    % template_detects         1x1              340060  struct
    %     PKS: {1×408 cell}
    %    LOCS: {1×408 cell}
    % mean_cc: {1×408 cell}
    %  std_cc: {1×408 cell}
    %  med_cc: {1×408 cell}
    %  mad_cc: {1×408 cell}
    %   nchan: [131 131 131 131 128 128 128 128 128 128 128 128 128 128 128 128 128 132 132 132 132 129 … ] (1×408 double)
    %     day: [345 345 345 345 345 345 345 345 345 345 345 345 345 345 345 345 345 345 345 345 345 345 … ] (1×408 double)
    %    hour: [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 … ] (1×408 double)

    % template_name_multi      1x17                 3502  cell
    ntemplates=length(template_name_multi);

    for hour=0:23
        nhour=hour+1;
        tzd0=datetime(2008,12,8+nday-1,hour,0,0);
        tzd=(tzd0:1/(24*60*60*200):tzd0+1/24)';
        zdetect=zeros(60*60*200+1,17);

        for i=1:ntemplates
            x(hour+1,:,i)=template_matched_output.med_output{i+ntemplates*hour};
        end
        xsum=sum(squeeze(x(hour+1,:,:)),2);
        %plot all detections with distingushible colors fo reach template
        c=distinguishable_colors(ntemplates);

        subplot(6,4,hour+1)

        idx=find(xsum>0);
        %plot(t(idx),xsum(idx));
        hold on
        %loop over templates for this day and hour
        for i=1:ntemplates
            %plot detections (with Dshift incorporated; see above)
            tdetects_plot_dtime=(template_detects.LOCS{i+ntemplates*hour}/200+Dshift(i))/(24*3600)+tzd0;
            plot(tdetects_plot_dtime,template_detects.PKS{i+ntemplates*hour},...
                '.','color',c(i,:),'markersize',20)
            %save the Dshift-corrected detector outputs for each template here
            zdetect(template_detects.LOCS{i+ntemplates*hour}+Dshift(i)*200,i)=sparse(template_detects.PKS{i+ntemplates*hour});
        end
        hold off
        xlim([tzd0 tzd0+1/24])
        xlabel('Seconds')
        ylabel('Detection')
        title(['Day: ',num2str(day),' Hour: ',num2str(hour)])
        bookfonts_TNR

        hold on
        subplot(6,4,hour+1)
        %triangular kernal to convolve the sum of the detector outputs across
        %templates

        if length(zdetect) > 60*60*200+1
            zdetect=zdetect(1:60*60*200+1,:);
        end
        %get moving point count,zn
        zcount=zdetect;
        zcount(zcount>0)=1;
        %moving normalized number of detections across detect_width (higher
        %value indicates that more templates had positive detections in timeseries_run_fast_dev.m)
        zn=sparse(movsum(sum(zcount,2),detect_width)/ntemplates);
        %detection function using detect_kernel
        zd=sparse(conv(sum(zdetect,2),detect_kernel,'same'));
        %plot detection function (higher value indicates higher correlation sums convolved across detect_kernel)

        plot(tzd,full(zd));
        %find detection peaks using this multi-template averaging kernel
        [P,L]=findpeaks(zd);
        indsel=find(zn(L)>=template_num_thresh);
        L=L(indsel);
        P=P(indsel);
        %plot templates contributing (normalized to 1 for all templates)
        plot(tzd(L),full(zn(L)),'k*');
        %save events detected
        %ev_count(nday,nhour)=numel(L);
        detect_peaks.P{nday,nhour}=P;
        detect_peaks.L{nday,nhour}=L;
        %save output detection time series
        detect_ts.zn{nday,nhour}=zn;
        detect_ts.zd{nday,nhour}=zd;
        %xlim([0 60*60])
        xlabel('Seconds')
        title(['Day: ',num2str(day),' Hour: ',num2str(hour),' N=',num2str(numel(L))])
        bookfonts_TNR
        hold off

        %end hour loop
    end
    %end day loop
disp(['Found ',num2str(sum(cellfun(@numel,{detect_peaks.L{nday,:}}))),' detections for day ',num2str(day)])
end
save Test_Results/med_outputs.mat detect_peaks detect_ts detect_kernel template_num_thresh;

%save all data from a multiday run of this program (often run with showfigs
%set to false
%save med_outputs_345_361.mat detect_kernel detect_peaks detect_ts template_num_thresh 

