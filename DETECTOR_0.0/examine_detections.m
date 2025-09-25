clear

ts_length=17*200*3600*24; %length of 17-day time series in 200 s/s samples
dt=1/(24*3600*200); %one sample in units of days for datetime variable time step
t0=datetime(2008,12,10); %start of data (day 345 of 2008)
td=(t0:dt:t0+17-dt)'; %datetime time base at full sampling rate
dets=sparse(zeros(ts_length,1)); %store detection peaks from see_med_ouptut.m here
evs=sparse(zeros(ts_length,1)); %store event times here (as impulses)

%create a full detection time series from output of see_med_output.m
load Test_Results/med_outputs.mat

%detection spikes (bartlett function convolved detections from
%see_med_output.m
zd_all=sparse([]);
%proportion of templates with detections
zn_all=sparse([]);

%sample index for populating detections (dets)
for day=345 %:361
    it0=(day-345)*24*200*3600;
    for nhour=1:24
    %add multihour shifts to to the index
    it1=it0+(nhour-1)*200*3600;
    nday=day-344;
    zd_day=detect_ts.zd{nday,nhour};
    zn_day=detect_ts.zn{nday,nhour};
    zd_all=[zd_all;zd_day(1:end-1)];
    zn_all=[zn_all;zn_day(1:end-1)];
    %Detection peak times are logged in terms of seconds since start of the relevant day
    dets(it1+detect_peaks.L{nday,nhour})=detect_peaks.P{nday,nhour};
    evs(it1+detect_peaks.L{nday,nhour})=1;
    %nhour
    end
    day;
end

%used number of templates for each event
ntemplates=sparse(zeros(ts_length,1));
ntemplates(evs>0)=zn_all(evs>0)*17;

%return

dets_per_min=full(movsum(evs,200*60));
dets_per_hour=full(movsum(evs,200*3600));
dets_per_day=full(movsum(evs,200*3600*24));

A=importdata('template_list.txt');
for i=1:17
    a=char(A(i));
    tday(i)=str2num(a(20:22));
    tmday(i)=10+tday(i)-345;
    thour(i)=str2num(a(24:25));
    tsec(i)=str2num(a(31:34));
    temptd(i)=datetime(2008,12,tmday(i),thour(i),0,tsec(i));
end

tde=td(evs>0);
detse=dets(evs>0);

ind=find(detse>3);
for i=1:17
    [~,ix]=min(abs(tde(ind)-temptd(i)));
    itdetemp(i)=ind(ix);
end

figure(10)
clf

c=jet(17);
colormap(c);
cvals=ntemplates(evs>0);
subplot(4,1,1)
scatter(tde,log10(detse),100*detse,cvals);
ylabel("log_{10} Detector Output")
hc=colorbar;
hc.Location='east';
hold on
plot(tde(itdetemp),log10(detse(itdetemp)),'k.','markersize',50);
hold off
yyaxis right
set(gca, 'YTick', [])
ylabel('N_{Templates}','color','k')
bookfonts_TNR

% subplot(4,1,1)
% plot(td(evs>0),dets(evs>0),'bo','markersize',5)
% ylabel("Detector Output")
% bookfonts_TNR

subplot(4,1,2)
plot(td(1:200:end),smooth(dets_per_hour(1:200:end)),'k-','linewidth',2);
ylabel("Events per Hour")
bookfonts_TNR

subplot(4,1,3)
plot(td(1:200:end),smooth(dets_per_day(1:200:end)),'k-','linewidth',2);
ylabel("Events per Day")
bookfonts_TNR

cs=cumsum(evs);
subplot(4,1,4)
plot(td(1:200:end),cs(1:200:end),'k-','linewidth',2);
ylabel("Cum. Events")
bookfonts_TNR

tdmin=min(td);
tdmax=max(td);
save workspace_examine_detections.mat evs dets tde detse tdmin tdmax dt itdetemp