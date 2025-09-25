%Create a template for match filtering

clear

dpath='/Volumes/Aster_WD_1/New_Home/Projects_new/Erebus/Erebus_2020/Data_Tomo_Erebus_orig';
A=importdata('../template_list.txt');
ntemplates=length(A);
for i=1:ntemplates
    tname=char(A{i});

daystr=tname(20:22); 
hourstr=tname(24:29);
secsinhour1=str2double(tname(31:34));
secsinhour2=str2double(tname(36:39));

fs=200;
%waveform bandpass
f1=1;
f2=10;
[b,a]=butter(4,[f1/(fs/2),f2/(fs/2)]);
%apply all three components
compstr={'ELZ','ELN','ELE'};


% hourstr=num2str(hour);
% if hour < 10
%     hourstr=['0',hourstr,'0000'];
% else
%     hourstr=[hourstr,'0000'];
% end
sensitivity=6.02383e8;
%station='ETS53';
awkcmd='''{print $2}''';
%create a list of complete data files for the desired day by doing an ls of
%the relevant data directory for files
cmd=['!/bin/ls ',dpath,'/Data_',daystr,'/*ETS*',hourstr,'.SAC | awk -F. ',awkcmd,' | sort | uniq > station_list_',daystr,'.txt'];
eval(cmd);

station_list=textread(['station_list_',daystr,'.txt'],'%s');

%[sachdr,data]=load_sac(['Data_345/Y4.',station,'..ELZ.M.2008.345.SAC']);
%Y4.ETS53..ELZ.M.2008.345.030000.SAC
%[sachdr,data]=load_sac(['Data_345/Y4.ETS42..ELZ.M.2008.345.SAC']);
nstas=length(station_list);
%three loops to read in each of the three seismogram components for each
j=0;
dlen=60*60*fs;

for i=1:nstas
    fname=[dpath,'/Data_',daystr,'/Y4.',char(station_list(i)),'..',char(compstr(1)),'.M.2008.',daystr,'.',hourstr,'.SAC'];
    if isfile(fname)
    j=j+1;
[sachdr(j),d]=load_sac(fname);
data(1:length(d),j)=d;
data=data-mean(data);
data(:,j)=filter(b,a,data(:,j));
data(:,j)=detrend(data(:,j))/sensitivity;
sd(j)=mad(data(:,j));
ndata(:,j)=data(:,j)/sd(j);
    end
end
% 
for i=nstas+1:2*nstas
    fname=[dpath,'/Data_',daystr,'/Y4.',char(station_list(i-nstas)),'..',char(compstr(2)),'.M.2008.',daystr,'.',hourstr,'.SAC'];
    if isfile(fname)
    j=j+1;
[sachdr(j),d]=load_sac(fname);
data(1:length(d),j)=d;
data=data-mean(data);
data(:,j)=filter(b,a,data(:,j));
data(:,j)=detrend(data(:,j))/sensitivity;
sd(j)=mad(data(:,j));
ndata(:,j)=data(:,j)/sd(j);
    end
end

for i=2*nstas+1:3*nstas
    fname=[dpath,'/Data_',daystr,'/Y4.',char(station_list(i-2*nstas)),'..',char(compstr(3)),'.M.2008.',daystr,'.',hourstr,'.SAC'];
    if isfile(fname)
    j=j+1;
[sachdr(j),d]=load_sac(fname);
data(1:length(d),j)=d;
data=data-mean(data);
data(:,j)=filter(b,a,data(:,j));
data(:,j)=detrend(data(:,j))/sensitivity;
%sd(j)=mad(data(:,j));
%ndata(:,j)=data(:,j)/sd(j);
    end
end

    %Nominal Lava Lake Location for station range
    LL_lat=-77.5274;
	LL_lon=167.1645;
	LL_elev=3572;
nchans=length(sachdr);
for i=1:nchans
    [azi(i,1),bazi(i,1),range(i,1),angle(i,1)] = edist(sachdr(i).stla,sachdr(i).stlo,LL_lat,LL_lon);
end

%for i=1:nchans
    %ndata_h(:,i)=abs(hilbert(ndata(:,i)));
    %ndata_h(:,i)=ndata_h(:,i)-mean(ndata_h(:,i));
%end
    
%extract a template

x1=secsinhour1*fs; 
x2=secsinhour2*fs; 
%length of template "cut"
for i=1:nchans
    template(:,i)=data(x1:x2,i);
    %normalize template energy
    %template(:,i)=template(:,i)/norm(template(:,i));
    %generate (de-meaned)envelope functions (analytic function magnitudes)
    %for additional matched filtering testing
    %template_h(:,i)=ndata_h(x1:x2,i);
    %template_h(:,i)=template_h(:,i)/norm(template_h(:,i));
end

%remove dead channels
sdh=std(template);
istore=find(sdh>mean(sdh)/1000);
template_sachdr=sachdr(istore);
template_range=range(istore);
template_nchans=numel(istore);
template=template(:,istore);
%template_h=template_h(:,istore);

%save a template and indexed sachdr information here for future use
%eval(['save template_',daystr,'_',hourstr,'_',sprintf('%4.4d',secsinhour1),'_',sprintf('%4.4d',secsinhour2),'.mat template_sachdr template template_h x1 x2 template_range template_nchans'])
eval(['save template_',daystr,'_',hourstr,'_',sprintf('%4.4d',secsinhour1),'_',sprintf('%4.4d',secsinhour2),'.mat template_sachdr template x1 x2 template_range template_nchans'])
end
