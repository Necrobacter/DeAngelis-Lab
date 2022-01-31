%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = '\\casfsb\biology_labs\Finzi Lab\Biogeochemistry\All Lab Data\Marc\Adrien\Harvard Forest dataset\Processed data\Manuscript - reviews by co-authors\Data for Harvard Forest Archive\Data files';


%%%%%%%%%%%%%%%%%%
% Data use comes from the Harvard Forest Archive, but files were often merged/splitted/slightly modified to make it easier to load in Matlab
%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fisher Meteorological Station
[a1,b1,raw] = xlsread([path,'\Fisher_2001-2004.xls'],'A3:AD10000');
[a2,b2,raw] = xlsread([path,'\Fisher_2001-2004.xls'],'A10001:AD20000');
[a3,b3,raw] = xlsread([path,'\Fisher_2001-2004.xls'],'A20001:AD34082');
data1_num = [a1 ; a2 ; a3];
data1_txt = [b1 ; b2 ; b3];
clear a1 a2 a3 b1 b2 b3;

Fisher_year1 = data1_num(:,1);
Fisher_jday1 = data1_num(:,2);
Fisher_hour1 = repmat([0:23]',1421,1);
Fisher_hour1 = Fisher_hour1(2:34081,1);

Fisher_year1_hh([1:2:length(Fisher_year1)*2],1) = Fisher_year1;
Fisher_year1_hh([2:2:length(Fisher_year1)*2],1) = Fisher_year1;
Fisher_year1_hh = [2001 ; Fisher_year1_hh(1:end-1)];
Fisher_jday1_hh([1:2:length(Fisher_jday1)*2],1) = Fisher_jday1;
Fisher_jday1_hh([2:2:length(Fisher_jday1)*2],1) = Fisher_jday1;
Fisher_jday1_hh = [42 ; Fisher_jday1_hh(1:end-1)];
Fisher_hour1_hh([1:2:length(Fisher_hour1)*2],1) = Fisher_hour1;
Fisher_hour1_hh([2:2:length(Fisher_hour1)*2],1) = Fisher_hour1;
Fisher_hour1_hh = [0 ; Fisher_hour1_hh(1:end-1)];
Fisher_minute1_hh = zeros(size(Fisher_hour1_hh));
Fisher_minute1_hh(1:2:length(Fisher_minute1_hh),1) = 30;

temp = data1_num(:,29);
Fisher_Ts_10cm1([1:2:length(temp)*2],1) = temp;
Fisher_Ts_10cm1([2:2:length(temp)*2],1) = temp;

Fisher_Ta1 = data1_num(:,3);
Fisher_PAR1 = data1_num(:,13);

%%%%%
[a1,b1,raw] = xlsread([path,'\Fisher_2005-2011.xlsx'],'A3:AD20000');
[a2,b2,raw] = xlsread([path,'\Fisher_2005-2011.xlsx'],'A20001:AD40000');
[a3,b3,raw] = xlsread([path,'\Fisher_2005-2011.xlsx'],'A40001:AD60000');
[a4,b4,raw] = xlsread([path,'\Fisher_2005-2011.xlsx'],'A60001:AD80000');
[a5,b5,raw] = xlsread([path,'\Fisher_2005-2011.xlsx'],'A80001:AD100000');
[a6,b6,raw] = xlsread([path,'\Fisher_2005-2011.xlsx'],'A100001:AD120000');
[a7,b7,raw] = xlsread([path,'\Fisher_2005-2011.xlsx'],'A120001:AD140000');
[a8,b8,raw] = xlsread([path,'\Fisher_2005-2011.xlsx'],'A140001:AD160000');
[a9,b9,raw] = xlsread([path,'\Fisher_2005-2011.xlsx'],'A160001:AD175298');
data2_num = [a1;a2;a3;a4;a5;a6;a7;a8;a9];
data2_txt = [b1;b2;b3;b4;b5;b6;b7;b8;b9];
clear a1 a2 a3 a4 a5 a6 a7 a8 a9 b1 b2 b3 b4 b5 b6 b7 b8 b9;
Fisher_year2 = data2_num(:,1);
Fisher_jday2 = data2_num(:,2);
temp = repmat([0:23]',1827,1);
temp = temp(1:43825,1);
Fisher_hour2([1:4:175300],1) = temp;
Fisher_hour2([2:4:175300],1) = temp;
Fisher_hour2([3:4:175300],1) = temp;
Fisher_hour2([4:4:175300],1) = temp;
Fisher_hour2 = Fisher_hour2(2:175297);

Fisher_year2_hh = Fisher_year2([2:2:length(Fisher_year2)],1);
Fisher_jday2_hh = Fisher_jday2([2:2:length(Fisher_jday2)],1);
Fisher_hour2_hh = Fisher_hour2([2:2:length(Fisher_hour2)],1);
Fisher_minute2_hh = zeros(size(Fisher_hour2_hh));
Fisher_minute2_hh(1:2:length(Fisher_minute2_hh),1) = 30;

temp = data2_num(:,29);
y=0;
for x=1:2:length(temp)
    y=y+1;
    Fisher_Ts_10cm2(y,1) = mynanmean(temp(x:x+1,1));
end

Fisher_Ta2 = zeros([length(data2_num)/4],1)*NaN;
Fisher_PAR2 = zeros([length(data2_num)/4],1)*NaN;
ind = find(isnan(data2_num(:,9)));
data2_num(ind,9) = 0;
y=0;
for x=1:length(Fisher_Ta2)
    y=y+1;
    Fisher_Ta2(y,1) = mean(data2_num([x*4-3:x*4],3));
    Fisher_PAR2(y,1) = mean(data2_num([x*4-3:x*4],13));
end

Fisher_timestamps_hh = [[Fisher_year1_hh ; Fisher_year2_hh] [Fisher_jday1_hh ; Fisher_jday2_hh] [Fisher_hour1_hh ; Fisher_hour2_hh] [Fisher_minute1_hh ; Fisher_minute2_hh]];
Fisher_Ts_10cm_hh = [Fisher_Ts_10cm1 ; Fisher_Ts_10cm2];
Fisher_Ta = [Fisher_Ta1 ; Fisher_Ta2];
Fisher_PAR = [Fisher_PAR1 ; Fisher_PAR2];
Fisher_timestamps_h = Fisher_timestamps_hh(2:2:end,1:3);

%Gap-fill Fisher Tower Tair (Ta) data
Fisher_Ta = FillSmallGapsByLinInterp(Fisher_Ta,6); %Gap-fill short gaps by linear interpolation


clear Fisher_year1_hh Fisher_year2_hh Fisher_jday1_hh Fisher_jday2_hh Fisher_hour1_hh Fisher_hour2_hh Fisher_minute1_hh Fisher_minute2_hh ...
    Fisher_year1 Fisher_year2 Fisher_jday1 Fisher_jday2 Fisher_hour1 Fisher_hour2 Fisher_minute1 Fisher_minute2 ...
    Fisher_Ts_10cm1 Fisher_Ts_10cm2 temp data1_num data1_txt raw data2_num data2_txt x y Fisher_PAR1 Fisher_PAR2 Fisher_Ta1 Fisher_Ta2;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load climate and eddy covariance data
%EMS Tower
data_filled = csvread([path,'\EMS - filled data2.csv'],2,0);

EMS_year = data_filled(:,1);
EMS_jday = data_filled(:,3);
EMS_hour = data_filled(:,4);
EMS_filled_Ta_27m_GF = data_filled(:,17);
EMS_filled_Ta_27m_GF([89528:89529 89531:89533 90388],1) = NaN; %Delete bad data
EMS_filled_Ta_27m_GF = FillSmallGapsByLinInterp(EMS_filled_Ta_27m_GF,4);%Gap-fill short gaps
EMS_filled_PAR_28m_GF = data_filled(:,29);

EMS_NEE = data_filled(:,7);
EMS_NEE_GF = data_filled(:,11);
EMS_R_GF = data_filled(:,12);
EMS_GEE_GF = data_filled(:,13);

clear data_filled

%Compute monthly sums from original (PI) gap-filling
NEE3=EMS_NEE_GF.*12.*3600./1000000;
GEE3=EMS_GEE_GF.*12.*3600./1000000;
R3=EMS_R_GF.*12.*3600./1000000;
y=0;
EMS_monthly_original_GF=[];
for x=EMS_year(1):EMS_year(end)
    for xx=1:12
        y=y+1;
        [dayi,dayj]=fn_month_to_jday(x,xx);
        ind = find(EMS_year==x & EMS_jday>=dayi & EMS_jday<=dayj);
        EMS_monthly_original_GF(y,:) = [x xx sum(NEE3(ind)) sum(R3(ind)) sum(GEE3(ind))];
    end
end
EMS_monthly_original_GF(1:10,:)=[]; %Use for Harvard Forest Archive (EMS-Re, file Re_Rs_monthly.csv)

y=0;
EMS_yearly_original_GF=[];
for x=1992:2009
    y=y+1;
    ind = find(EMS_monthly_original_GF(:,1)==x);
    EMS_yearly_original_GF(y,:) = [x sum(EMS_monthly_original_GF(ind,4))]; %Use for Harvard Forest Archive (EMS-Re, file Re_Rs_annual.csv)
end

y=0;
EMS_daily_original_GF=[];
for x=EMS_year(1):EMS_year(end)
    for xx=1:366
        ind = find(EMS_year==x & EMS_jday==xx);
        if ~isempty(ind)
            y=y+1;
            try
                ind = ind+1;
                EMS_daily_original_GF(y,:) = [x xx sum(NEE3(ind)) sum(R3(ind)) sum(GEE3(ind))];
            catch
                ind = ind(2:end)-1;
                EMS_daily_original_GF(y,:) = [x xx sum(NEE3(ind)) sum(R3(ind)) sum(GEE3(ind))]; %Use for Harvard Forest Archive (column 4 is EMS-Re, file Re_Rs_daily.csv)
            end
        end
    end
end
clear NEE3 GEE3 R3;

%EMS part #2
data_final = csvread([path,'\EMS - final data2.csv'],1,0);
EMS_Ts_20cm = data_final(:,35);
EMS_Ts_20cm([74676 81445 84755 100478 101390 109383 119298:119440 120135 120230 131830 136386 140058 141499],1)=NaN; %Delete bad data
EMS_Ts_20cm(EMS_Ts_20cm>=20)=NaN; %Delete bad data
EMS_Ts_surf = data_final(:,39);
EMS_Ts_surf([74676],1)=NaN; %Delete bad data

EMS_year_hh = zeros(length(EMS_Ts_20cm)*2,1)*NaN;
EMS_year_hh([1:2:length(EMS_Ts_20cm)*2],1) = EMS_year;
EMS_year_hh([2:2:length(EMS_Ts_20cm)*2],1) = EMS_year;
EMS_jday_hh = zeros(length(EMS_Ts_20cm)*2,1)*NaN;
EMS_jday_hh([1:2:length(EMS_Ts_20cm)*2],1) = EMS_jday;
EMS_jday_hh([2:2:length(EMS_Ts_20cm)*2],1) = EMS_jday;
EMS_hour_hh = zeros(length(EMS_Ts_20cm)*2,1)*NaN;
EMS_hour_hh([1:2:length(EMS_Ts_20cm)*2],1) = EMS_hour;
EMS_hour_hh([2:2:length(EMS_Ts_20cm)*2],1) = EMS_hour;
EMS_minute_hh = zeros(length(EMS_Ts_20cm)*2,1);
EMS_minute_hh([2:2:length(EMS_Ts_20cm)*2],1) = 30;
EMS_Ts_20cm_hh = zeros(length(EMS_Ts_20cm)*2,1)*NaN;
EMS_Ts_20cm_hh([1:2:length(EMS_Ts_20cm)*2],1) = EMS_Ts_20cm;
EMS_Ts_20cm_hh([2:2:length(EMS_Ts_20cm)*2],1) = EMS_Ts_20cm;
EMS_Ts_surf_hh = zeros(length(EMS_Ts_20cm)*2,1)*NaN;
EMS_Ts_surf_hh([1:2:length(EMS_Ts_20cm)*2],1) = EMS_Ts_surf;
EMS_Ts_surf_hh([2:2:length(EMS_Ts_20cm)*2],1) = EMS_Ts_surf;

EMS_timestamps = [EMS_year EMS_jday EMS_hour];
EMS_timestamps_hh = [EMS_year_hh EMS_jday_hh EMS_hour_hh EMS_minute_hh];

clear data_final EMS_final_time_days EMS_Ts_20cm EMS_Ts_surf EMS_year_hh EMS_jday_hh EMS_hour_hh EMS_minute_hh ...
    EMS_year EMS_jday EMS_hour EMS_minute




%HEM Tower
data1 = csvread([path,'\HEM Tower - eddy flux 2000-01.csv'],2,0);

%Add missing rows and delete those doubled
data1 = [data1(1:35,:) ; [2000 306.75 ones(1,11)*NaN] ; ...
    data1(36:479,:) ; data1(481:4061,:) ; ...
    [2001 24.63 ones(1,11)*NaN] ; [2001 24.65 ones(1,11)*NaN] ; ...
    data1(4062:4193,:) ; data1(4195:4444,:) ; data1(4446:7293,:) ; ...
    [2001 91.96 data1(7340,3:end)] ; [2001 91.98 data1(7342,3:end)] ; ...
    data1(7294:7339,:) ; data1(7341,:) ; data1(7343:12684,:) ; data1(12686:13807,:) ; ...
    [2001 227.65 ones(1,11)*NaN] ; [2001 227.67 ones(1,11)*NaN] ; ...
    data1(13808:14034,:) ; [2001 232.42 ones(1,11)*NaN] ; data1(14035:17326,:) ; ...
    data1(17328,:) ; data1(17330:end,:)];

%Compute NEE
HEM1_NEE_measured = data1(:,11);
HEM1_NEE_estimated = data1(:,12);
HEM1_NEE_flag = data1(:,13);
HEM1_NEE_temp = HEM1_NEE_flag*NaN;
HEM1_NEE_measured2 = HEM1_NEE_measured*NaN;
ind = find(HEM1_NEE_flag==1);
HEM1_NEE_measured2(ind,1) = HEM1_NEE_measured(ind,1);
HEM1_NEE_temp(ind,1) = HEM1_NEE_measured(ind,1);
ind = find(HEM1_NEE_flag==2);
HEM1_NEE_temp(ind,1) = HEM1_NEE_estimated(ind,1);


%Convert to hourly data
ind1 = 1:2:17542;
data1b = zeros(8771,13)*NaN;
HEM1_NEE_GF = zeros(8771,1)*NaN;
HEM1_NEE = zeros(8771,1)*NaN;
for x=1:8771
    ind = ind1(x);
    data1b(x,:) = mean(data1(ind:ind+1,:),1);
    HEM1_NEE_GF(x,1) = mean(HEM1_NEE_temp(ind:ind+1,:),1);
    HEM1_NEE(x,1) = mean(HEM1_NEE_measured2(ind:ind+1,:),1);
end


Climate_HEM1_year = data1(:,1);
Climate_HEM1_day_time = data1(:,2);
Climate_HEM1_day = floor(Climate_HEM1_day_time);
Climate_HEM1_hour = repmat([sort(repmat([0:23]',2,1))],366,1);
Climate_HEM1_hour = Climate_HEM1_hour(2:17544,1);
Climate_HEM1_year_hh = Climate_HEM1_year;
Climate_HEM1_day_hh = Climate_HEM1_day;
Climate_HEM1_hour_hh = Climate_HEM1_hour;
Climate_HEM1_minute_hh = zeros(17543,1);
Climate_HEM1_minute_hh([1:2:17543],1) = 30;
Climate_HEM1_year = Climate_HEM1_year(ind1+1,:);
Climate_HEM1_day = Climate_HEM1_day(ind1+1,:);
Climate_HEM1_hour = Climate_HEM1_hour(ind1+1,:);

Climate_HEM1_Ts_10cm = data1(:,9);


%HEM part #2
data2 = csvread([path,'\HEM Tower - eddy flux 2004.csv'],2,0);

HEM2_NEE_measured = data2(:,1)*NaN;
HEM2_NEE_temp = data2(:,1)*NaN;
ind = find(data2(:,13)==0);
HEM2_NEE_temp(ind,1) = data2(ind,12);
ind = find(data2(:,13)==1);
HEM2_NEE_temp(ind,1) = data2(ind,11);
HEM2_NEE_measured(ind,1) = data2(ind,11);

%Convert to hourly data
ind1 = 2:2:9422;
data2b = zeros(4711,25)*NaN;
HEM2_NEE = zeros(4711,1)*NaN;
HEM2_NEE_GF = zeros(4711,1)*NaN;
for x=1:4711
    ind = ind1(x);
    data2b(x,:) = mean(data2(ind:ind+1,:),1);
    HEM2_NEE(x,1) = mean(HEM2_NEE_measured(ind:ind+1,:),1);
    HEM2_NEE_GF(x,1) = mean(HEM2_NEE_temp(ind:ind+1,:),1);
end

data2(end,1:2) = [2005 1];
Climate_HEM2_year = data2(:,1);
Climate_HEM2_day_time = data2(:,2);
Climate_HEM2_day = floor(Climate_HEM2_day_time);
Climate_HEM2_hour = repmat([sort(repmat([0:23]',2,1))],198,1);
Climate_HEM2_hour = Climate_HEM2_hour(35:9457,1);
Climate_HEM2_year_hh = Climate_HEM2_year;
Climate_HEM2_day_hh = Climate_HEM2_day;
Climate_HEM2_hour_hh = Climate_HEM2_hour;
Climate_HEM2_minute_hh = zeros(9423,1);
Climate_HEM2_minute_hh([2:2:9423],1) = 30;
Climate_HEM2_year = Climate_HEM2_year(ind1+1,:);
Climate_HEM2_day = Climate_HEM2_day(ind1+1,:);
Climate_HEM2_hour = Climate_HEM2_hour(ind1+1,:);

Climate_HEM2_Ts_10cm = data2(:,24);



%HEM part #3
data3 = csvread([path,'\HEM Tower - eddy flux since 2005.csv'],2,0);

HEM3_NEE_measured = data3(:,1)*NaN;
HEM3_NEE_temp = data3(:,1)*NaN;
ind = find(data3(:,15)==0);
HEM3_NEE_temp(ind,1) = data3(ind,14);
ind = find(data3(:,15)==1);
HEM3_NEE_temp(ind,1) = data3(ind,13);
HEM3_NEE_measured(ind,1) = data3(ind,13);

%Convert to hourly data
ind1 = 1:2:87648;
data3b = zeros(43824,26)*NaN;
HEM3_NEE = zeros(43824,1)*NaN;
HEM3_NEE_GF = zeros(43824,1)*NaN;
for x=1:43824
    ind = ind1(x);
    data3b(x,:) = mean(data3(ind:ind+1,:),1);
    HEM3_NEE(x,1) = mean(HEM3_NEE_measured(ind:ind+1,:),1);
    HEM3_NEE_GF(x,1) = mean(HEM3_NEE_temp(ind:ind+1,:),1);
end

Climate_HEM3_year = data3(:,1);
Climate_HEM3_day_time = data3(:,2);
Climate_HEM3_day = floor(Climate_HEM3_day_time);
Climate_HEM3_hour = repmat([sort(repmat([0:23]',2,1))],1827,1);
Climate_HEM3_hour = Climate_HEM3_hour(2:87649,1);
Climate_HEM3_year_hh = Climate_HEM3_year;
Climate_HEM3_day_hh = Climate_HEM3_day;
Climate_HEM3_hour_hh = Climate_HEM3_hour;
Climate_HEM3_minute_hh = zeros(87648,1);
Climate_HEM3_minute_hh([1:2:87648],1) = 30;
Climate_HEM3_year = Climate_HEM3_year(ind1+1,:);
Climate_HEM3_day = Climate_HEM3_day(ind1+1,:);
Climate_HEM3_hour = Climate_HEM3_hour(ind1+1,:);

Climate_HEM3_Ts_10cm = data3(:,25);


%Merge all raw data files together
year(1:365*24,1) = 2001;
year(end+1:end+365*24,1) = 2002;
year(end+1:end+365*24,1) = 2003;
year(end+1:end+365*24,1) = 2004;
day = repmat(sort(repmat([1:365]',24,1)),4,1);
hour = repmat([0:23]',365*4,1);
timestamps = [year day hour];
timestamps = timestamps(7309:30354,:);

year=[];
year(1:365*48,1) = 2001;
year(end+1:end+365*48,1) = 2002;
year(end+1:end+365*48,1) = 2003;
year(end+1:end+365*48,1) = 2004;
day = repmat(sort(repmat([1:365]',48,1)),4,1);
hour = repmat(sort([0:23 0:23]'),365*4,1);
minutes = zeros(70080,1);
minutes([2:2:70080],1) = 30;
timestamps_hh = [year day hour minutes];
timestamps_hh = timestamps_hh(14617:60706,:);

HEM_year = [Climate_HEM1_year ; timestamps(:,1) ; Climate_HEM2_year ; Climate_HEM3_year];
HEM_day = [Climate_HEM1_day ; timestamps(:,2) ; Climate_HEM2_day ; Climate_HEM3_day];
HEM_hour = [Climate_HEM1_hour ; timestamps(:,3) ; Climate_HEM2_hour ; Climate_HEM3_hour];
HEM_timestamps = [HEM_year HEM_day HEM_hour];

HEM_year_hh = [Climate_HEM1_year_hh ; timestamps_hh(:,1) ; Climate_HEM2_year_hh ; Climate_HEM3_year_hh];
HEM_day_hh = [Climate_HEM1_day_hh ; timestamps_hh(:,2) ; Climate_HEM2_day_hh ; Climate_HEM3_day_hh];
HEM_hour_hh = [Climate_HEM1_hour_hh ; timestamps_hh(:,3) ; Climate_HEM2_hour_hh ; Climate_HEM3_hour_hh];
HEM_minute_hh = [Climate_HEM1_minute_hh ; timestamps_hh(:,4) ; Climate_HEM2_minute_hh ; Climate_HEM3_minute_hh];
HEM_timestamps_hh = [HEM_year_hh HEM_day_hh HEM_hour_hh HEM_minute_hh];

HEM_NEE_GF = [HEM1_NEE_GF ; ones(23046,1)*NaN ; HEM2_NEE_GF ; HEM3_NEE_GF];
HEM_NEE = [HEM1_NEE ; ones(23046,1)*NaN ; HEM2_NEE ; HEM3_NEE];

HEM_PAR = [data1b(:,6) ; ones(23046,1)*NaN ; data2b(:,25) ; data3b(:,26)];
HEM_Ta_tc = [data1b(:,7) ; ones(23046,1)*NaN ; data2b(:,22) ; ones(43824,1)*NaN]; 
HEM_Ts_10cm_hh = [data1(:,9) ; ones(46090,1)*NaN ; data2(:,24) ; data3(:,25)];
HEM_Ts_10cm_hh(479:480,1) = NaN;
HEM_Ta_sonic = [ones(31817,1)*NaN ; data2b(:,19) ; data3b(:,21)];
HEM_Ta = [ones(31817,1)*NaN ; data2b(:,20) ; data3b(:,22)];
HEM_R = [ones(31817,1)*NaN ; data2b(:,14) ; data3b(:,16)];
HEM_GEE = [ones(31817,1)*NaN ; data2b(:,16) ; data3b(:,17)];

clear year day hour timestamps Climate_HEM1_year Climate_HEM2_year Climate_HEM3_year Climate_HEM1_day Climate_HEM2_day Climate_HEM3_day ...
    Climate_HEM1_hour Climate_HEM2_hour Climate_HEM3_hour Climate_HEM1_day_time Climate_HEM2_day_time Climate_HEM3_day_time data1 data2 data3 data1b data2b data3b ...
    EC_HEM1_NEE EC_HEM1_NEE_estimated EC_HEM1_NEE_flag EC_HEM1_NEE_measured EC_HEM1_NEE_temp ind ind1 x ...
    timestamps_hh HEM_year_hh HEM_day_hh HEM_hour_hh HEM_minute_hh Climate_HEM1_year_hh Climate_HEM2_year_hh Climate_HEM3_year_hh ...
    Climate_HEM1_day_hh Climate_HEM2_day_hh Climate_HEM3_day_hh Climate_HEM1_hour_hh Climate_HEM2_hour_hh Climate_HEM3_hour_hh ...
    Climate_HEM1_minute_hh Climate_HEM2_minute_hh Climate_HEM3_minute_hh ...
    minutes HEM_year HEM_day HEM_hour Climate_HEM1_Ts_10cm Climate_HEM2_Ts_10cm Climate_HEM3_Ts_10cm ...
    HEM1_NEE HEM1_NEE_estimated HEM1_NEE_measured HEM1_NEE_flag HEM1_NEE_temp HEM1_NEE_measured2 HEM2_NEE_measured ...
    HEM3_NEE_measured HEM2_NEE_temp HEM3_NEE_temp HEM1_NEE_GF HEM2_NEE_GF HEM3_NEE_GF HEM2_NEE HEM3_NEE



%LPH Tower
data1 = csvread([path,'\LPH Tower - eddy flux 2002-04.csv'],2,0);

%Convert to hourly data
ind1 = 2:2:45235;
data1b = zeros(22617,26)*NaN;
for x=1:22617
    ind = ind1(x);
    data1b(x,:) = mean(data1(ind:ind+1,:),1);
end

Climate_LPH1_year = data1(:,1);
Climate_LPH1_year_hh = data1(:,1);
Climate_LPH1_day_time = data1(:,2);
Climate_LPH1_day = floor(Climate_LPH1_day_time);
Climate_LPH1_day_hh = Climate_LPH1_day;
Climate_LPH1_hour = repmat([sort(repmat([0:23]',2,1))],944,1);
Climate_LPH1_hour = Climate_LPH1_hour(31:45265,1);
Climate_LPH1_hour_hh = Climate_LPH1_hour;
Climate_LPH1_minute_hh = zeros(size(Climate_LPH1_year_hh));
Climate_LPH1_minute_hh([2:2:45235],1) = 30;
Climate_LPH1_year = Climate_LPH1_year(ind1+1,:);
Climate_LPH1_day = Climate_LPH1_day(ind1+1,:);
Climate_LPH1_hour = Climate_LPH1_hour(ind1+1,:);

Climate_LPH1_Ts_10cm_hh = data1(:,22);

 
%LPH part #2
data2 = csvread([path,'\LPH Tower - eddy flux since 2005.csv'],1,0);

%Convert to hourly data
ind1 = 1:2:87648;
data2b = zeros(43824,25)*NaN;
for x=1:43824
    ind = ind1(x);
    data2b(x,:) = mean(data2(ind:ind+1,:),1);
end

Climate_LPH2_year = data2(:,1);
Climate_LPH2_year_hh = data2(:,1);
Climate_LPH2_day_time = data2(:,2);
Climate_LPH2_day = floor(Climate_LPH2_day_time);
Climate_LPH2_day_hh = floor(Climate_LPH2_day_time);
Climate_LPH2_hour = repmat([sort(repmat([0:23]',2,1))],1827,1);
Climate_LPH2_hour = Climate_LPH2_hour(2:87649,1);
Climate_LPH2_hour_hh = Climate_LPH2_hour;
Climate_LPH2_minute_hh = zeros(size(Climate_LPH2_year_hh));
Climate_LPH2_minute_hh([1:2:87648],1) = 30;
Climate_LPH2_year = Climate_LPH2_year(ind1+1,:);
Climate_LPH2_day = Climate_LPH2_day(ind1+1,:);
Climate_LPH2_hour = Climate_LPH2_hour(ind1+1,:);

Climate_LPH2_Ts_10cm_hh = data2(:,21);


%Merge all data together
LPH_year = [Climate_LPH1_year ; Climate_LPH2_year];
LPH_day = [Climate_LPH1_day ; Climate_LPH2_day];
LPH_hour = [Climate_LPH1_hour ; Climate_LPH2_hour];
LPH_timestamps = [LPH_year LPH_day LPH_hour];
LPH_year_hh = [Climate_LPH1_year_hh ; Climate_LPH2_year_hh];
LPH_day_hh = [Climate_LPH1_day_hh ; Climate_LPH2_day_hh];
LPH_hour_hh = [Climate_LPH1_hour_hh ; Climate_LPH2_hour_hh];
LPH_minute_hh = [Climate_LPH1_minute_hh ; Climate_LPH2_minute_hh];
LPH_timestamps_hh = [LPH_year_hh LPH_day_hh LPH_hour_hh LPH_minute_hh];
LPH_Ta = [data1b(:,19) ; data2b(:,18)];
LPH_Ta(61674:61702,1) = NaN;
LPH_Ts_10cm = [data1b(:,22) ; data2b(:,21)];
LPH_Ts_10cm(9949:9964,1) = NaN;
LPH_Ts_10cm_hh = [Climate_LPH1_Ts_10cm_hh ; Climate_LPH2_Ts_10cm_hh];
LPH_Ts_10cm_hh(19899:19928,1) = NaN;
LPH_Ts_10cm_hh(21264:21266,1) = NaN;
LPH_Ts_10cm_hh(33266:33303,1) = NaN;
LPH_NEE = [data1b(:,11) ; data2b(:,12)];
ind = find(LPH_NEE>=-0.05 & LPH_NEE<=0.02);
ind(ind<57734)=[];
LPH_NEE(ind,1) = NaN; %Delete bad data
LPH_NEE_GF = [data1b(:,12) ; data2b(:,13)];
LPH_NEE_GF(ind,1) = NaN; %Delete bad data

LPH_R = [data1b(:,14) ; ones(43824,1)*NaN];
LPH_GEE = [data1b(:,15) ; ones(43824,1)*NaN];

clear Climate_LPH1_year Climate_LPH2_year Climate_LPH1_day_time Climate_LPH2_day_time Climate_LPH1_day Climate_LPH2_day Climate_LPH1_hour Climate_LPH2_hour ...
    ind ind1 data1 data2 data1b data2b x Climate_LPH1_Ts_10cm_hh Climate_LPH2_Ts_10cm_hh ...
    LPH_year Climate_LPH1_year Climate_LPH2_year LPH_day Climate_LPH1_day Climate_LPH2_day LPH_hour Climate_LPH1_hour Climate_LPH2_hour ...
    LPH_year_hh Climate_LPH1_year_hh Climate_LPH2_year_hh LPH_day_hh Climate_LPH1_day_hh Climate_LPH2_day_hh LPH_hour_hh ...
    Climate_LPH1_hour_hh Climate_LPH2_hour_hh LPH_minute_hh Climate_LPH1_minute_hh Climate_LPH2_minute_hh;



%%%%%%%%%%%
%Gap-fill soil temperature data
[c,ind] = setdiff(EMS_timestamps_hh,LPH_timestamps_hh,'rows');
Ts_timestamps = [EMS_timestamps_hh(ind,:) ; LPH_timestamps_hh];
Ts = FillSmallGapsByLinInterp(LPH_Ts_10cm_hh,4);  %Fill small gaps by linear interp (<=2h)
Ts(85000:98000,1) = NaN; %Delete year ~2007 since data looks weird (too warm compared to other sites)
Ts = [ones(size(ind))*NaN ; Ts];

[c,ia,ib] = intersect(EMS_timestamps_hh,Ts_timestamps,'rows');
EMS_timestamps_hh_temp = Ts_timestamps*NaN;
EMS_timestamps_hh_temp(ib,:) = EMS_timestamps_hh(ia,:);
EMS_Ts_20cm_temp = Ts*NaN;
EMS_Ts_20cm_temp(ib,:) = EMS_Ts_20cm_hh(ia,:);
EMS_Ts_surf_temp = Ts*NaN;
EMS_Ts_surf_temp(ib,:) = EMS_Ts_surf_hh(ia,:);

[c,ia,ib] = intersect(HEM_timestamps_hh,Ts_timestamps,'rows');
HEM_timestamps_hh_temp = Ts_timestamps*NaN;
HEM_timestamps_hh_temp(ib,:) = HEM_timestamps_hh(ia,:);
HEM_Ts_10cm_temp = Ts*NaN;
HEM_Ts_10cm_temp(ib,:) = HEM_Ts_10cm_hh(ia,:);

[c,ia,ib] = intersect(Fisher_timestamps_hh,Ts_timestamps,'rows');
Fisher_timestamps_hh_temp = Ts_timestamps*NaN;
Fisher_timestamps_hh_temp(ib,:) = Fisher_timestamps_hh(ia,:);
Fisher_Ts_10cm_temp = Ts*NaN;
Fisher_Ts_10cm_temp(ib,:) = Fisher_Ts_10cm_hh(ia,:);


%Calculate regression coefficients between reference Ts series and other Ts
%series used to gap-fill
ind=find(~isnan(Ts) & ~isnan(EMS_Ts_20cm_temp));
t=Ts(ind);
t2=EMS_Ts_20cm_temp(ind);
p_EMS20 = polyfit(EMS_Ts_20cm_temp(ind),Ts(ind),2); %r2=0.9226 p=0.0188    0.9088   -0.9934  f(x) = p1*x^2 + p2*x + p3

ind=find(~isnan(Ts) & ~isnan(HEM_Ts_10cm_temp));
t=Ts(ind);
t2=HEM_Ts_10cm_temp(ind);
p_HEM = polyfit(HEM_Ts_10cm_temp(ind),Ts(ind),2); %r2=0.9623 p=-0.0059    1.1212   -0.1540

ind=find(~isnan(Ts) & ~isnan(Fisher_Ts_10cm_temp));
t=Ts(ind);
t2=Fisher_Ts_10cm_temp(ind);
p_Fisher = polyfit(Fisher_Ts_10cm_temp(ind),Ts(ind),2); %r2=0.9656 p=-0.0016    0.8065    0.2276


%Gap-fill Ts
Ts(isnan(Ts),1) = p_HEM(1).*HEM_Ts_10cm_temp(isnan(Ts),1).^2 + p_HEM(2).*HEM_Ts_10cm_temp(isnan(Ts),1) + p_HEM(3);
Ts(isnan(Ts),1) = p_Fisher(1).*Fisher_Ts_10cm_temp(isnan(Ts),1).^2 + p_Fisher(2).*Fisher_Ts_10cm_temp(isnan(Ts),1) + p_Fisher(3);
Ts(isnan(Ts),1) = p_EMS20(1).*EMS_Ts_20cm_temp(isnan(Ts),1).^2 + p_EMS20(2).*EMS_Ts_20cm_temp(isnan(Ts),1) + p_EMS20(3);


%Calculate regression coefficients between reference Ts series and
%Ts_surface series used to gap-fill
ind=find(~isnan(Ts) & ~isnan(EMS_Ts_surf_temp));
t=Ts(ind);
t2=EMS_Ts_surf_temp(ind);
p_EMSsurf = polyfit(EMS_Ts_surf_temp(ind),Ts(ind),2); %r2=0.9867 p=0.0119    0.9316   -1.2343

%Gap-fill Ts
Ts(isnan(Ts),1) = p_EMSsurf(1).*EMS_Ts_surf_temp(isnan(Ts),1).^2 + p_EMSsurf(2).*EMS_Ts_surf_temp(isnan(Ts),1) + p_EMSsurf(3); 

Ts = FillSmallGapsByLinInterp(Ts,2000);  %Fill remaining gaps


%Add month to Ts_timestamps
month = ones(length(Ts_timestamps),1);
for x=1:length(Ts_timestamps)
    month(x,1) = fn_jday_to_month(Ts_timestamps(x,1),Ts_timestamps(x,2));
end
Ts_timestamps = [Ts_timestamps(:,1) month Ts_timestamps(:,2:4)];


%Keep data for supplementary figure (figure to show regressions between
%LPH-Ts and other Ts used to gap-filled and create Ts_ref)
temp=[Ts_timestamps Ts Fisher_Ts_10cm_temp HEM_Ts_10cm_temp EMS_Ts_20cm_temp EMS_Ts_surf_temp];   
HF_Archive_Ts_hourly = [temp(:,[1 3:10])]; %Use in Harvard Forest Archive (file Ts_ref.csv)


%%%%%%%%%%%%%%%%%%%%%%%
%Create hourly Ts vector and timestamps
ind = 1:2:length(Ts_timestamps);
Ts_gfNEE = Ts(ind,1);
timestamps_gfNEE = Ts_timestamps(ind,1:4);


%%%%%%%%%%%%%%%%%%%%%%%
%Gap-fill air temperature and PAR data
Ta_gfNEE = [EMS_filled_Ta_27m_GF ; ones([length(timestamps_gfNEE)-length(EMS_filled_Ta_27m_GF)],1)*NaN];
PAR_gfNEE = [EMS_filled_PAR_28m_GF ; ones([length(timestamps_gfNEE)-length(EMS_filled_PAR_28m_GF)],1)*NaN];
PAR_gfNEE(153085:153131) = NaN; %Delete bad data

[c,ia,ib] = intersect(HEM_timestamps,timestamps_gfNEE(:,[1 3 4]),'rows');
HEM_Ta_temp = Ta_gfNEE*NaN;
HEM_Ta_temp(ib,:) = HEM_Ta(ia,:);

%Keep HEM data to gap-fill using Fluxnet-Canada algorithm
HEM_PAR_temp = Ta_gfNEE*NaN;
HEM_PAR_temp(ib,:) = HEM_PAR(ia,:);
HEM_NEE_gfNEE = Ta_gfNEE*NaN;
HEM_NEE_gfNEE(ib,:) = HEM_NEE(ia,:);
HEM_NEE_gfNEE_GF = Ta_gfNEE*NaN;
HEM_NEE_gfNEE_GF(ib,:) = HEM_NEE_GF(ia,:);
HEM_GEE_gfNEE = Ta_gfNEE*NaN;
HEM_GEE_gfNEE(ib,:) = HEM_GEE(ia,:);
HEM_R_gfNEE = Ta_gfNEE*NaN;
HEM_R_gfNEE(ib,:) = HEM_R(ia,:);


[c,ia,ib] = intersect(Fisher_timestamps_h,timestamps_gfNEE(:,[1 3 4]),'rows');
Fisher_Ta_temp = Ta_gfNEE*NaN;
Fisher_Ta_temp(ib,:) = Fisher_Ta(ia,:);
Fisher_PAR_temp = Ta_gfNEE*NaN;
Fisher_PAR_temp(ib,:) = Fisher_PAR(ia,:);


%Gap-fill Ta
ind=find(~isnan(Ta_gfNEE) & ~isnan(HEM_Ta_temp));
p_HEM_Ta = polyfit(HEM_Ta_temp(ind),Ta_gfNEE(ind),2); %r2=0.9876

ind=find(~isnan(Ta_gfNEE) & ~isnan(Fisher_Ta_temp));
p_Fisher_Ta = polyfit(Fisher_Ta_temp(ind),Ta_gfNEE(ind),2); %r2=0.9859

Ta_gfNEE(isnan(Ta_gfNEE),1) = p_HEM_Ta(1).*HEM_Ta_temp(isnan(Ta_gfNEE),1).^2 + p_HEM_Ta(2).*HEM_Ta_temp(isnan(Ta_gfNEE),1) + p_HEM_Ta(3);
Ta_gfNEE(isnan(Ta_gfNEE),1) = p_Fisher_Ta(1).*Fisher_Ta_temp(isnan(Ta_gfNEE),1).^2 + p_Fisher_Ta(2).*Fisher_Ta_temp(isnan(Ta_gfNEE),1) + p_Fisher_Ta(3);

Ta_gfNEE = FillSmallGapsByLinInterp(Ta_gfNEE,4);  %Fill remaining gap

%Gap-fill PAR
ind=find(~isnan(PAR_gfNEE) & ~isnan(Fisher_PAR_temp));
p_Fisher_PAR = polyfit(Fisher_PAR_temp(ind),PAR_gfNEE(ind),1); %r2=0.8577

PAR_gfNEE(isnan(PAR_gfNEE),1) = p_Fisher_PAR(1).*Fisher_PAR_temp(isnan(PAR_gfNEE),1) + p_Fisher_PAR(2);

PAR_gfNEE = FillSmallGapsByLinInterp(PAR_gfNEE,4);  %Fill remaining gap


clear c ind ia ib t t2 p_Fisher p_HEM p_EMS20 p_EMSsurf EMS_Ts_20cm_hh EMS_Ts_surf_hh HEM_Ts_10cm_hh ...
    LPH_Ts_10cm_hh Fisher_Ts_10cm_hh EMS_Ts_20cm_temp EMS_Ts_surf_temp HEM_Ts_10cm_temp LPH_Ts_10cm_temp Fisher_Ts_10cm_temp ...
    EMS_timestamps_hh_temp HEM_timestamps_hh_temp LPH_timestamps_hh_temp Fisher_timestamps_hh_temp ...
    EMS_timestamps_hh HEM_timestamps_hh LPH_timestamps_hh Fisher_timestamps_hh month p_Fisher_PAR p_Fisher_Ta p_HEM_Ta x ...
    EMS_GEE_GF EMS_NEE_GF EMS_R_GF EMS_NEE HEM_GEE HEM_R HEM_NEE HEM_NEE_GF LPH_GEE LPH_R LPH_NEE LPH_NEE_GF





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Soil respiration
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load main soil respiration file
data = csvread([path,'\HF_soil_CO2_v8_with_WarmN_V2.csv'],1,0);
% HF_soil_CO2_v8_with_WarmN_V2.csv metadata
% 
% PI
% 1=Davidson, 2=Ellison, 3=Hadley, 4=Melillo, 5=Munger, 6=Varner, 7=Frey, 8=Nadelhoffer
% 
% Site	
% 1=Barre Woods, 2=Charlton, 3=Dry Down, 4=EMS, 5=EMS-autochambers, 6=Farm, 7=Hardwood, 8=Hemlock, 9=LPH, 10=Montauk,
% 11=NWF, 12=NWM, 13=NWN, 14=Prospect Hill, 15=Simes, 16=SWF, 17=SWM, 18=SWN, 19=Black gum swamp
% 
% Study number (as listed in the paper) = PI-Site
% S1=1-4
% S2=1-3
% S3=1-15
% S4=1-2
% S5=1-6
% S6=1-7
% S7=1-10
% S8=1-11
% S9=1-12
% S10=1-13
% S11=1-16
% S12=1-17
% S13=1-18
% S14=2-15
% S15=7-14
% S16=7-XX - Not included in main file - must be appended later in the program
% S17=3-8
% S18=3-9
% S19=4-1
% S20=4-14
% S21=5-4
% S22=8-XX - Not included in main file - must be appended later in the program
% S23=6-5
% 
% 
% Site.desc (Vegetation type)
% 1=Deciduous, 2=Hemlock, 3=Mixed deciduous-coniferous, 4=Red pine plantation, 5=Sphagnum bog, 6=swamp, 7=Deciduous wet forest, 8=Deciduous-logged, 9=Mixed deciduous-coniferous-logged
% 
% Soil series	
% 1=Agawam fine sandy loam, 2=Becket-Skerry association, 3=Brookfield loam, 4=Bucksport and Wonsqueak mucks, 5=Canton, 6=Canton-Scituate association, 7=Charlton,
% 8=Montauk, 9=Peat, 10=Peru-Marlow association, 11=Scituate, 12=Whitman, 13=Becket-Monadnock association, 14=Hinckley loamy sand, 15=Lyman-Tunbridge-Berkshire association,
% 16=Pillsbury-Peacham association, 17=Tunbridge-Lyman-Berkshire association, 18=Gloucester
% 
% Soil.drainage.class	
% 1=ED, 2=WD, 3=MWD, 4=M, 5=P, 6=PD, 7=VP, 8=VPD, 9=SED
% 
% obs.exp	
% 1= exp, 2=obs
% 
% Auto.manual	
% 1=autochamber, 2=manual
% 
% treat.level	
% 0=control (NaN), 1=heated, 2=control, 3=disturbance control, 4=drydown, 5=girdled, 6=hardwood, 7=logged, 8=logged, 9=fertilized (5g N m-2 y-1), 10=fert+heated
% 11=double litter inputs, 12=no litter inputs, 13=no root inputs (trenched), 14= no root and litter inputs,
% 15=O and A horizon replaced with B-horizon soil, 16=fertilized (15g N m-2 y-1), 17=fertilized (5g N m-2 y-1) + sulphur
%%%%%%%%%%%%%%%%%%%%%%%%%%%


ind = find(data(:,1)==1991 & data(:,4)<301);
data(ind,:) = []; %Delete data for which we have no soil temperature series
ind = find(data(:,1)==2006 & data(:,2)==3 & data(:,6)==3 & data(:,7)==8);
data(ind,:) = []; %Delete data for which soil temperature measurements seem bad (too high)
ind = find(data(:,1)==2002 & data(:,6)==3 & data(:,7)==9);
data(ind,:) = []; %Delete bad soil respiration measurements. There does not seem to be a link between Rs and Ts.
indi = find(data(:,1)==2005 & data(:,4)==320 & data(:,5)==1033 & data(:,6)==6 & data(:,7)==5);
indj = find(data(:,1)==2005 & data(:,4)==353 & data(:,5)==1243 & data(:,6)==6 & data(:,7)==5);
data(indi:indj,:) = []; %Delete bad soil respiration measurements. Way too low compared to other winters. Snow on the ground? Creates problems with bias correction - annual sum becomes way too high.
ind = find(data(:,1)==2006 & data(:,6)==5 & data(:,7)==4);
data(ind,:) = []; %Delete data; not enough measurements to get good regression
ind = find(isnan(data(:,21)));
data(ind,21) = 2; %Replace NaNs by Manual measurements


%Recalculate day and month from jday data (there were errors - let's consider julian days are ok)
for x=1:length(data)
    [d,m] = fn_jday_to_day(data(x,1),data(x,4));
    data(x,2) = m;
    data(x,3) = d;
end

%Compute hour and minutes, and correct other timestamps accordingly
hour =floor(data(:,5)./100);
minute = data(:,5)-(100.*hour);
minute(minute<=30) = 30;
ind = find(minute>30);
minute(ind,1) = 0;
hour(ind,1) = hour(ind,1)+1;
ind2 = find(hour>=24);
hour(ind2,1) = 0;
data(ind2,3) = data(ind2,3)+1;
data(ind2,4) = data(ind2,4)+1;
for x=1991:2006
    for xx=1:12
        nb_days = fn_nb_days_in_month(x,xx);
        ind3 = find(data(:,1)==x & data(:,2)==xx & data(:,3)>nb_days);
        data(ind3,2) = data(ind3,2)+1;
        data(ind3,3) = 1;
    end
end
ind4 = find(data(:,2)>12);
data(ind4,2) = 1;
data(ind4,1) = data(ind4,1)+1;
        
data = [data(:,1:4) hour minute data(:,6:end)];
ind = find(isnan(data(:,20)));
data(ind,20) = 0; %Replace NaNs by 0 for treatment number

%Create treatment numbers for Munger experiment
ind=find(ismember(data(:,13),[101 102 105 113 114 115 117 118]) & data(:,1)<=2000);
data(ind,:) = []; %Delete data from logged collars before logging was done
ind=find(ismember(data(:,13),[101 102 105 113 114 115 117 118]) & data(:,1)>=2002);
data(ind,20) = 8;
ind=find(ismember(data(:,13),[103 104 106 107 108 109 110 111 112 116]));
data(ind,20) = 2;



%DIRT experiment (was not included in main soil respiration file)
[num_Rs,txt,raw] = xlsread([path,'\DIRT_soil_respiration.xls'],'hf007-06-soil-respiration'); %Metadata is included in the file
[num_Ts,txt,raw] = xlsread([path,'\DIRT_soil_temperature.xls'],'hf007-05-soil-temperature'); %Metadata is included in the file

DIRT_Rs = [num_Rs(:,1:2) num_Rs(:,3:end).*1000000./(3600*1000*12)];
DIRT_Ts = num_Ts;

%Compute daily reference Ts
temp = unique(Ts_timestamps(:,[1 3]),'rows');
ind = find(temp(:,1)>2001 | temp(:,1)<1992);
temp(ind,:)=[];
Ts_ref_daily=ones(length(temp),3)*NaN;
for x=1:length(temp)
    y=temp(x,1);
    d=temp(x,2);
    ind = find(Ts_timestamps(:,1)==y & Ts_timestamps(:,3)==d);
    Ts_ref_daily(x,:) = [temp(x,:) mean(Ts(ind,1))];
end

%Scale Ts_ref_daily to DIRT Ts for each treatment and collar
[c,ia,ib]=intersect(DIRT_Ts(:,1:2),Ts_ref_daily(:,1:2),'rows');
x=DIRT_Ts(ia,3:end);
y=Ts_ref_daily(ib,3);
DIRT_Ts_gf=Ts_ref_daily(:,1:2);
% for xx=1:6
for xx=1:21
    x1=x(:,xx); %We used Ts at 5cm depth - 15cm is also available, and sometimes 10 or 12cm are available too
    ind=find(~isnan(x1) & ~isnan(y));
    p = polyfit(x(ind,xx),y(ind),1);
    DIRT_Ts_gf(:,xx+2) = Ts_ref_daily(:,3).*p(1) + p(2);
end


%ChronicN experiment (was not included in main soil respiration file)
[num,txt,raw] = xlsread([path,'\ChronicN_Rs_data_2009.xls'],'Calculated Flux');
ind = find(isnan(num(:,24)));
num(ind,:) = []; %Delete when timestamps are missing
temp = num(:,24);
temp1 = floor(temp./100);
temp2 = temp-100.*temp1;
ind = find(temp2>0 & temp2<30);
temp2(ind) = 30;
ind = find(temp2>30);
temp2(ind) = 0;
temp1(ind) = temp1(ind)+1;

ChronicN_data = [ones(length(num),1)*2009 num(:,23) temp1 temp2 num(:,[5 36 8])]; %Jday, hour, minute, collar, Rs, Ts
ind = find(isnan(ChronicN_data(:,6)) | isnan(ChronicN_data(:,7)));
ChronicN_data(ind,:)=[]; %Delete when Rs or Ts is missing

%Scale reference Ts to ChronicN Ts for each treatment
Ts_ref_temp = [];
for x=1:length(ChronicN_data)
    [c,ia,ib]=intersect(ChronicN_data(x,1:4),Ts_timestamps(:,[1 3:5]),'rows');
    Ts_ref_temp(x,1) = Ts(ib);
end


%Compute bootstrap to estimate uncertainty (except for DIRT experiment
%which is done later)
fn_compute_Rs_ChronicN_bias_corr_bootstrap(ChronicN_data,Ts_ref_temp,Ts,Ts_timestamps,path);

%Compute Rs vs Ts regression coefficients per treatment & vegetation cover type
ChronicN_Ts_reg=[];
ChronicN_Ts_pred=[];
ChronicN_Ts_pred2=[];
for x=1:8
    temp = x*100;
    ind = find(ChronicN_data(:,5)>=temp & ChronicN_data(:,5)<=temp+5);
    
    p = polyfit(ChronicN_data(ind,7),Ts_ref_temp(ind),2);
    ChronicN_Ts_pred(:,x) = p(1).*Ts_ref_temp.^2 + p(2).*Ts_ref_temp + p(3);
    ChronicN_Ts_pred2(:,x) = p(1).*Ts.^2 + p(2).*Ts + p(3);
     
    F = [Ts_ref_temp(ind).^0 Ts_ref_temp(ind) Ts_ref_temp(ind).^2];     % make design matrix [1,x]
    c = F\ChronicN_data(ind,7);                  % get least-squares fit
    res = ChronicN_data(ind,7) - F*c;            % calculate residuals
    r2 = 1 - var(res)/var(ChronicN_data(ind,7)); % calculate R^2
    
    n = length(ind);
    ChronicN_Ts_reg = [ChronicN_Ts_reg ; [p r2 n]];
end
%Reorder treatments (Hardwood control/highN/lowN/lowN+S - Pine control/highN/lowN/lowN+S)
ChronicN_Ts_reg = ChronicN_Ts_reg([1 3 2 4 5 7 6 8],:);

%Estimate Rs
Q10_final_ChronicN=[];
var_pred_ChronicN=[];
Rs_hh=[];
for x=1:8
    temp = x*100;
    ind = find(ChronicN_data(:,5)>=temp & ChronicN_data(:,5)<=temp+5);
    
    temp = ChronicN_data(ind,6);
    ind0 = find(temp<=0);
    temp(ind0) = 1;
    Rs_flux2 = log(temp);
    Rs_flux2(ind0) = 0;
    
    [p,Q10,R10,r2]=get_Q10_R10_r2(ChronicN_Ts_pred(ind,x),Rs_flux2,'lin'); %All collars together
    Q10_final_ChronicN = [Q10_final_ChronicN ; [2009 x p Q10 R10 r2]];
    
    %Estimate Rs for each day, including bias correction
    y_pred = p(1).*ChronicN_Ts_pred(ind,x) + p(2);
    var_pred = sum(((Rs_flux2 - y_pred).^2)./(n-2));
    Rs_hh(:,x) = (R10.*Q10.^((ChronicN_Ts_pred2(:,x)-10)./10)).*exp(0.5.*var_pred);
    var_pred_ChronicN = [var_pred_ChronicN ; [var_pred]];
end

%Keep only year 2009
ind = find(Ts_timestamps(:,1)==2009);
Rs_hh = Rs_hh(ind,:);
temp = Ts_timestamps(ind,:);

%Compute daily and then monthly mean (gC m-2 d-1)
Rs_hh2 = Rs_hh.*(3600*24*12/1000000); %Convert from umolCO2 m-2 s-1 to gC m-2 d-1
y=0;
ChronicN_Rs_monthly=[];
for x=1:12
    [i,j] = fn_month_to_jday(2009,x);
    ind = find(temp(:,3)>=i & temp(:,3)<=j);
    y=y+1;
    ChronicN_Rs_monthly(y,:) = [2009 x mean(Rs_hh2(ind,1:8),1)];
end
ChronicN_Rs_monthly2 = ChronicN_Rs_monthly(:,1:2);
for x=1:12
    i = fn_nb_days_in_month(2009,x);
    ChronicN_Rs_monthly2(x,3:10) = [ChronicN_Rs_monthly(x,3:10).*i];
end


%%%%%%%
%Estimate Rs for DIRT experiment
option=0;
[Rs_DIRT_monthly,Q10_final_DIRT,Q10_final_yearly_DIRT,Q10_final_all_years_DIRT,Q10_final_DIRT2,var_pred_for_Rs_MAT_calculation_DIRT,Q10_final_yearly_DIRT2] = fn_compute_Rs_DIRT_bias_corr(DIRT_Rs,DIRT_Ts_gf,option);
var_pred_for_Rs_MAT_calculation_DIRT = var_pred_for_Rs_MAT_calculation_DIRT(7:12,:); %Keep only coefficients for control plots

%Compute bootstrap to estimate uncertainty (DIRT experiment only; the
%others were done earlier)
temp = fn_compute_Rs_DIRT_bias_corr_bootstrap(DIRT_Rs,DIRT_Ts_gf,option);


%%%%%%%%%%%%%%%%%%%%%%%%
Rs_experiments = unique([data(:,7) data(:,8) data(:,22)],'rows');
%Rs_experiments metadata

% Column 1 = PI
% 1=Davidson, 2=Ellison, 3=Hadley, 4=Melillo, 5=Munger, 6=Varner, 7=Frey, 8=Nadelhoffer
% 
% Column 2 = Site	
% 1=Barre Woods, 2=Charlton, 3=Dry Down, 4=EMS, 5=EMS-autochambers, 6=Farm, 7=Hardwood, 8=Hemlock, 9=LPH, 10=Montauk,
% 11=NWF, 12=NWM, 13=NWN, 14=Prospect Hill, 15=Simes, 16=SWF, 17=SWM, 18=SWN, 19=Black gum swamp
% 
% Study number (as listed in the paper) = PI-Site
% S1=1-4
% S2=1-3
% S3=1-15
% S4=1-2
% S5=1-6
% S6=1-7
% S7=1-10
% S8=1-11
% S9=1-12
% S10=1-13
% S11=1-16
% S12=1-17
% S13=1-18
% S14=2-15
% S15=7-14
% S16=7-XX - Not included in main file - must be appended later in the program
% S17=3-8
% S18=3-9
% S19=4-1
% S20=4-14
% S21=5-4
% S22=8-XX - Not included in main file - must be appended later in the program
% S23=6-5

%Column 3 = Auto.manual	
% 1=autochamber, 2=manual
%%%%%%%%%%%%%%%%%%%%%%%%

%For each experiment, get all data needed for calculations
option = 0; %0=no figures, 1=figures
Rs_data= [sort(repmat([1991:2009]',12,1)) repmat([1:12]',2009-1991+1,1)];
Rs_data = [Rs_data ones(length(Rs_data),75)*NaN];
col = 2;
Rs_data_winter=[];
for exp_ind=1:length(Rs_experiments)
    [Rs_data_processed,temp_Rs_data_winter] = fn_compute_Rs_Harvard_Forest_bias_corr(exp_ind,data,Ts,Ts_timestamps,Rs_experiments,option,path);
    Rs_data_winter = [Rs_data_winter ; temp_Rs_data_winter];
    if ismember(Rs_experiments(exp_ind,1:2),[1 2 ; 1 4 ; 1 6 ; 1 7 ; 1 10 ; 1 11 ; 1 12 ; 1 13 ; 1 16 ; 1 17 ; 1 18 ; 3 8 ; 3 9 ; 6 5],'rows') %1 treatment
        [c,ia,ib] = intersect(Rs_data(:,1:2),Rs_data_processed(:,1:2),'rows');
        Rs_data(ia,col+1:col+3) = Rs_data_processed(ib,[7:9]);
        col = col+3;
    elseif ismember(Rs_experiments(exp_ind,1:2),[1 3 ; 4 1 ; 1 15 ; 4 14 ; 5 4 ; 7 14],'rows') %2, 3 or 4 treatments
        treatments = unique(Rs_data_processed(:,5));
        for x=1:length(treatments)
            ind = find(Rs_data_processed(:,5)==treatments(x));
            temp = Rs_data_processed(ind,[1 2 7:9]);
            [c,ia,ib] = intersect(Rs_data(:,1:2),temp(:,1:2),'rows');
            Rs_data(ia,col+1:col+3) = temp(ib,[3:5]);
            col = col+3;
        end
    end
end
fn_delete_Excel_default_sheets([path,'\'],'Harvard_Forest_soil_respiration_final.xls');
fn_delete_Excel_default_sheets([path,'\'],'Harvard_Forest_soil_respiration_Q10_final.xls');


%For each experiment (control only), get all data needed for calculations
option = 0; %0=no figures, 1=figures
Rs_data_controls= [sort(repmat([1991:2009]',12,1)) repmat([1:12]',2009-1991+1,1)];
Rs_data_controls = [Rs_data_controls ones(length(Rs_data_controls),57)*NaN];
Rs_data_controls_daily=[];
All_Ts=[];
col = 2;
for exp_ind=1:length(Rs_experiments)
    [Rs_controls_data_processed,Rs_controls_daily_data_processed,Ts_temp] = fn_compute_Rs_controls_Harvard_Forest_bias_corr(exp_ind,data,Ts,Ts_timestamps,Rs_experiments,option,path);
    if ~ismember(Rs_experiments(exp_ind,1:2),[2 15 ; 3 19],'rows') %Do nothing for experiments with no timestamps
        [c,ia,ib] = intersect(Rs_data_controls(:,1:2),Rs_controls_data_processed(:,1:2),'rows');
        Rs_data_controls(ia,col+1:col+3) = Rs_controls_data_processed(ib,[7:9]);
        col = col+3;
        
        All_Ts = [All_Ts ; Ts_temp];
        
        if exp_ind==1
            Rs_data_controls_daily = [Rs_controls_daily_data_processed ones(length(Rs_controls_daily_data_processed),18)*NaN];
            y=3;
        else
            y=y+1;
            Rs_data_controls_daily(:,y) = Rs_controls_daily_data_processed(:,3);
        end
    end
end
fn_delete_Excel_default_sheets([path,'\'],'Harvard_Forest_soil_respiration_controls_final.xls');
fn_delete_Excel_default_sheets([path,'\'],'Harvard_Forest_soil_respiration_controls_Q10_final.xls');

%Delete missing data and separate between autochambers and manual
ind = find(isnan(All_Ts(:,2)));
All_Ts(ind,:) = [];
ind = find(All_Ts(:,1)==3 | All_Ts(:,1)==21);
All_Ts_auto = All_Ts(ind,:);
All_Ts_manual = All_Ts;
All_Ts_manual(ind,:)=[];

%Compute number of Ts measurements (manual or autochambers) per 1degC bin
[All_Ts_auto_bins,All_Ts_auto_bins_limits,All_Ts_auto_nb_data] = fn_bin_data(1,All_Ts_auto(:,2),ones(length(All_Ts_auto),1));
temp = [All_Ts_auto_bins_limits All_Ts_auto_nb_data]; %Use for Harvard Forest Archive (file Ts_number_of_measurements.csv, column 3)
[All_Ts_manual_bins,All_Ts_manual_bins_limits,All_Ts_manual_nb_data] = fn_bin_data(1,All_Ts_manual(:,2),ones(length(All_Ts_manual),1));
temp2 = [All_Ts_manual_bins_limits All_Ts_manual_nb_data]; %Use for Harvard Forest Archive (file Ts_number_of_measurements.csv, column 4)


%%%%%%%%
%Merge all Nmin and Rs data for final Excel file (slightly modified
%copy/paste of the lines just above)
Nmin_Rs_data = [sort(repmat([1991:2009]',12,1)) repmat([1:12]',2009-1991+1,1)];
Nmin_Rs_data = [Nmin_Rs_data zeros([length(Nmin_Rs_data)],126)*NaN];
[c,ia,ib] = intersect(Nmin_Rs_data(:,1:2),Rs_data(:,1:2),'rows');
Nmin_Rs_data(ia,3:92) = Rs_data(ib,3:end);
Nmin_Rs_data([1:6 end-2:end],:)=[];

clear Climate_* LPH_* Fisher_* Shaler_* Precip_* EMS_filled_PAR_28m_GF EMS_filled_Ta_27m_GF EMS_timestamps c ib ia num_Rs num_Ts raw txt hour minute HEM_NEE_gfNEE_GF HEM_PAR HEM_T* HEM_timestamps;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute treatment/control ratios for all experiments with treatments
%(REFERENCE IS THE EXPERIMENTS' CONTROLS, NOT HF CONTROLS)
[c,ia,ib] = intersect(Rs_DIRT_monthly(:,1:2),Nmin_Rs_data(:,1:2),'rows');
temp = ones(size(Nmin_Rs_data,1),size(Rs_DIRT_monthly,2)-3)*NaN;
temp(ib,:) = Rs_DIRT_monthly(ia,4:end);
Nmin_Rs_data2 = [Nmin_Rs_data(:,1:end-1) temp];
Nmin_Rs_data2 = [Nmin_Rs_data2 ones(size(Nmin_Rs_data2,1),8)*NaN];
[c,ia,ib] = intersect(Nmin_Rs_data2(:,1:2),ChronicN_Rs_monthly(:,1:2),'rows');
Nmin_Rs_data2(ia,end-7:end) = ChronicN_Rs_monthly(ib,[3 5 4 6 7 9 8 10]); %Reorder treatments (Hardwood: control, high N, low N, low N + sulphur - Pine: control, high N, low N, low N + sulphur)

%Convert data to monthly instead of daily (multiply by nb days in month)
controls = Nmin_Rs_data2(:,[6 33 33 60 66 66 72 81 81 81 132 132 132 132 132 151 151 151 155 155 155]); %drydown / Davidson-Simes / Barre Woods / ProspectHill / Munger / Warm+N / DIRT / ChronicN
treatments = Nmin_Rs_data2(:,[9 36 39 57 63 69 75 84 87 90 128 136 140 144 148 152 153 154 156 157 158]);
controls(controls==0)=NaN;
treatments(treatments==0)=NaN;

%Compute monthly sums
for x=1:size(controls,1)
    nb_days = fn_nb_days_in_month(Nmin_Rs_data(x,1),Nmin_Rs_data(x,2));
    controls(x,:) = controls(x,:).*nb_days;
    treatments(x,:) = treatments(x,:).*nb_days;
end
controls(controls==0)=NaN;
treatments(treatments==0)=NaN;
 
%Compute ratios
y=0;
ratio_summer=[];
for x=1991:2009
    y=y+1;
    ind2 = find(Nmin_Rs_data(:,1)==x & Nmin_Rs_data(:,2)>=4 & Nmin_Rs_data(:,2)<=10);
    ratio_summer(y,:) = [x nansum(treatments(ind2,:))./nansum(controls(ind2,:))];
end
ratio_summer(isinf(ratio_summer))=NaN;
ratio_summer(ratio_summer==0)=NaN;

%%%%%
%Format ratio data
ratio_summer = [ratio_summer(:,1:4) ratio_summer(:,8) ratio_summer(:,[5:7 9:11 17:22 12:16])]; %Change study order to drydown / Davidson-Simes / Munger / Barre Woods / ProspectHill / Warm+N / ChronicN / DIRT
ratio_summer(1,:)=[];
%Separate years before/after treatments in different columns
ratio_summer = ratio_summer(:,[1:2 2 3 3 4 4 5 6 6 7:22]);
ratio_summer(12:13,2)=NaN;
ratio_summer(10:11,3)=NaN;
ratio_summer(14:18,4)=NaN;
ratio_summer(12:13,5)=NaN;
ratio_summer(14:18,6)=NaN;
ratio_summer(13,7)=NaN;
ratio_summer(12:18,9)=NaN;
ratio_summer(11,10)=NaN;
ratio_summer(:,[18 21]) = [];
ratio_summer = [ratio_summer(:,1:13) ratio_summer(:,15) ratio_summer(:,14) ratio_summer(:,16:end)];
ratio_summer2 = ratio_summer; %Use for Harvard Forest Archive (file Rs_treatment_Rs_control_ratio.csv)

%Compute CV of mean treatment/control ratios
ratio_summer(:,[3 4 6 9])=[];  %Use for Harvard Forest Archive (CV of treatments, so that's why we deleted data before/after treatments) (file CV_treatments_variability.csv)
temp = nanmean(ratio_summer(:,2:end),1); %Use for Harvard Forest Archive (last row of file CV_treatments_variability.csv)


%%%%%%%%%%%%%
%Estimate uncertainty on ratios
Rs_experiments2 = Rs_experiments([1:13 15:16 18:end],:); %Same metadata as Rs_experiments
for exp_ind=1:length(Rs_experiments2)
    temp = fn_compute_Rs_Harvard_Forest_bias_corr_bootstrap(exp_ind,data,Ts,Ts_timestamps,Rs_experiments2,option,path);
end

%Load data and compute growing season totals
[data_out] = fn_load_ratio_bootstrap_data(path); %To copy/save in file Ratio_bootstrap_data.csv

%%%%%%%%%%%%%%
%Compute ratios and prepare matrix to use
data_ratio_bootstrap = csvread([path,'\Ratio_bootstrap_data.csv'],1,0);
% Ratio_bootstrap_data.csv metadata
% 
% Study = study number as in the paper. Study #160 = same as study #16, but #16 is for hardwoods forests and #160 for red pine plantations
% 
% Treatment
% 0=control 1=heated, 3=disturbance control, 4=drydown, 5=girdled, 7=logged, 8=logged, 9=fertilized (5g N m-2 y-1), 10=fert+heated
% 11=double litter inputs, 12=no litter inputs, 13=no root inputs (trenched), 14= no root and litter inputs,
% 15=O and A horizon replaced with B-horizon soil, 16=fertilized (15g N m-2 y-1), 17=fertilized (5g N m-2 y-1) + sulphur
%%%%%%%%%%%%%%
studies_temp = unique(data_ratio_bootstrap(:,1));
temp2=[];
for x=1:length(studies_temp)
    ind = find(data_ratio_bootstrap(:,1)==studies_temp(x));
    temp_data = data_ratio_bootstrap(ind,:);
    year_temp = unique(temp_data(:,3));
    for xx=1:length(year_temp)
        ind = find(temp_data(:,3)==year_temp(xx));
        temp_data2 = temp_data(ind,:);
        ind0 = find(temp_data2(:,2)==0);
        treat_temp = unique(temp_data2(:,2));
        treat_temp(treat_temp==0)=[];
        for xxx=1:length(treat_temp)
            ind = find(temp_data2(:,2)==treat_temp(xxx));
            temp = temp_data2(ind,4:end)./temp_data2(ind0,4:end);
            temp2 = [temp2 ; [studies_temp(x) treat_temp(xxx) year_temp(xx) temp]];
        end
    end
end
temp2 = sortrows(temp2);
data_ratio_bootstrap = temp2(:,4:end)';
%Change study order to drydown / Davidson-Simes / Munger / Barre Woods / ProspectHill / Warm+N / ChronicN / DIRT
%remove sulfur treatments, and insert columns of NaNs
data_ratio_bootstrap = [data_ratio_bootstrap(:,1:4) ones(200,5)*NaN data_ratio_bootstrap(:,5:17) ones(200,5)*NaN data_ratio_bootstrap(:,63:64) ones(200,5)*NaN ...
    data_ratio_bootstrap(:,33:40) ones(200,5)*NaN data_ratio_bootstrap(:,41:62) ones(200,5)*NaN data_ratio_bootstrap(:,[18:21 26:29 22:25]) ones(200,5)*NaN ...
    data_ratio_bootstrap(:,[31 30 96 95]) ones(200,5)*NaN data_ratio_bootstrap(:,65:94)];
temp = data_ratio_bootstrap;
temp(:,[5:9 23:27 30:34 43:47 70:74 87:91 96:100])=[]; %Use for Harvard Forest Archive (file Ratio_bootstrap_data.csv)
ratio_summer2 = [repmat(ratio_summer2(:,1),size(ratio_summer2,2)-1,1) reshape(ratio_summer2(:,2:end),[],1)];
ind = find(isnan(ratio_summer2(:,2)));
ratio_summer2(ind,:)=[];
ratio_summer2 = [ratio_summer2(:,1)]'; %Use for Harvard Forest Archive (row #2 of file Ratio_bootstrap_data.csv)



%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%GET R CODE FROM ADRIEN
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%Compute R10 and Q10 per collar (mainly using R code, but format resulting
%data here)
rawdata = csvread([path,'\q10r10_Adrien_rawdata_102612.csv'],1,0); %Data computed using R code


veg_types = unique(rawdata(:,end));
min_num = 0;
max_num = max(rawdata(:,3));
nb_data = [0:round(max_num.*3)]';
nb_data = [nb_data ones(length(nb_data),length(veg_types))*NaN];
for x=1:5
    ind = find(rawdata(:,4)==x);
    x_data = rawdata(ind,3);
    y_data = x_data;
    [data_out,bins_limits,nb_data_temp] = fn_bin_data([1/3],x_data,y_data);
    bins_limits = round(bins_limits.*3);
    for xx=1:length(nb_data_temp)
        ind = find(nb_data(:,1)==bins_limits(xx,1));
        nb_data(ind,x+1) = nb_data_temp(xx,1);
    end
end
nb_data(:,1) = nb_data(:,1)./3;
nb_data(isnan(nb_data))=0;

temp = sum(nb_data(:,2:end-1),1);
interval = max(temp);
interval2 = 100;
nb_lines = 5*interval + 6*interval2;
temp = ones(nb_lines,5)*NaN;
temp2 = ones(nb_lines,5)*NaN;
y=interval2;
for x=1:5
    ind = find(rawdata(:,4)==x);
    ind2 = fn_space_data(ind,interval);
    temp(y+ind2,x) = rawdata(ind,3);
    temp2(y+ind2,x) = rawdata(ind,2);
    y = y+interval+interval2;
end
temp3 = temp*NaN;
temp4 = temp2*NaN;
for x=1:5
    ind = find(~isnan(temp(:,x)));
    temp3(1:length(ind),x) = temp(ind,x);
    ind = find(~isnan(temp2(:,x)));
    temp4(1:length(ind),x) = temp2(ind,x);
end
ind = length(find(~isnan(temp3(:,1))));
temp3 = temp3(1:ind,:); %Use for Harvard Forest Archive (file R10.csv)
temp4 = temp4(1:ind,:); %Use for Harvard Forest Archive (file Q10.csv)



%%%%%%%%%%%%%%%%%%%%%%%%
% Tree basal area
[num1,txt1,raw1] = xlsread([path,'\tow_06_dbh.xls']); %Live trees
[num2,txt2,raw2] = xlsread([path,'\tow_06_mort_dbh.xls']); %Dead trees
[num3,txt3,raw3] = xlsread([path,'\tow_06_rcrt_dbh.xls']); %Recruitment
[num4,txt4,raw4] = xlsread([path,'\tow_sapling_06.xls']); %Saplings
dbh = [num1(:,end) ; num2(:,end) ; num3(:,end) ; num4(:,end)];
basal_area = pi.*(((dbh./2)./100).^2);
plots = [txt1(2:end,1) ; txt2(2:end,1) ; txt3(2:end,1) ; txt4(2:end,1)];
plots2 = unique(plots);

for x=1:length(plots2)
    ind = find(strcmp(plots,plots2(x,:)));
    basal_area_per_plot(x,1) = sum(basal_area(ind,1));
end
basal_area_per_hectare = basal_area_per_plot.*10000./(pi.*10^2);
basal_area_per_hectare = mean(basal_area_per_hectare);
basal_area_cover = basal_area_per_hectare./10000; %cover by basal area (100% = 1)



%%%%%%%%%%%%%%%%%%%%%%%
% Surface area covered by rocks
[num,txt,raw] = xlsread([path,'\rocks (Land Use and Forest Dynamics at Harvard Forest).xls']);
rocks = num(:,6);

%Nb of plots per rock cover category
nb_plots=[];
for x=0:7
    nb_plots(x+1,1) = length(find(rocks==x));
end

%Mean percent cover of rocks (using mid-point of each category)
midpoint = [0 0.5 2 4 10 20 37.5 62.5]';
rock_cover = (sum(nb_plots.*midpoint)./sum(nb_plots))/100; %cover by rocks (100% = 1)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate percent cover by vegetation type in footprint (90% contribution)
%.xls files contain GIA data. Footprint extent was computed using Kjlun et
%al. (2004) model for period with and without snow on the ground. Forest
%stands within the footprints were identified using a GIS map of Harvard
%Forest available in the Harvard forest data archive.
%EMS
% No snow on the ground
% footprint_area_nosnow = xlsread([path,'\Harvard_Forest_ArcView\90prct_nosnow_footprint.xls'],'D2:D2');
[num,txt,raw] = xlsread([path,'\90prct_nosnow_data.xls'],'A2:AE63');
stand_area_nosnow = sum(num(:,end));
indC = find(strcmp('C',txt(:,7)));
indD = find(strcmp('D',txt(:,7)));
indM = find(strcmp('MD',txt(:,7)) | strcmp('MC',txt(:,7)));
indP = find(strcmp('P',txt(:,7)));
indW = find(strcmp('W',txt(:,7)));
veg_type_area_nosnow = [sum(num(indD,end)) ; sum(num(indC,end)) ; sum(num(indM,end)) ; sum(num(indP,end)) ; sum(num(indW,end))];%deciduous, conifers, mixed, pines, wetland
veg_type_area_percent_nosnow = veg_type_area_nosnow./sum(veg_type_area_nosnow); %Percent area of footprint for each vegetation cover type %Use in Harvard Forest Archive (file Vegetation_type_cover_in_footprint.csv, row 4, EMS site)

% Snow on the ground
% footprint_area_snow = xlsread([path,'\Harvard_Forest_ArcView\90prct_snow_footprint.xls'],'D2:D2');
[num,txt,raw] = xlsread([path,'\90prct_nosnow_data.xls'],'A2:AE51');
stand_area_snow = sum(num(:,end));
indC = find(strcmp('C',txt(:,7)));
indD = find(strcmp('D',txt(:,7)));
indM = find(strcmp('MD',txt(:,7)) | strcmp('MC',txt(:,7)));
indP = find(strcmp('P',txt(:,7)));
indW = find(strcmp('W',txt(:,7)));
veg_type_area_snow = [sum(num(indD,end)) ; sum(num(indC,end)) ; sum(num(indM,end)) ; sum(num(indP,end)) ; sum(num(indW,end))];%deciduous, conifers, mixed, pines, wetland
veg_type_area_percent_snow = veg_type_area_snow./sum(veg_type_area_snow); %Percent area of footprint for each vegetation cover type %Use in Harvard Forest Archive (file Vegetation_type_cover_in_footprint.csv, row 2, EMS site)

% Intermittent snow cover (compute average of snow and nosnow)
veg_type_area_percent_inter_snow = mean([veg_type_area_percent_snow veg_type_area_percent_nosnow],2); %Use in Harvard Forest Archive (file Vegetation_type_cover_in_footprint.csv, row 3, EMS site)


%HEM
% No snow on the ground
[num,txt,raw] = xlsread([path,'\90prct_nosnow_data_HEM.xls'],'A2:AE40');
stand_area_nosnow_hem = sum(num(:,end));
indC = find(strcmp('C',txt(:,7)));
indD = find(strcmp('D',txt(:,7)));
indM = find(strcmp('MD',txt(:,7)) | strcmp('MC',txt(:,7)));
indP = find(strcmp('P',txt(:,7)));
indW = find(strcmp('W',txt(:,7)));
veg_type_area_nosnow_hem = [sum(num(indD,end)) ; sum(num(indC,end)) ; sum(num(indM,end)) ; sum(num(indP,end)) ; sum(num(indW,end))];%deciduous, conifers, mixed, pines, wetland
veg_type_area_percent_nosnow_hem = veg_type_area_nosnow_hem./sum(veg_type_area_nosnow_hem); %Percent area of footprint for each vegetation cover type %Use in Harvard Forest Archive (file Vegetation_type_cover_in_footprint.csv, row 4, HEM site)

% Snow on the ground
[num,txt,raw] = xlsread([path,'\90prct_nosnow_data_HEM.xls'],'A2:AE38');
stand_area_snow_hem = sum(num(:,end));
indC = find(strcmp('C',txt(:,7)));
indD = find(strcmp('D',txt(:,7)));
indM = find(strcmp('MD',txt(:,7)) | strcmp('MC',txt(:,7)));
indP = find(strcmp('P',txt(:,7)));
indW = find(strcmp('W',txt(:,7)));
veg_type_area_snow_hem = [sum(num(indD,end)) ; sum(num(indC,end)) ; sum(num(indM,end)) ; sum(num(indP,end)) ; sum(num(indW,end))];%deciduous, conifers, mixed, pines, wetland
veg_type_area_percent_snow_hem = veg_type_area_snow_hem./sum(veg_type_area_snow_hem); %Percent area of footprint for each vegetation cover type %Use in Harvard Forest Archive (file Vegetation_type_cover_in_footprint.csv, row 2, HEM site)

% Intermittent snow cover (compute average of snow and nosnow)
veg_type_area_percent_inter_snow_hem = mean([veg_type_area_percent_snow_hem veg_type_area_percent_nosnow_hem],2); %Use in Harvard Forest Archive (file Vegetation_type_cover_in_footprint.csv, row 3, HEM site)




%%%%%%%%%%%%%%%%%%
%Compute daily and monthly Re for HEM tower (EMS done earlier in the script)
HEM_R_gfNEE_backup = HEM_R_gfNEE;
HEM_NEE_gfNEE_backup = HEM_NEE_gfNEE;
HEM_GEE_gfNEE_backup = HEM_GEE_gfNEE;
HEM_R_gfNEE = HEM_R_gfNEE.*(3600*24*12/1000000); %Convert umol m-2 s-1 to g m-2 d-1
HEM_NEE_gfNEE = HEM_NEE_gfNEE.*(3600*24*12/1000000); %Convert umol m-2 s-1 to g m-2 d-1
HEM_GEE_gfNEE = HEM_GEE_gfNEE.*(3600*24*12/1000000); %Convert umol m-2 s-1 to g m-2 d-1

%Daily and monthly sums with original PI gap-filling
y=0;
HEM_daily_original_GF=[];
for x=2004:2007 %no Re data in 2000-01, and none available after 2007
    for xx=1:366
        ind = find(timestamps_gfNEE(:,1)==x & timestamps_gfNEE(:,3)==xx);
        if ~isempty(ind)
            y=y+1;
            try
                ind = ind+1;
                HEM_daily_original_GF(y,:) = [x xx nanmean(HEM_R_gfNEE(ind)) nanmean(HEM_NEE_gfNEE(ind)) nanmean(HEM_GEE_gfNEE(ind))];
            catch
                ind = ind(2:end)-1;
                HEM_daily_original_GF(y,:) = [x xx nanmean(HEM_R_gfNEE(ind)) nanmean(HEM_NEE_gfNEE(ind)) nanmean(HEM_GEE_gfNEE(ind))];
            end
        end
    end
end
y=0;
HEM_monthly_original_GF=[];
for x=2004:2007
    for xx=1:12
        nb_days = fn_nb_days_in_month(x,xx);
        y=y+1;
        [dayi,dayj]=fn_month_to_jday(x,xx);
        ind = find(HEM_daily_original_GF(:,1)==x & HEM_daily_original_GF(:,2)>=dayi & HEM_daily_original_GF(:,2)<=dayj);
        ind2 = find(HEM_daily_original_GF(:,1)==x & HEM_daily_original_GF(:,2)>=dayi & HEM_daily_original_GF(:,2)<=dayj & isnan(HEM_daily_original_GF(:,3)));
        HEM_monthly_original_GF(y,:) = [x xx nanmean(HEM_daily_original_GF(ind,4)).*nb_days nanmean(HEM_daily_original_GF(ind,3)).*nb_days nanmean(HEM_daily_original_GF(ind,5)).*nb_days length(ind2)];
    end
end
%Delete data when there were not enough days with data in month
ind = find(HEM_monthly_original_GF(:,end)>10);
HEM_monthly_original_GF(ind,4) = NaN;
HEM_monthly_original_GF(:,end)=[]; %Use in Harvard Forest Archive (column 4 is Re-HEM_PI in file Re_Rs_monthly.csv)

%%%
%HEM tower R - fill gaps using Fluxnet-Canada (FCRN) gap-filling algorithm
[HEM_Rgf_FCRN,HEM_FCRN_timestamps] = fn_HEM_Re_gapfill(Ts,Ts_timestamps);
HEM_Rgf_FCRN = HEM_Rgf_FCRN.*(3600*24*12/1000000); %Convert umol m-2 s-1 to g m-2 d-1

%Daily and monthly sums with FCRN-gap-filled data
y=0;
HEM_daily_original_GF_FCRN=[];
for x=2004:2009
    for xx=1:366
        ind = find(HEM_FCRN_timestamps(:,1)==x & HEM_FCRN_timestamps(:,2)==xx);
        if ~isempty(ind)
            y=y+1;
            try
                ind = ind+1;
                HEM_daily_original_GF_FCRN(y,:) = [x xx nanmean(HEM_Rgf_FCRN(ind)) nanmean(HEM_Rgf_FCRN(ind))*NaN nanmean(HEM_Rgf_FCRN(ind))*NaN];
            catch
                ind = ind(2:end)-1;
                HEM_daily_original_GF_FCRN(y,:) = [x xx nanmean(HEM_Rgf_FCRN(ind)) nanmean(HEM_Rgf_FCRN(ind))*NaN nanmean(HEM_Rgf_FCRN(ind))*NaN];
            end
        elseif ~(xx==366 & (x==2005 | x==2006 | x==2007 | x==2009))
            y=y+1;
            ind = ind+1;
            HEM_daily_original_GF_FCRN(y,:) = [x xx NaN NaN NaN];
        end
    end
end
y=0;
HEM_monthly_original_GF_FCRN=[];
for x=2004:2009
    for xx=1:12
        nb_days = fn_nb_days_in_month(x,xx);
        y=y+1;
        [dayi,dayj]=fn_month_to_jday(x,xx);
        ind = find(HEM_daily_original_GF_FCRN(:,1)==x & HEM_daily_original_GF_FCRN(:,2)>=dayi & HEM_daily_original_GF_FCRN(:,2)<=dayj);
        ind2 = find(HEM_daily_original_GF_FCRN(:,1)==x & HEM_daily_original_GF_FCRN(:,2)>=dayi & HEM_daily_original_GF_FCRN(:,2)<=dayj & isnan(HEM_daily_original_GF_FCRN(:,3)));
        HEM_monthly_original_GF_FCRN(y,:) = [x xx nanmean(HEM_daily_original_GF_FCRN(ind,4)).*nb_days nanmean(HEM_daily_original_GF_FCRN(ind,3)).*nb_days nanmean(HEM_daily_original_GF_FCRN(ind,5)).*nb_days length(ind2)];
    end
end
ind = find(HEM_monthly_original_GF_FCRN(:,1)==2004 & HEM_monthly_original_GF_FCRN(:,2)==6);
HEM_monthly_original_GF_FCRN(ind,4)=NaN;
HEM_monthly_original_GF_FCRN(:,end)=[]; %Use for Harvard Forest Archive (column 4 is Re-HEM_FCRN in file Re_Rs_monthly.csv)

%Merge both series, and overwrite bad original gap-filling with FCRN-gap-filled data
%Compute daily, monthly and yearly series
temp = HEM_daily_original_GF;
ind = find(temp(:,1)==2006 & temp(:,2)>=188 & temp(:,2)<=279); %Period with bad original gap-filling
temp(ind,3:end) = NaN;
HEM_daily_original_GF_FCRN2 = HEM_daily_original_GF_FCRN*NaN;
HEM_daily_original_GF_FCRN2(:,1:3) = HEM_daily_original_GF_FCRN(:,1:3);
for x=1:length(temp)
    if ~isnan(temp(x,3))
        ind = find(HEM_daily_original_GF_FCRN2(:,1)==temp(x,1) & HEM_daily_original_GF_FCRN2(:,2)==temp(x,2));
        HEM_daily_original_GF_FCRN2(ind,3) = temp(x,3); %Use for Harvard Forest Archive (column 3 is HEM-Re in file Re_Rs_daily.csv)
    end
end
y=0;
HEM_monthly_original_GF_FCRN2=[];
for x=2004:2009
    for xx=1:12
        nb_days = fn_nb_days_in_month(x,xx);
        y=y+1;
        [dayi,dayj]=fn_month_to_jday(x,xx);
        ind = find(HEM_daily_original_GF_FCRN2(:,1)==x & HEM_daily_original_GF_FCRN2(:,2)>=dayi & HEM_daily_original_GF_FCRN2(:,2)<=dayj);
        ind2 = find(HEM_daily_original_GF_FCRN2(:,1)==x & HEM_daily_original_GF_FCRN2(:,2)>=dayi & HEM_daily_original_GF_FCRN2(:,2)<=dayj & isnan(HEM_daily_original_GF_FCRN2(:,3)));
        HEM_monthly_original_GF_FCRN2(y,:) = [x xx nanmean(HEM_daily_original_GF_FCRN2(ind,4)).*nb_days nanmean(HEM_daily_original_GF_FCRN2(ind,3)).*nb_days nanmean(HEM_daily_original_GF_FCRN2(ind,5)).*nb_days length(ind2)];
    end
end
ind = find(HEM_monthly_original_GF_FCRN2(:,1)==2004 & HEM_monthly_original_GF_FCRN2(:,2)==6);
HEM_monthly_original_GF_FCRN2(ind,4)=NaN;
HEM_monthly_original_GF_FCRN2(:,end)=[];

y=0;
HEM_yearly_original_GF_FCRN=[];
for x=1992:2009
    y=y+1;
    ind = find(HEM_monthly_original_GF_FCRN2(:,1)==x);
    HEM_yearly_original_GF_FCRN(y,:) = [x sum(HEM_monthly_original_GF_FCRN2(ind,4))];
end
HEM_yearly_original_GF_FCRN(HEM_yearly_original_GF_FCRN==0) = NaN; %Use for Harvard Forest Archive (HEM-Re in file Re_Rs_annual.csv)
 
 
%%%%%%%%%%%%%%%%%%
%Compute monthly Rs for each experiment and vegetation type (controls only) - 18-year time series

%Get all data needed for calculations
ind = find(data(:,7)==6 & data(:,23)<=4);
data(ind,16) = 1; % upland
ind = find(data(:,7)==6 & data(:,23)>=5);
data(ind,16) = 7; % wetland margin
data2 = data;
ind = find(isnan(data2(:,5)));
data2(ind,:)=[]; %Delete experiments with no timestamps
ind = find(data2(:,8)==19);
data2(ind,:)=[]; %Delete experiment with only 3 data points
%Also, there is no DIRT data (no timestamps available)


%Compute number of studies available for each month of each year
y=0;
for x=1991:2009
    for xx=1:12
        y=y+1;
        ind = find(data2(:,1)==x & data2(:,2)==xx);
        temp = data2(ind,[7 8]);
        temp2 = unique(temp,'rows');
        nb_studies(y,:) = [x xx size(temp2,1)]; %Use for Harvard Forest archive (file Re_Rs_monthly.csv)
    end
end



Rs_experiments_veg_type = unique([data2(:,7) data2(:,8) data2(:,22) data2(:,16)],'rows');
ind = find(Rs_experiments_veg_type(:,4)>=8);
Rs_experiments_veg_type(ind,:)=[]; %Delete logged forest types
option = 0; %0=no figures, 1=figures


%Compute daily and monthly Rs for each experiment and vegetation type (controls only) - 18-year time series
Rs_hh_for_Fig2 = [ones(length(Ts_timestamps),length(Rs_experiments_veg_type))*NaN];
for x_to_do=1:length(Rs_experiments_veg_type)
    [Deciduous_Rs_raw_B_reserved,temp_Rs_Deciduous,Deciduous_Rs_raw_B,Ts_ref_Deciduous,Rs_hh_Deciduous,R10_Dec,Q10_Dec] = Script_Rs_estimate_for_figures_all_studies_separated_bias_corr(x_to_do,Rs_experiments_veg_type,data2,Ts,Ts_timestamps,option);
    Rs_hh_for_Fig2(:,x_to_do) = Rs_hh_Deciduous;
end

Rs_data_winter_control=[]; %Jan-Mar only
for x=1:length(Rs_experiments_veg_type)
    ind = find(data2(:,7)==Rs_experiments_veg_type(x,1) & data2(:,8)==Rs_experiments_veg_type(x,2) & data2(:,16)==Rs_experiments_veg_type(x,4) & (data2(:,20)==0 | data2(:,20)==2) & (data2(:,2)<=3));
    temp_data2 = data2(ind,:);
    for xx=1:size(temp_data2,1)
        ind2 = find(Ts_timestamps(:,1)==temp_data2(xx,1) & Ts_timestamps(:,3)==temp_data2(xx,4) & Ts_timestamps(:,4)==temp_data2(xx,5) & Ts_timestamps(:,5)==temp_data2(xx,6));
        Rs_data_winter_control = [Rs_data_winter_control ; [Rs_experiments_veg_type(x,[1 2 4]) Ts_timestamps(ind2,:) temp_data2(xx,28) Rs_hh_for_Fig2(ind2,x) Ts(ind2,1)]]; %PI study treatment year month jday hh mn measured_Rs modeled_Rs Ts %Use columns 9-10 for Harvard Forest Archive (file Rs_winter_measured_vs_modeled.csv)
    end
end

days = unique(Ts_timestamps(:,[1 3]),'rows');
days(end,:) = [];
Rs_data_controls_veg_type_daily = [days ones(length(days),size(Rs_hh_for_Fig2,2))*NaN];
for x=1:length(days)
    ind = find(Ts_timestamps(:,1)==days(x,1) & Ts_timestamps(:,3)==days(x,2));
    ind = ind+1;
    Rs_data_controls_veg_type_daily(x,3:end) = sum(Rs_hh_for_Fig2(ind,:),1);
end
Rs_data_controls_veg_type_daily(:,3:end) = Rs_data_controls_veg_type_daily(:,3:end).*(1800.*12./1000000); %convert units to gC m-2 d-1


%Calculate average daily sums per vegetation type
ind = find(Rs_experiments_veg_type(:,4)==1);
Deciduous_Rs_daily = mean(Rs_data_controls_veg_type_daily(:,ind+2),2);
ind = find(Rs_experiments_veg_type(:,4)==2);
Conifers_Rs_daily = mean(Rs_data_controls_veg_type_daily(:,ind+2),2);
ind = find(Rs_experiments_veg_type(:,4)==3);
Mixed_Rs_daily = mean(Rs_data_controls_veg_type_daily(:,ind+2),2);
ind = find(Rs_experiments_veg_type(:,4)==4);
Pine_Rs_daily = mean(Rs_data_controls_veg_type_daily(:,ind+2),2);
ind = find(Rs_experiments_veg_type(:,4)==5 | Rs_experiments_veg_type(:,4)==6 | Rs_experiments_veg_type(:,4)==7);
Wetland_Rs_daily = mean(Rs_data_controls_veg_type_daily(:,ind+2),2);


%Compute Rs annual (and mean) per vegetation type
temp = [Deciduous_Rs_daily Conifers_Rs_daily Mixed_Rs_daily Pine_Rs_daily Wetland_Rs_daily]; %Use in Harvard Forest Archive (file Rs_vegetation_type_daily.csv)
Rs_annual_veg_type = [];
Rs_annual_veg_type_mean = [];
for x=1992:2009
    ind = find(Rs_data_controls_veg_type_daily(:,1)==x);
    Rs_annual_veg_type = [Rs_annual_veg_type ; sum(temp(ind,:),1)]; %Use in Harvard Forest Archive (columns 2-6 in file Rs_vegetation_type_annual.csv)
end
Rs_annual_veg_type = Rs_annual_veg_type.*[1-(basal_area_cover + rock_cover)]; %Adjust for rocks and basal area %Use in Harvard Forest Archive (columns 7-11 in file Rs_vegetation_type_annual.csv)
Rs_annual_veg_type_mean = [mean(Rs_annual_veg_type,1) ; std(Rs_annual_veg_type,1)]';

%Weigh Rs according to vegetation type cover in flux tower footprint
%Multiply by percent cover for snow, nosnow, intermittent snow periods
%snow=1; no snow=2; intermittent snow=3
snow_code=ones(length(Rs_data_controls_veg_type_daily),1)*3;
ind = find(Rs_data_controls_veg_type_daily(:,2)<=95 | Rs_data_controls_veg_type_daily(:,2)>=340); %Snowy period
snow_code(ind,1) = 1;
ind = find(Rs_data_controls_veg_type_daily(:,2)>=113 & Rs_data_controls_veg_type_daily(:,2)<=289); %Snow-free period
snow_code(ind,1) = 2;

ind1 = find(snow_code==1); %snow
ind2 = find(snow_code==2); %nosnow
ind3 = find(snow_code==3); %intermittent snow
%EMS
Rs_daily(ind1,1:5) = [Deciduous_Rs_daily(ind1,1).*veg_type_area_percent_snow(1,1) Conifers_Rs_daily(ind1,1).*veg_type_area_percent_snow(2,1) ...
    Mixed_Rs_daily(ind1,1).*veg_type_area_percent_snow(3,1) Pine_Rs_daily(ind1,1).*veg_type_area_percent_snow(4,1) ...
    Wetland_Rs_daily(ind1,1).*veg_type_area_percent_snow(5,1)];
Rs_daily(ind2,1:5) = [Deciduous_Rs_daily(ind2,1).*veg_type_area_percent_nosnow(1,1) Conifers_Rs_daily(ind2,1).*veg_type_area_percent_nosnow(2,1) ...
    Mixed_Rs_daily(ind2,1).*veg_type_area_percent_nosnow(3,1) Pine_Rs_daily(ind2,1).*veg_type_area_percent_nosnow(4,1) ...
    Wetland_Rs_daily(ind2,1).*veg_type_area_percent_nosnow(5,1)];
Rs_daily(ind3,1:5) = [Deciduous_Rs_daily(ind3,1).*veg_type_area_percent_inter_snow(1,1) Conifers_Rs_daily(ind3,1).*veg_type_area_percent_inter_snow(2,1) ...
    Mixed_Rs_daily(ind3,1).*veg_type_area_percent_inter_snow(3,1) Pine_Rs_daily(ind3,1).*veg_type_area_percent_inter_snow(4,1) ...
    Wetland_Rs_daily(ind3,1).*veg_type_area_percent_inter_snow(5,1)];

%HEM
Rs_daily_hem(ind1,1:5) = [Deciduous_Rs_daily(ind1,1).*veg_type_area_percent_snow_hem(1,1) Conifers_Rs_daily(ind1,1).*veg_type_area_percent_snow_hem(2,1) ...
    Mixed_Rs_daily(ind1,1).*veg_type_area_percent_snow_hem(3,1) Pine_Rs_daily(ind1,1).*veg_type_area_percent_snow_hem(4,1) ...
    Wetland_Rs_daily(ind1,1).*veg_type_area_percent_snow_hem(5,1)];
Rs_daily_hem(ind2,1:5) = [Deciduous_Rs_daily(ind2,1).*veg_type_area_percent_nosnow_hem(1,1) Conifers_Rs_daily(ind2,1).*veg_type_area_percent_nosnow_hem(2,1) ...
    Mixed_Rs_daily(ind2,1).*veg_type_area_percent_nosnow_hem(3,1) Pine_Rs_daily(ind2,1).*veg_type_area_percent_nosnow_hem(4,1) ...
    Wetland_Rs_daily(ind2,1).*veg_type_area_percent_nosnow_hem(5,1)];
Rs_daily_hem(ind3,1:5) = [Deciduous_Rs_daily(ind3,1).*veg_type_area_percent_inter_snow_hem(1,1) Conifers_Rs_daily(ind3,1).*veg_type_area_percent_inter_snow_hem(2,1) ...
    Mixed_Rs_daily(ind3,1).*veg_type_area_percent_inter_snow_hem(3,1) Pine_Rs_daily(ind3,1).*veg_type_area_percent_inter_snow_hem(4,1) ...
    Wetland_Rs_daily(ind3,1).*veg_type_area_percent_inter_snow_hem(5,1)];


%Compute Rs ref control for all cover types together
Rs_daily = sum(Rs_daily,2); %EMS
Rs_daily_hem = sum(Rs_daily_hem,2); %HEM

%Adjust for rocks and basal area
Rs_daily2 = Rs_daily.*[1-(basal_area_cover + rock_cover)]; %EMS %Use for Harvard Forest Archive (EMS-Rs in file Re_Rs_daily.csv)
Rs_daily2_hem_backup = Rs_daily_hem.*[1-(basal_area_cover + rock_cover)]; %HEM

%HEM - Keep only years with data
indi = find(Rs_data_controls_veg_type_daily(:,1)==2004 & Rs_data_controls_veg_type_daily(:,2)==1);
indj = find(Rs_data_controls_veg_type_daily(:,1)==2007 & Rs_data_controls_veg_type_daily(:,2)==365);
Rs_daily2_hem = Rs_daily2_hem_backup(indi:indj,:);
indi = find(Rs_data_controls_veg_type_daily(:,1)==2004 & Rs_data_controls_veg_type_daily(:,2)==1);
indj = find(Rs_data_controls_veg_type_daily(:,1)==2009 & Rs_data_controls_veg_type_daily(:,2)==365);
Rs_daily2_hem_FCRN = Rs_daily2_hem_backup(indi:indj,:); %Use for Harvard Forest Archive (HEM-Rs in file Re_Rs_daily.csv)


%Compute daily Tsoil and Tair 
Ts_daily=[];
for x=1:366
    ind = find(timestamps_gfNEE(:,3)==x);
    Ts_daily(x,1) = mean(Ts_gfNEE(ind)); %Use in Harvard Forest Archive (file Ts_Ta_daily.csv)
    Ta_daily(x,1) = mean(Ta_gfNEE(ind)); %Use in Harvard Forest Archive (file Ts_Ta_daily.csv)
end

%%%%%%%%%%
%Compute Rs ref MEDIAN for each day of the year
Rs_daily_median=[];
Re_daily_median=[];
Rs_daily_hem_median=[];
Rs_daily_hem_FCRN_median=[];
Re_daily_hem_median=[];
Re_daily_hem_FCRN_median=[];
for x=1:366
    ind = find(Rs_data_controls_veg_type_daily(:,2)==x & Rs_data_controls_veg_type_daily(:,1)>=1996 & Rs_data_controls_veg_type_daily(:,1)<=2009); %1996-2009
    ind4 = find(HEM_daily_original_GF_FCRN2(:,2)==x);
    Rs_daily_median(x,1) = median(Rs_daily2(ind,1)); %Rs ref control %Use for Harvard Forest Archive (EMS-Rs in file Re_Rs_Raboveground_daily_median.csv)
    Re_daily_median(x,1) = median(EMS_daily_original_GF(ind,4)); %Re EMS tower %Use for Harvard Forest Archive (EMS-Re in file Re_Rs_Raboveground_daily_median.csv)
    Rs_daily_hem_FCRN_median(x,1) = median(Rs_daily2_hem_FCRN(ind4,1)); %Rs ref control %Use for Harvard Forest Archive (HEM-Rs in file Re_Rs_Raboveground_daily_median.csv)
    Re_daily_hem_FCRN_median(x,1) = nanmedian(HEM_daily_original_GF_FCRN2(ind4,3)); %Re HEM tower_FCRN %Use for Harvard Forest Archive (HEM-Re in file Re_Rs_Raboveground_daily_median.csv)
end

%Compute Raboveground
R_aboveground_median = Re_daily_median - Rs_daily_median; %EMS %Use for Harvard Forest Archive (EMS-Raboveground in file Re_Rs_Raboveground_daily_median.csv)
R_aboveground_hem_FCRN_median = Re_daily_hem_FCRN_median - Rs_daily_hem_FCRN_median; %HEM-FCRN %Use for Harvard Forest Archive (HEM-Raboveground in file Re_Rs_Raboveground_daily_median.csv)



%%%%%%%%%%%%
%Compute monthly Rs totals
%EMS
months = ones(length(Rs_daily2),1)*NaN;
for x=1:length(Rs_daily2)
    months(x,1) = fn_jday_to_month(Rs_data_controls_veg_type_daily(x,1),Rs_data_controls_veg_type_daily(x,2));
end
months2 = unique([Rs_data_controls_veg_type_daily(:,1) months],'rows');
months2(1,:)=[];

Rs_monthly_for_Fig2a=[];
for x=1:length(months2)
    ind = find(Rs_data_controls_veg_type_daily(:,1)==months2(x,1) & months(:,1)==months2(x,2));
    nb_days = fn_nb_days_in_month(months2(x,1),months2(x,2));
    Rs_monthly_for_Fig2a(x,1) = mean(Rs_daily2(ind,1)).*nb_days; %Use for Harvard Forest Archive (EMS-Rs in file Re_Rs_monthly.csv)
end

Rs_yearly=[];
y=0;
for x=1992:2009
    y=y+1;
    ind = find(Rs_data_controls_veg_type_daily(:,1)==x);
    Rs_yearly(y,:) = [x sum(Rs_daily2(ind,1))]; %Use for Harvard Forest Archive (EMS-Rs in file Re_Rs_annual.csv)
end
%Rs/Re ratio
Rs_Re_ratio_yearly = Rs_yearly(:,2)./EMS_yearly_original_GF(:,2); %Use for Harvard Forest Archive (EMS-Rs/EMS-Re ratio in file Re_Rs_annual.csv)



%HEM with FCRN gap-filling
months = ones(length(Rs_daily2_hem_FCRN),1)*NaN;
for x=1:length(Rs_daily2_hem_FCRN)
    months(x,1) = fn_jday_to_month(HEM_daily_original_GF_FCRN2(x,1),HEM_daily_original_GF_FCRN2(x,2));
end
months2 = unique([HEM_daily_original_GF_FCRN2(:,1) months],'rows');

Rs_monthly_for_Fig2c_FCRN=[];
for x=1:length(months2)
    ind = find(HEM_daily_original_GF_FCRN2(:,1)==months2(x,1) & months(:,1)==months2(x,2));
    nb_days = fn_nb_days_in_month(months2(x,1),months2(x,2));
    Rs_monthly_for_Fig2c_FCRN(x,1) = mean(Rs_daily2_hem_FCRN(ind,1)).*nb_days; %Use for Harvard Forest Archive (HEM-Rs in file Re_Rs_monthly.csv)
end

Rs_yearly_FCRN=[];
y=0;
for x=1992:2009
    y=y+1;
    ind = find(HEM_daily_original_GF_FCRN2(:,1)==x);
    Rs_yearly_FCRN(y,:) = [x sum(Rs_daily2_hem_FCRN(ind,1))];
end
Rs_yearly_FCRN(Rs_yearly_FCRN==0)=NaN; %Use for Harvard Forest Archive (HEM-Rs in file Re_Rs_annual.csv)
%Rs/Re ratio
Rs_Re_ratio_yearly_FCRN = Rs_yearly_FCRN(:,2)./HEM_yearly_original_GF_FCRN(:,2); %Use for Harvard Forest Archive (HEM-Rs/HEM-Re ratio in file Re_Rs_annual.csv)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rs from control plots from all experiments
all_controls = Nmin_Rs_data2(:,[3 6 12 15 18 21 24 27 30 33 42 45 48 51 54 60 66 72 78 81 132 151 155]);
all_controls(all_controls==0)=NaN;

%Compute monthly sums
for x=1:size(controls,1)
    nb_days = fn_nb_days_in_month(Nmin_Rs_data(x,1),Nmin_Rs_data(x,2));
    all_controls(x,:) = all_controls(x,:).*nb_days;
end
all_controls(all_controls==0)=NaN;
    
y=0;
Rs_control_summer=[];
Rs_treatment_summer=[];
Rs_all_control_summer=[];
for x=1991:2009
    y=y+1;
    ind2 = find(Nmin_Rs_data(:,1)==x & Nmin_Rs_data(:,2)>=4 & Nmin_Rs_data(:,2)<=10);
    Rs_control_summer(y,:) = [x sum(controls(ind2,:))]; %Summer with all months present only
    Rs_treatment_summer(y,:) = [x sum(treatments(ind2,:))]; %Summer with all months present only
    Rs_all_control_summer(y,:) = [x sum(all_controls(ind2,:))]; %Summer with all months present only
end
Rs_control_summer(isinf(Rs_control_summer))=NaN;
Rs_treatment_summer(isinf(Rs_treatment_summer))=NaN;
Rs_all_control_summer(isinf(Rs_all_control_summer))=NaN;
Rs_control_summer(Rs_control_summer==0)=NaN;
Rs_treatment_summer(Rs_treatment_summer==0)=NaN;
Rs_all_control_summer(Rs_all_control_summer==0)=NaN;

%%%%%
Rs_control_summer = [Rs_control_summer(:,1:4) Rs_control_summer(:,8) Rs_control_summer(:,[5:7 9:end])]; %Change study order to drydown / Davidson-Simes / Munger / Barre Woods / ProspectHill / Warm+N / DIRT / ChronicN
Rs_treatment_summer = [Rs_treatment_summer(:,1:4) Rs_treatment_summer(:,8) Rs_treatment_summer(:,[5:7 9:end])]; %Change study order to drydown / Davidson-Simes / Munger / Barre Woods / ProspectHill / Warm+N / DIRT / ChronicN
Rs_control_summer(1,:)=[];
Rs_treatment_summer(1,:)=[];
Rs_control_summer(:,[19 22])=[]; %Delete sulphur treatments from ChronicN, which we don't use
Rs_treatment_summer(:,[19 22])=[]; %Delete sulphur treatments from ChronicN, which we don't use
Rs_control_summer(12:13,2:3)=NaN;
Rs_control_summer(11:13,4)=NaN;
Rs_control_summer(11,6)=NaN;
Rs_treatment_summer(12:13,2:3)=NaN;
Rs_treatment_summer(11:13,4)=NaN;
Rs_treatment_summer(11,6)=NaN;

%Treatments and corresponding controls
temp = [repmat(Rs_control_summer(:,1),size(Rs_control_summer,2)-1,1) reshape(Rs_treatment_summer(:,2:end),[],1) reshape(Rs_control_summer(:,2:end),[],1)];
treat = reshape(repmat([4 5 8 8 7 7 3 7 9 10 11 12 13 14 15 9 9 9 9],18,1),[],1);
treat2 = reshape(repmat([1 2 2 2 3 3 7 3 4 6 6 5 5 5 5 4 4 4 4],18,1),[],1); %1=drydown 2=girdled/logged 3=warm 4=fert 5=DIRT (except double litter) 6=fert+warm/double litter 7=dist.ctrl
study = reshape(repmat([3 15 15 4 1 14 14 14 14 14 23 23 23 23 23 24 24 24 24],18,1),[],1); %see list of all variable to get codes
PI = reshape(repmat([1 1 1 5 4 4 4 7 7 7 8 8 8 8 8 7 7 7 7],18,1),[],1);
ind = find(isnan(temp(:,2)) | isnan(temp(:,3)));
temp2 = [temp(:,1) PI study treat treat2 temp(:,2:3)];
temp2(ind,:) = [];

%All controls
Rs_all_control_summer(1,:)=[];
temp = [repmat(Rs_all_control_summer(:,1),size(Rs_all_control_summer,2)-1,1) reshape(Rs_all_control_summer(:,2:end),[],1)];
treat = ones(length(temp),1)*2; %all controls
study = reshape(repmat([2 3 5 6 7 10 11 12 13 15 16 17 18 8 9 1 14 4 5 9 23 24 24],18,1),[],1);
PI = reshape(repmat([1 1 1 1 1 1 1 1 1 1 1 1 1 3 3 4 4 5 6 7 8 7 7],18,1),[],1);
ind = find(isnan(temp(:,2)));
temp = [temp(:,1) PI study treat treat.*0 temp(:,2) temp(:,2).*NaN];
temp(ind,:) = [];
SAS_data = [temp2 ; temp];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Prepare data for comparison between climate and treatment effect

%Compute Rs for each control (study/site and vegetation cover) with annual Ts regressions
option = 0; %0=no figures, 1=figures
Rs_data_veg_type = [sort(repmat([1991:2009]',12,1)) repmat([1:12]',2009-1991+1,1)];
Rs_data_veg_type = [Rs_data_veg_type ones(length(Rs_data_veg_type),length(Rs_experiments_veg_type).*3)*NaN];
col = 2;
var_pred_for_Rs_MAT_calculation = [[1991:2009]' ones([2009-1991+1],length(Rs_experiments_veg_type)).*NaN];
for exp_ind=1:length(Rs_experiments_veg_type)
    [Rs_data_processed_veg_type,var_pred_temp] = fn_compute_Rs_veg_type_Harvard_Forest_bias_corr(exp_ind,data2,Ts,Ts_timestamps,Rs_experiments_veg_type,option,path);
    [c,ia,ib] = intersect(Rs_data_veg_type(:,1:2),Rs_data_processed_veg_type(:,1:2),'rows');
    Rs_data_veg_type(ia,col+1:col+3) = Rs_data_processed_veg_type(ib,[7:9]);
    col = col+3;
    [c,ia,ib] = intersect(var_pred_for_Rs_MAT_calculation(:,1),var_pred_temp(:,1));
    var_pred_for_Rs_MAT_calculation(ia,exp_ind+1) = var_pred_temp(ib,2);
end
fn_delete_Excel_default_sheets([path,'\'],'Harvard_Forest_soil_respiration_veg_type_final.xls');
fn_delete_Excel_default_sheets([path,'\'],'Harvard_Forest_soil_respiration_veg_type_Q10_final.xls');

Rs_data_veg_type2 = [Rs_data_veg_type(:,1:2) Rs_data_veg_type(:,4:3:size(Rs_data_veg_type,2))];
Ts_data_veg_type2 = [Rs_data_veg_type(:,1:2) Rs_data_veg_type(:,5:3:size(Rs_data_veg_type,2))];
%Append DIRT data
[c,ia,ib] = intersect(Rs_DIRT_monthly(:,1:2),Rs_data_veg_type2(:,1:2),'rows');
temp = ones(size(Rs_data_veg_type2,1),1)*NaN;
temp(ib,:) = Rs_DIRT_monthly(ia,8);
Rs_data_veg_type2 = [Rs_data_veg_type2 temp];
var_pred_for_Rs_MAT_calculation = [var_pred_for_Rs_MAT_calculation ones(size(var_pred_for_Rs_MAT_calculation,1),1).*NaN];
[c,ia,ib] = intersect(var_pred_for_Rs_MAT_calculation(:,1),var_pred_for_Rs_MAT_calculation_DIRT(:,1));
var_pred_for_Rs_MAT_calculation(ia,end) = var_pred_for_Rs_MAT_calculation_DIRT(ib,2);
var_pred_for_Rs_MAT_calculation(1,:) = [];
%Append ChronicN data
var_pred_for_Rs_MAT_calculation = [var_pred_for_Rs_MAT_calculation ones(size(var_pred_for_Rs_MAT_calculation,1),2)*NaN];
var_pred_for_Rs_MAT_calculation(end,end-1:end) = [var_pred_ChronicN(1,1) var_pred_ChronicN(5,1)];
    
%Monthly sums
for x=1:size(Rs_data_veg_type2,1)
    nb_days = fn_nb_days_in_month(Rs_data_veg_type2(x,1),Rs_data_veg_type2(x,2));
    Rs_data_veg_type2(x,3:end) = Rs_data_veg_type2(x,3:end).*nb_days;
end

%Summer totals
Rs_control_summer_veg_type=[];
for x=1992:2009
    ind = find(Rs_data_veg_type2(:,1)==x & Rs_data_veg_type2(:,2)>=4 & Rs_data_veg_type2(:,2)<=10);
    Rs_control_summer_veg_type(x-1991,:) = [x sum(Rs_data_veg_type2(ind,3:end),1)];
end

%For SigmaPlot
%Summer Rs (all controls), average of all years
summer_controls_sums = [nanmean(Rs_control_summer_veg_type(:,2:end),1)]';

%Summer Rs (all treatments), average of all years
ind = find(SAS_data(:,4)>2);
temp = SAS_data(ind,:);
experiments = unique(temp(:,2:4),'rows');
summer_treatments_sums=[];
for x=1:length(experiments)
    ind = find(temp(:,2)==experiments(x,1) & temp(:,3)==experiments(x,2) & temp(:,4)==experiments(x,3));
    summer_treatments_sums(x,1) = mean(temp(ind,6));
end



%%%%%%%%%%%%%%%%%%%%%%
%Coefficient of variation of summer Rs (Harvard Forest reference control) -
%for interannual variability in climate and phenology
summer_Rs_ref_sums_CV = mynanstd(Rs_control_summer_veg_type(:,2:end),1)./mynanmean(Rs_control_summer_veg_type(:,2:end),1);
temp = [Rs_control_summer_veg_type [ones(17,2)*NaN ; [sum(ChronicN_Rs_monthly2(4:10,3)) sum(ChronicN_Rs_monthly2(4:10,7))]]]; %Append ChronicN (hardwood, pine)
temp(:,4)=[]; %Use for Harvard Forest Archive (file CV_interannual_and_spatial_variability.csv)
temp2 = nanmean(temp,1); %Use for Harvard Forest Archive (row 20 (mean) in file CV_interannual_and_spatial_variability.csv)
temp = mynanstd(temp(:,2:end),1)./mynanmean(temp(:,2:end),1); %Use for Harvard Forest Archive (row 21 (CV) in file CV_interannual_and_spatial_variability.csv)



%%%%%%%%%%%%%%%%%%%%%%
%Do as Bahn 2010 (Biogeosciences paper)

%Compute Rs annual sums
Rs_data_veg_type_annual = [];
y=0;
for x=1992:2009
    y=y+1;
    ind = find(Rs_data_veg_type2(:,1)==x);
    Rs_data_veg_type_annual(y,:) = [x sum(Rs_data_veg_type2(ind,3:end),1)];
end
%Add ChronicN data
Rs_data_veg_type_annual = [Rs_data_veg_type_annual ones(size(Rs_data_veg_type_annual,1),2)*NaN];
Rs_data_veg_type_annual(end,end-1:end) = [sum(ChronicN_Rs_monthly2(:,3)) sum(ChronicN_Rs_monthly2(:,7))];

%Compute DIRT Ts for controls
temp = [DIRT_Ts_gf(:,1:2) mean(DIRT_Ts_gf(:,6:11),2)];
DIRT_Ts_monthly=[];
y=0;
for x=1992:2001
    for xx=1:12
        y=y+1;
        [i,j] = fn_month_to_jday(x,xx);
        ind = find(temp(:,1)==x & temp(:,2)>=i & temp(:,2)<=j);
        DIRT_Ts_monthly(y,:) = [x xx mean(temp(ind,3))];
    end
end
[c,ia,ib] = intersect(DIRT_Ts_monthly(:,1:2),Ts_data_veg_type2(:,1:2),'rows');
temp = ones(size(Ts_data_veg_type2,1),1)*NaN;
temp(ib,:) = DIRT_Ts_monthly(ia,3);
Ts_data_veg_type2 = [Ts_data_veg_type2 temp];

%Compute ChronicN Ts for controls
ChronicN_Ts_monthly=[];
for x=1:12
    [i,j] = fn_month_to_jday(2009,x);
    ind = find(Ts_timestamps(:,1)==2009 & Ts_timestamps(:,3)>=i & Ts_timestamps(:,3)<=j);
    ChronicN_Ts_monthly(x,:) = [2009 x mean(ChronicN_Ts_pred2(ind,3)) mean(ChronicN_Ts_pred2(ind,7))];
end
[c,ia,ib] = intersect(ChronicN_Ts_monthly(:,1:2),Ts_data_veg_type2(:,1:2),'rows');
temp = ones(size(Ts_data_veg_type2,1),2)*NaN;
temp(ib,:) = ChronicN_Ts_monthly(ia,3:4);
Ts_data_veg_type2 = [Ts_data_veg_type2 temp];

%Compute Ts annual
Ts_data_veg_type_annual = [];
y=0;
for x=1992:2009
    y=y+1;
    nb_days=[];
    for xx=1:12
        nb_days(xx,1) = fn_nb_days_in_month(x,xx);
    end
    ind = find(Ts_data_veg_type2(:,1)==x);
    temp = Ts_data_veg_type2(ind,3:end).*repmat(nb_days,1,30);
    Ts_data_veg_type_annual(y,:) = [x sum(temp,1)./sum(nb_days)];
end
Ts_data_veg_type_annual(5:8,end) = NaN;

%Compute Rs at mean annual soil temperature
[Q10_controls] = fn_load_Q10_controls_annual(path);
ind = find(Q10_final_DIRT2(:,2)==2);
Q10_controls = [Q10_controls ; [Q10_final_DIRT2(ind,1) repmat([8 23 1],length(ind),1) Q10_final_DIRT2(ind,5:6)]];
Q10_controls = [Q10_controls ; [2009 9 14 1 Q10_final_ChronicN(1,5:6)]];
Q10_controls = [Q10_controls ; [2009 9 14 4 Q10_final_ChronicN(5,5:6)]];
Rs_MAT = Ts_data_veg_type_annual(:,1);
temp=unique(Q10_controls(:,2:4),'rows');
for x=1:length(temp)
    ind = find(Q10_controls(:,2)==temp(x,1) & Q10_controls(:,3)==temp(x,2) & Q10_controls(:,4)==temp(x,3));
    for xx=1:length(ind)
        ind1 = find(Ts_data_veg_type_annual(:,1)==Q10_controls(ind(xx),1));
        Rs_MAT(ind1,x+1) = [(Q10_controls(ind(xx),6).*Q10_controls(ind(xx),5).^((Ts_data_veg_type_annual(ind1,x+1)-10)./10)).*exp(0.5.*var_pred_for_Rs_MAT_calculation(ind1,x+1))]; %Apply bias correction
    end
end
Rs_MAT(Rs_MAT==0) = NaN;


%%%%
%Compute Rs at mean annual temperature using actual measurements only
Rs_MAT_collars=[];
Ts_range = 0.5;
for exp_ind=1:length(Rs_experiments_veg_type) %DIRT not used since measurements are made daily only (with soda lime)
    MAT_temp = Ts_data_veg_type_annual(:,[1 exp_ind+1]);
    [Rs_MAT_collars_temp] = fn_Rs_at_MAT_per_collar(exp_ind,data,Ts,Ts_timestamps,Rs_experiments_veg_type,MAT_temp,Ts_range);
    Rs_MAT_collars = [Rs_MAT_collars ; Rs_MAT_collars_temp];
end

%Add column with Rs annual totals
Rs_MAT_collars = [Rs_MAT_collars ones(length(Rs_MAT_collars),1)*NaN];
for x=1:length(Rs_experiments_veg_type)
    ind = find(Rs_MAT_collars(:,2)==Rs_experiments_veg_type(x,1) & Rs_MAT_collars(:,3)==Rs_experiments_veg_type(x,2) & Rs_MAT_collars(:,6)==Rs_experiments_veg_type(x,4));
    for xx=1:length(ind)
        ind2 = find(Rs_data_veg_type_annual(:,1)==Rs_MAT_collars(ind(xx),1));
        Rs_MAT_collars(ind(xx),7) = Rs_data_veg_type_annual(ind2,x+1);
    end
end
ind = find(Rs_MAT_collars(:,2)==6 & Rs_MAT_collars(:,3)==5);
temp = Rs_MAT_collars; %Use columns 1-3,5,7 for Harvard Forest Archive (columns 2-3 are equivalent to column 2 in file Rs_MAT_manual_vs_Rs_annual.csv, columns 5 and 7 are columns 3-4 in the file)


%%%%%%
%Format data
temp = [reshape(Rs_MAT(:,2:end),[],1) reshape(Rs_data_veg_type_annual(:,2:end),[],1)]; %Rs MAT, Rs annual
Rs_experiments_veg_type_temp = [Rs_experiments_veg_type(:,[1 2 4]) ; [8 23 1] ; [7 24 1] ; [7 24 4]];
temp=[];
for x=2:size(Rs_data_veg_type_annual,2)
    temp = [temp ; [Rs_MAT(:,1) repmat(Rs_experiments_veg_type_temp(x-1,:),18,1) reshape(Rs_MAT(:,x),[],1) reshape(Rs_data_veg_type_annual(:,x),[],1)]];
end
ind = find(~isnan(temp(:,end)));
temp2 = temp(ind,:); %Use for Harvard Forest Archive (columns 2-3 are equivalent to column 2 in file Rs_MAT_modeled_vs_Rs_annual.csv, column 4 is equivalent to column 3 in file, columns 5-6 are columns 4-5 in file)


%%%%%%%%%%%%%%%%%%%
%Phenology
[num,txt,raw] = xlsread('\\casfsb\biology_labs\Finzi Lab\Biogeochemistry\All Lab Data\Marc\Adrien\Harvard Forest dataset\Processed data\Bud_break_all_data.xls','Sorted');

years = ones(length(txt)-1,1);
for x=2:length(txt)
    temp = txt{x,1};
    temp = str2num(temp(1,end-3:end));
    years(x-1,1) = temp;
end
trees2 = [];
ACRU=[];
QURU=[];
TSCA=[];
trees = unique(txt(2:end,3));
for x=1:length(trees)
    ind = find(strcmp(txt(2:end,3),trees(x,1)));
    data_temp = num(ind,:);
    year_temp = years(ind,1);
    for xx=year_temp(1):year_temp(end)
        ind = find(year_temp==xx);
        data_temp2 = data_temp(ind,:);
        ind1 = find(data_temp2(:,3)>=50);
        ind1 = ind1(1);
        ind2 = find(data_temp2(:,4)>=75);
        ind2 = ind2(1);
        ind3 = find(data_temp2(:,5)>=90);
        if isempty(ind3)
            ind3 = size(data_temp2,1);
        else
            ind3 = ind3(1);
        end
        trees2 = [trees2 ; trees(x)];
        test=cell2mat(trees(x));
        if strcmp(test(1,1:4),'ACRU')
            ACRU = [ACRU ; [xx data_temp2(ind1,1) data_temp2(ind2,1) data_temp2(ind3,1)]];
        elseif strcmp(test(1,1:4),'QURU')
            QURU = [QURU ; [xx data_temp2(ind1,1) data_temp2(ind2,1) data_temp2(ind3,1)]];
        elseif strcmp(test(1,1:4),'TSCA')
            TSCA = [TSCA ; [xx data_temp2(ind1,1) data_temp2(ind2,1) data_temp2(ind3,1)]];
        end
    end
end
    
y=0;
ACRU2=[];
QURU2=[];
TSCA2=[];
for x=1992:2009
    y=y+1;
    ind = find(ACRU(:,1)==x);
    ACRU2 = [ACRU2 ; [x mean(ACRU(ind,2)) mean(ACRU(ind,3)) mean(ACRU(ind,4))]]; %Use for Harvard Forest Archive (columns 2 and 4 are columns 5-6 in file Phenology.csv)
    ind = find(QURU(:,1)==x);
    QURU2 = [QURU2 ; [x mean(QURU(ind,2)) mean(QURU(ind,3)) mean(QURU(ind,4))]]; %Use for Harvard Forest Archive (columns 2 and 4 are columns 2-3 in file Phenology.csv)
    ind = find(TSCA(:,1)==x);
    TSCA2 = [TSCA2 ; [x mean(TSCA(ind,2)) mean(TSCA(ind,3)) mean(TSCA(ind,4))]]; %Use for Harvard Forest Archive (columns 2 and 4 are columns 11-12 in file Phenology.csv)
end

hardwood = [ACRU2+QURU2]./2; %Use for Harvard Forest Archive (columns 2 and 4 are columns 8-9 in file Phenology.csv)
bud_break_hardwood_avg = mean(hardwood(:,2)); %Use for Harvard Forest Archive (file Phenology.csv)
full_leaf_hardwood_avg = mean(hardwood(:,4)); %Use for Harvard Forest Archive (file Phenology.csv)

bud_break_hemlock_avg = nanmean(TSCA2(:,2)); %Use for Harvard Forest Archive (file Phenology.csv)
full_leaf_hemlock_avg = nanmean(TSCA2(:,4)); %Use for Harvard Forest Archive (file Phenology.csv)


%%%%%%%%%%%%%%%%%%%
[num,txt,raw] = xlsread('\\casfsb\biology_labs\Finzi Lab\Biogeochemistry\All Lab Data\Marc\Adrien\Harvard Forest dataset\Processed data\Leaf_fall.xls','Sorted');

years = ones(length(txt)-1,1);
for x=2:length(txt)
    temp = txt{x,1};
    temp = str2num(temp(1,end-3:end));
    years(x-1,1) = temp;
end
trees2 = [];
ACRU_fall=[];
QURU_fall=[];
TSCA_fall=[];
trees = unique(txt(2:end,3));
for x=1:length(trees)
    ind = find(strcmp(txt(2:end,3),trees(x,1)));
    data_temp = num(ind,:);
    year_temp = years(ind,1);
    for xx=year_temp(1):year_temp(end)
        ind = find(year_temp==xx);
        data_temp2 = data_temp(ind,:);
        ind1 = find(data_temp2(:,3)>=20);
        try
            ind1 = ind1(1);
        end
        ind2 = find(data_temp2(:,4)>=80);
        try
            ind2 = ind2(1);
        end
        trees2 = [trees2 ; trees(x)];
        test=cell2mat(trees(x));
        if (isempty(ind1) & isempty(ind2))
            if strcmp(test(1,1:4),'ACRU')
                ACRU_fall = [ACRU_fall ; [xx NaN NaN]];
            elseif strcmp(test(1,1:4),'QURU')
                QURU_fall = [QURU_fall ; [xx NaN NaN]];
            elseif strcmp(test(1,1:4),'TSCA')
                TSCA_fall = [TSCA_fall ; [xx NaN NaN]];
            end
        elseif isempty(ind1)
            if strcmp(test(1,1:4),'ACRU')
                ACRU_fall = [ACRU_fall ; [xx NaN data_temp2(ind2,1)]];
            elseif strcmp(test(1,1:4),'QURU')
                QURU_fall = [QURU_fall ; [xx NaN data_temp2(ind2,1)]];
            elseif strcmp(test(1,1:4),'TSCA')
                TSCA_fall = [TSCA_fall ; [xx NaN data_temp2(ind2,1)]];
            end
        else
            if strcmp(test(1,1:4),'ACRU')
                ACRU_fall = [ACRU_fall ; [xx data_temp2(ind1,1) data_temp2(ind2,1)]];
            elseif strcmp(test(1,1:4),'QURU')
                QURU_fall = [QURU_fall ; [xx data_temp2(ind1,1) data_temp2(ind2,1)]];
            elseif strcmp(test(1,1:4),'TSCA')
                TSCA_fall = [TSCA_fall ; [xx data_temp2(ind1,1) data_temp2(ind2,1)]];
            end
        end
    end
end
    
y=0;
ACRU_fall2=[];
QURU_fall2=[];
TSCA_fall2=[];
for x=1992:2009
    y=y+1;
    ind = find(ACRU_fall(:,1)==x);
    ACRU_fall2 = [ACRU_fall2 ; [x nanmean(ACRU_fall(ind,2)) nanmean(ACRU_fall(ind,3))]]; %Use for Harvard Forest Archive (column 2 is column 7 in file Phenology.csv)
    ind = find(QURU_fall(:,1)==x);
    QURU_fall2 = [QURU_fall2 ; [x nanmean(QURU_fall(ind,2)) nanmean(QURU_fall(ind,3))]]; %Use for Harvard Forest Archive (column 2 is column 4 in file Phenology.csv)
    ind = find(TSCA_fall(:,1)==x);
    TSCA_fall2 = [TSCA_fall2 ; [x nanmean(TSCA_fall(ind,2)) nanmean(TSCA_fall(ind,3))]]; %Use for Harvard Forest Archive (column 2 is column 13 in file Phenology.csv)
end

hardwood_fall = [ACRU_fall2+QURU_fall2]./2; %Use for Harvard Forest Archive (column 2 is column 10 in file Phenology.csv)
leaf_color_hardwood_avg = nanmean(hardwood_fall(:,2)); %Use for Harvard Forest Archive (file Phenology.csv)
leaf_fall_hardwood_avg = nanmean(hardwood_fall(:,3)); %Use for Harvard Forest Archive (file Phenology.csv)

leaf_color_hemlock_avg = nanmean(TSCA_fall2(:,2)); %Use for Harvard Forest Archive (file Phenology.csv)
leaf_fall_hemlock_avg = nanmean(TSCA_fall2(:,3)); %Use for Harvard Forest Archive (file Phenology.csv)


%%%%%%%%%%%%%%%%%%%
%PPFD to determine snow cover
data = csvread('\\casfsb\biology_labs\Finzi Lab\Biogeochemistry\All Lab Data\Marc\Adrien\Harvard Forest dataset\Processed data\EMS_radiation.csv',1,0);

pardown = data(:,13);
parup = data(:,12);
pardown(1:85587) = pardown(1:85587).*1000;
parup(1:85587) = parup(1:85587).*1000;
dayind = find(pardown>8);
pardown = pardown(dayind,1);
parup = parup(dayind,1);
year = data(dayind,1);
day = data(dayind,3);

y=0;
for x=year(1):year(end)
    for xx=1:365
        y=y+1;
        ind = find(year==x & day==xx);
        parup2(y,1) = nanmean(parup(ind,1)); %Use in Harvard Forest Archive (file PPFD_snow.csv)
        pardown2(y,1) = nanmean(pardown(ind,1)); %Use in Harvard Forest Archive (file PPFD_snow.csv)
        day2(y,1) = xx;
        year2(y,1) = x;
    end
end
for x=1:365
    ind = find(day2==x);
    parup3(x,1) = nanmean(parup2(ind));
    pardown3(x,1) = nanmean(pardown2(ind));
end

ratio = parup2./pardown2; %Use in Harvard Forest Archive (file PPFD_snow.csv)


%%%%%%%%%%%%%%%%%%%%%%%%%
%Wind speed and u*
%EMS Tower
data_filled = csvread([path,'\EMS - filled data2.csv'],2,0);

EMS_year = data_filled(:,1);
EMS_jday = data_filled(:,3);
EMS_hour = data_filled(:,4);
ustar = data_filled(:,10);

data_final = csvread([path,'\EMS - final data2.csv'],1,0);

windspeed = data_final(:,5);
EMS_timestamps = [EMS_year EMS_jday EMS_hour];

wind=[];
ustar2=[];
for x=1:366
    ind = find(EMS_timestamps(:,2)==x);
    wind(x,1) = nanmean(windspeed(ind,1)); %Use for Harvard Forest Archive (file Wind_speed_friction_velocity.csv)
    ustar2(x,1) = nanmean(ustar(ind,1))./100; %Use for Harvard Forest Archive (file Wind_speed_friction_velocity.csv)
end