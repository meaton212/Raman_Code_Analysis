clearvars; 
close all; 

%% Input data:
% The program expects multiple (up to 10) files in .txt or .dat format. In each file, the Raman shift
% should be included in the first column, followed by the Intensity data in
% consecutive columns (as many spectra available, and not neccesarily the same number of spectra per file).

% The location and names of the files should be indicated in path and file_name
% variables, respectively. The total number of files should be included in
% total.

path = '/Users/natms/Documents/Trabajo/papers:manuscripts/2023 Raman analysis/code/';
file_name={ 'p-G 2-2';
    'CuMINT';};
total=2; % Total number of files
type = '.txt';
name={'6,5-SWCNT';
    'CuMINT'}; % Indicate here the names of the data for legend, titles and data output. No spaces allowed.


%% Model: Choose the method to analyse spectral features 

rbm=0;% Set to 1 if RBM analysis desired
lorentz=0; % lorentz fits for G, D and 2D (=0 for simple identification; =1 for lorentzian peak fitting)
nt=0; % Splitting G band in G+ and G-? (=0 for no; =1 for yes), note only valid if lorentz=1
map=1; % Set to 1 if spatial mapping desired

%% Normalization (mandatory): 
% Choose normalization spectral range. The spectra will be normalised to
% the maximun within the specified spectral range: 
normLow=1500; %Lower limit in cm-1
normHigh=1700; %Upper limit in cm-1

%% Peak identification: Intensity and Shifts (mandatory)
% Select appropriate spectral range where the full peaks are resolved. 

% Peak 1: G band * Note that peak 1 can be optionally fitted into 2 Lorentzian
band1Low=1470; %Lower limit in cm-1
band1High=1650; %Upper limit in cm-1
% Peak 2: D band
band2Low=1230; %Lower limit in cm-1
band2High=1360; %Upper limit in cm-1
% Peak 3: 2D band
band3Low=2450; %Lower limit in cm-1
band3High=2800; %Upper limit in cm-1

% RBM modes (only if RBM=1)
RBMregion_Low=100; % Lower limit in cm-1
RBMregion_High=500; % Upper limit in cm-1
Prom=0.01; % This value sets the max limit at which peaks will be considered. 
%Local maxima with a maximum height lower than Prom will be discarded.

%% Fitting parameters (only if lorentz=1)

% If lorentzian fits are selected, the user is required to indicate
% initial guess for the peak position (cm-1), FWHM (cm-1) and intensity

%Peak 1: G
Gmin_init=[1530 25 0.3; 1530 20 0.20]; %[center FWHM intensity_max] for G- (Only if nt=1)
Gplus_init=[1578 30 1; 1580 30 1]; %[center FWHM intensity_max] for G+ (or G if nt=0)
%Peak 2: D
D_init=[1300 40 0.08; 1305 40 0.06]; %[center FWHM intensity_max] for D band
%Peak 3: 2D
init_2D=[2600 50 0.4; 2640 50 0.3]; %[center FWHM intensity_max] for 2D band

%% Mapping options (only if map=1)
% Select here if mapping is wanted, and indicate the details of the map

col=32 ; % number of columns / pixels in x
raws=32 ; % number of raws / pixels in y

colmum=32.5; % micrometers of the map in the x axis
rawsmum=24; % micrometers of the map in the y axis

%% Ploting options/Output plots:
% Choose the desired output plots (1 is yes and 0 is no). 
raw=1; % Make figure af all raw spectra.
norm=1; % Plot all normalised spectra.
range=1; % Plot spectral regions chosen for G, D and 2D bands in intensity calculation
peaks=1; % plot found peaks in the RBM region, Note that if the
% number of spectra is very high, the computing is going to slow down
% significantly.
correlations=1; % Plot correlations between spectral features
width=3; % width of Raman shift histograms 
width_fw=5; % width of FWHM histograms 


%% Code Starts here:

if lorentz==0
    nt=0;
end

for z=1:total

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Import data %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IN = importdata([path file_name{z} type]);
Shift(:) = IN(:,1); % Raman shift is the first column of the file
Intensity(:,:)=IN(:,2:end); % Raman intensities are all data from column 2 til end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%% Set color scales for plots %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmap = {gray(length(Intensity(1,:))+20);
winter(length(Intensity(1,:))+5);
summer(length(Intensity(1,:))+5);
autumn(length(Intensity(1,:))+5);
parula(length(Intensity(1,:))+5);
cool(length(Intensity(1,:))+5);
spring(length(Intensity(1,:))+5);
copper(length(Intensity(1,:))+5);
pink(length(Intensity(1,:))+5);
hot(length(Intensity(1,:))+5)}; % Sets color maps for different files

color=cmap{z};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%% Set distribution of subfigures %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if total==1
    div=[1 1];
elseif total==2
    div=[1 2];
elseif total==3
    div=[1 3];
elseif total==4
    div=[2 2];
elseif total==5
    div=[2 3];
elseif total==6
    div=[2 3];
elseif total==10
    div=[2 5];
else total==7
    div=[3 3];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Normalization %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, indLow]=min(abs(Shift-normLow));
[~, indHigh]=min(abs(Shift-normHigh));

Intensity_norm=zeros(size(Intensity));

for n=1:length(Intensity(1,:))
Intensity_norm(:,n)=(Intensity(:,n)-min(Intensity(:,n))); % substract min
Intensity_norm(:,n)=Intensity_norm(:,n)./max(Intensity_norm(indLow:indHigh,n)); % divide by norm
end

Intensity_av=mean(Intensity_norm'); % Calculate average spectra

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Spectral features calculation %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define spectral range for G, D and 2D
[~ ,ind1Low]=min(abs(Shift-band1Low));
[~ ,ind1High]=min(abs(Shift-band1High));

[~, ind2Low]=min(abs(Shift-band2Low));
[~, ind2High]=min(abs(Shift-band2High));

[~, ind3Low]=min(abs(Shift-band3Low));
[~, ind3High]=min(abs(Shift-band3High));

% Define Intensity and shift vectors

Intensity_1range=zeros(ind1High-ind1Low+1,length(Intensity_norm(1,:)));
Shift_range1=Shift(ind1Low:ind1High);
Intensity_2range=zeros(ind2High-ind2Low+1,length(Intensity_norm(1,:)));
Shift_range2=Shift(ind2Low:ind2High);
Intensity_3range=zeros(ind3High-ind3Low+1,length(Intensity_norm(1,:)));
Shift_range3=Shift(ind3Low:ind3High);
Int_G=zeros(1,length(Intensity_norm(1,:)));
Int_D=zeros(1,length(Intensity_norm(1,:)));
Int_2D=zeros(1,length(Intensity_norm(1,:)));
center_G=zeros(1,length(Intensity_norm(1,:)));
center_D=zeros(1,length(Intensity_norm(1,:)));
center_2D=zeros(1,length(Intensity_norm(1,:)));


for n=1:length(Intensity_norm(1,:))   
   
   Intensity_1range(:,n)=Intensity_norm(ind1Low:ind1High,n);
   Intensity_2range(:,n)=Intensity_norm(ind2Low:ind2High,n);
   Intensity_3range(:,n)=Intensity_norm(ind3Low:ind3High,n);
    
   %%
   %%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%% Method 1 %%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%
   if lorentz==0;
   % peak 1/G: 
   Int_G(n)=max(Intensity_1range(:,n))-min(Intensity_1range(:,n)); % Calculate intensity
   [~, indMax1]=max(Intensity_1range(:,n));
   center_G(n)=Shift_range1(indMax1); % Calculate position
   
   % peak 2/D: 
   Int_D(n)=max(Intensity_2range(:,n))-min(Intensity_2range(:,n));
   [~, indMax2]=max(Intensity_2range(:,n));
   center_D(n)=Shift_range2(indMax2);
   
   % peak 3/2D: 
   Int_2D(n)=max(Intensity_3range(:,n))-min(Intensity_3range(:,n));
   [~, indMax3]=max(Intensity_3range(:,n));
   center_2D(n)=Shift_range3(indMax3);
   
   end
  %% 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%% Method 2: Lorentzian fitting %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
   if lorentz==1;
    
    fit_func=@(gamma,x)(gamma(1)./((x-gamma(2)).^2+gamma(3))+gamma(4));% Single lorentz fitting function
    
    % Peak 2: D band
    InitGuess_D=[D_init(z,3)*(D_init(z,2)/2)^2 D_init(z,1) (D_init(z,2)/2)^2 0];
    Intensity_2range(:,n)=Intensity_norm(ind2Low:ind2High,n);
    gamma01=InitGuess_D(:);
    gamma_D = zeros(length(gamma01), 1);
    gamma_D(:) = nlinfit(Shift_range2',Intensity_2range(:,n),fit_func,gamma01);  
    
    I_D_fit(:,n)=(gamma_D(1)./((Shift_range2(:)'-gamma_D(2)).^2+gamma_D(3))+gamma_D(4));
    center_D(n)=gamma_D(2);
    FWHM_D(n)=2*sqrt(gamma_D(3));
    Int_D(n)=gamma_D(1)/gamma_D(3);  
       
   % Peak 3: 2D band
    InitGuess_2D=[init_2D(z,3)*(init_2D(z,2)/2)^2 init_2D(z,1) (init_2D(z,2)/2)^2 0];
    Intensity_3range(:,n)=Intensity_norm(ind3Low:ind3High,n);
    gamma01=InitGuess_2D(:);
    gamma_2D = zeros(length(gamma01), 1);
    gamma_2D(:) = nlinfit(Shift_range3',Intensity_3range(:,n),fit_func,gamma01);  % 104:end
    
    I_2D_fit(:,n)=(gamma_2D(1)./((Shift_range3(:)'-gamma_2D(2)).^2+gamma_2D(3))+gamma_2D(4));
    center_2D(n)=gamma_2D(2);
    FWHM_2D(n)=2*sqrt(gamma_2D(3));
    Int_2D(n)=gamma_2D(1)/gamma_2D(3);
    
    if nt==0
    % Peak 1: G band, single peak fitting
    InitGuess_G=[Gplus_init(z,3)*(Gplus_init(z,2)/2)^2 Gplus_init(z,1) (Gplus_init(z,2)/2)^2 0];    
    Intensity_1range(:,n)=Intensity_norm(ind1Low:ind1High,n);
    gamma01=InitGuess_G(:);
    gamma_G = zeros(length(gamma01), 1);
    gamma_G(:) = nlinfit(Shift_range1',Intensity_1range(:,n),fit_func,gamma01);  
    
    I_G_fit(:,n)=(gamma_G(1)./((Shift_range1(:)'-gamma_G(2)).^2+gamma_G(3))+gamma_G(4));
    center_G(n)=gamma_G(2);
    FWHM_G(n)=2*sqrt(gamma_G(3));
    Int_G(n)=gamma_G(1)/gamma_G(3); 
    
    else
    % Peak 1: G band split G+ and G-, double peak fitting
    InitGuess_G=[Gmin_init(z,3)*(Gmin_init(z,2)/2)^2 Gmin_init(z,1) (Gmin_init(z,2)/2)^2 Gplus_init(z,3)*(Gplus_init(z,2)/2)^2 Gplus_init(z,1) (Gplus_init(z,2)/2)^2 0.1];
    Intensity_1range(:,n)=Intensity_norm(ind1Low:ind1High,n);
    fit_func=@(gamma,x)(gamma(1)./((x-gamma(2)).^2+gamma(3))+gamma(4)./((x-gamma(5)).^2+gamma(6))+gamma(7));% lorentz
    gamma01=InitGuess_G(:);
    gamma_G = zeros(length(gamma01), 1);
    gamma_G(:) = nlinfit(Shift_range1',Intensity_1range(:,n),fit_func,gamma01);  % 104:end
    
    I_G_fit(n,:)=(gamma_G(1)./((Shift_range1(:)'-gamma_G(2)).^2+gamma_G(3))+gamma_G(4)./((Shift_range1(:)'-gamma_G(5)).^2+gamma_G(6))+gamma_G(7));
 
        if gamma_G(2)<gamma_G(5)
        center_Gmin(n)=gamma_G(2);
        FWHM_Gmin(n)=2*sqrt(gamma_G(3));
        Int_Gmin(n)=gamma_G(1)/gamma_G(3);
        center_Gplus(n)=gamma_G(5);
        FWHM_Gplus(n)=2*sqrt(gamma_G(6));
        Int_Gplus(n)=gamma_G(4)/gamma_G(6);
        else
        center_Gplus(n)=gamma_G(2);
        FWHM_Gplus(n)=2*sqrt(gamma_G(3));
        Int_Gplus(n)=gamma_G(1)/gamma_G(3);
        center_Gmin(n)=gamma_G(5);
        FWHM_Gmin(n)=2*sqrt(gamma_G(6));
        Int_Gmin(n)=gamma_G(4)/gamma_G(6);
        end
 
    end
   end
end

%% Calculate mean values and standard deviation

if nt==1
I21=Int_D./Int_Gplus; % Calculate intensity ratio
I21_av=mean(I21);
I21_error=std(I21);

I2dg=Int_2D./Int_Gplus; % Calculate intensity ratio
I2dg_av=mean(I2dg);
I2dg_error=std(I2dg);

pos_Gplus_av=mean(center_Gplus);
pos_Gplus_err=std(center_Gplus);
pos_Gmin_av=mean(center_Gmin);
pos_Gmin_err=std(center_Gmin);
pos_2D_av=mean(center_2D);
pos_2D_err=std(center_2D);
pos_D_av=mean(center_D);
pos_D_err=std(center_D);

else
I21=Int_D./Int_G; % Calculate intensity ratio
I21_av=mean(I21);
I21_error=std(I21);

I2dg=Int_2D./Int_G; % Calculate intensity ratio
I2dg_av=mean(I2dg);
I2dg_error=std(I2dg);

pos_G_av=mean(center_G);
pos_G_err=std(center_G);
pos_2D_av=mean(center_2D);
pos_2D_err=std(center_2D);
pos_D_av=mean(center_D);
pos_D_err=std(center_D);
end

if lorentz==1 & nt==0
FWHM_G_av=mean(FWHM_G);
FWHM_G_err=std(FWHM_G);
FWHM_2D_av=mean(FWHM_2D);
FWHM_2D_err=std(FWHM_2D);
FWHM_D_av=mean(FWHM_D);
FWHM_D_err=std(FWHM_D);   
end

if lorentz==1 & nt==1
FWHM_Gplus_av=mean(FWHM_Gplus);
FWHM_Gplus_err=std(FWHM_Gplus);   
FWHM_Gmin_av=mean(FWHM_Gmin);
FWHM_Gmin_err=std(FWHM_Gmin); 
FWHM_2D_av=mean(FWHM_2D);
FWHM_2D_err=std(FWHM_2D);
FWHM_D_av=mean(FWHM_D);
FWHM_D_err=std(FWHM_D);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% RBM modes Shifts  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if rbm
[~, indRBMLow]=min(abs(Shift-RBMregion_Low));
[~, indRBMHigh]=min(abs(Shift-RBMregion_High));
PeaksInt=[];
PeaksLoc=[];
IntRBM={};
LocRBM={};

for n=1:length(Intensity_norm(1,:))

if peaks
figure(401),
subplot(div(1),div(2),z);hold on
findpeaks(Intensity_norm(indRBMLow:indRBMHigh,n),Shift(indRBMLow:indRBMHigh),'MinPeakProminence',Prom)
xlabel('Raman shift / cm^{-1}','FontSize',12);
ylabel('Intensity/a.u.','FontSize',12);
title([name{z} ' :RBM region peaks']);
set(gca,'Box','on');
%plotbrowser('on');
end

[pksRBM,locsRBM]=findpeaks(Intensity_norm(indRBMLow:indRBMHigh,n),Shift(indRBMLow:indRBMHigh),'MinPeakProminence',Prom);
PeaksInt=[PeaksInt [pksRBM]'];
PeaksLoc=[PeaksLoc [locsRBM]];

IntRBM(n,:)={[pksRBM]};
LocRBM(n,:)={[locsRBM]};
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Output data %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save data in txt file

if rbm
    
if lorentz==0 & nt==1
T=table(Int_Gmin',center_Gmin',FWHM_Gmin',Int_Gplus',center_Gplus',FWHM_Gplus', Int_D',center_D', Int_2D',center_2D',IntRBM,LocRBM,'VariableNames',{'Intensity_G-','Shift_G-','FWHM_G-','Intensity_G+','Shift_G+','FWHM_G+','Intensity_D','Shift_D','Intensity_2D','Shift_2D','Intensity_RBM','Shift_RBM'});
end

if lorentz==0 & nt==0
T=table(Int_G',center_G', Int_D',center_D', Int_2D',center_2D',IntRBM,LocRBM,'VariableNames',{'Intensity_G','Shift_G','Intensity_D','Shift_D','Intensity_2D','Shift_2D','Intensity_RBM','Shift_RBM'});
end

if lorentz==1 & nt==1
T=table(Int_Gmin',center_Gmin',FWHM_Gmin',Int_Gplus',center_Gplus',FWHM_Gplus', Int_D',center_D',FWHM_D', Int_2D',center_2D',FWHM_2D',IntRBM,LocRBM,'VariableNames',{'Intensity_G-','Shift_G-','FWHM_G-','Intensity_G+','Shift_G+','FWHM_G+','Intensity_D','Shift_D','FWHM_D','Intensity_2D','Shift_2D','FWHM_2D','Intensity_RBM','Shift_RBM'});
end

if lorentz==1 & nt==0
T=table(Int_G',center_G',FWHM_G', Int_D',center_D',FWHM_D', Int_2D',center_2D',FWHM_2D',IntRBM,LocRBM,'VariableNames',{'Intensity_G','Shift_G','FWHM_G','Intensity_D','Shift_D','FWHM_D','Intensity_2D','Shift_2D','FWHM_2D','Intensity_RBM','Shift_RBM'});
end

else

if lorentz==0 & nt==1
T=table(Int_Gmin',center_Gmin',FWHM_Gmin',Int_Gplus',center_Gplus',FWHM_Gplus', Int_D',center_D', Int_2D',center_2D','VariableNames',{'Intensity_G-','Shift_G-','FWHM_G-','Intensity_G+','Shift_G+','FWHM_G+','Intensity_D','Shift_D','Intensity_2D','Shift_2D'});
end

if lorentz==0 & nt==0
T=table(Int_G',center_G', Int_D',center_D', Int_2D',center_2D','VariableNames',{'Intensity_G','Shift_G','Intensity_D','Shift_D','Intensity_2D','Shift_2D'});
end

if lorentz==1 & nt==1
T=table(Int_Gmin',center_Gmin',FWHM_Gmin',Int_Gplus',center_Gplus',FWHM_Gplus', Int_D',center_D',FWHM_D', Int_2D',center_2D',FWHM_2D','VariableNames',{'Intensity_G-','Shift_G-','FWHM_G-','Intensity_G+','Shift_G+','FWHM_G+','Intensity_D','Shift_D','FWHM_D','Intensity_2D','Shift_2D','FWHM_2D'});
end

if lorentz==1 & nt==0
T=table(Int_G',center_G',FWHM_G', Int_D',center_D',FWHM_D', Int_2D',center_2D',FWHM_2D','VariableNames',{'Intensity_G','Shift_G','FWHM_G','Intensity_D','Shift_D','FWHM_D','Intensity_2D','Shift_2D','FWHM_2D'});
end
    
end

nameT=strcat(name(z),'_results','.csv');
namefinal=nameT{1};
writetable(T,namefinal);%,'Delimiter',' ');

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Figures %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

%% General
% Raw spectra
if raw
   figure(1), hold on;
   subplot(div(1),div(2),z);
   %colororder(color);
   set(gca, 'ColorOrder', color, 'NextPlot', 'replacechildren');
   plot(Shift,Intensity(:,:));
   xlabel('Raman shift / cm^{-1}','FontSize',12);
   ylabel('Intensity / a.u.','FontSize',12);
   title([name{z} ' :Raw data']);
   %plotbrowser('on');
   set(gca,'Box','on');
   %legend show
end

% Normalised spectra
if norm
figure(2), hold on;
subplot(div(1),div(2),z);
set(gca, 'ColorOrder', color, 'NextPlot', 'replacechildren');
plot(Shift,Intensity_norm(:,1:end));
xlabel('Raman shift / cm^{-1}','FontSize',12);
ylabel('Intensity / a.u.','FontSize',12);
title([name{z} ' :Normalised spectra']);
set(gca,'Box','on');
%plotbrowser('on');
%legend show
end

% Average spectra
figure(3), hold on;
plot(Shift,Intensity_av(),'color',color(round(end/2),:),'DisplayName',name{z});
xlabel('Raman shift / cm^{-1}','FontSize',12);
ylabel('Intensity / a.u.','FontSize',12);
title('Average spectra');
set(gca,'Box','on');
%plotbrowser('on');
legend show

% Spectral range for Intensity ratio
if range
   figure(10),
   subplot(div(1),div(2),z);hold on
   set(gca, 'ColorOrder', color, 'NextPlot', 'replacechildren');
   plot(Shift(ind1Low:ind1High),Intensity_1range(:,:));hold on
   plot(Shift(ind2Low:ind2High),Intensity_2range(:,:));hold on
   plot(Shift(ind3Low:ind3High),Intensity_3range(:,:));
   xlabel('Raman shift / cm^{-1}','FontSize',12);
   ylabel('Intensity / a.u.','FontSize',12);
   title([name{z} ': Range used for Intensity ratio calculation'])
   set(gca,'Box','on');
end

%% Lorentz fit results: plots the spectra together with the best fit

if lorentz
    
figure(101);
    subplot(div(1),div(2),z); hold on
    set(gca, 'ColorOrder', color, 'NextPlot', 'replacechildren');
    plot(Shift_range1(:),Intensity_1range(:,:));hold on,
    plot(Shift_range1(:),I_G_fit(:,:),'r');
    xlabel('Raman shift / cm^{-1}','FontSize',12)
    ylabel('Intensity / a.u.','FontSize',12)
    set(gca,'Box','on');     
    title([name(z) 'Fitting results G'],'FontSize',12);
    

figure(201);
    subplot(div(1),div(2),z); hold on
    set(gca, 'ColorOrder', color, 'NextPlot', 'replacechildren');
    plot(Shift_range2(:),Intensity_2range(:,:));hold on,
    plot(Shift_range2(:),I_D_fit(:,:),'r');
    xlabel('Raman shift / cm^{-1}','FontSize',12)
    ylabel('Intensity / a.u.','FontSize',12)
    set(gca,'Box','on'); 
    title([name(z) 'Fitting results D'],'FontSize',12);
    
figure(301);
    subplot(div(1),div(2),z); hold on
    set(gca, 'ColorOrder', color, 'NextPlot', 'replacechildren');
    plot(Shift_range3(:),Intensity_3range(:,:));hold on,
    plot(Shift_range3(:),I_2D_fit(:,:),'r');
    xlabel('Raman shift / cm^{-1}','FontSize',12)
    ylabel('Intensity / a.u.','FontSize',12)
    set(gca,'Box','on'); 
    title([name(z) 'Fitting results 2D'],'FontSize',12);

end

%% Histograms

% Histogram intensity ratio
binsInt=round((max(I21)-min(I21))/0.01);
figure(11), hold on;
histogram(I21,binsInt,'FaceColor',color(round(end/2),:),'DisplayName',[name{z} ': {I_d/I_g}=' num2str(I21_av) '\pm' num2str(I21_error)]);
xlabel('I_d/I_g','FontSize',12);
ylabel('counts','FontSize',12);
title('Intensity ratio: {I_d/I_g}','FontSize',12);
set(gca,'Box','on');
legend show


% Raman shift G mode

if nt
figure(100), hold on;
g1=histogram(center_Gplus,'FaceColor',color(round(end/3),:),'DisplayName',[name{z} 'shift G^+=' num2str(pos_Gplus_av) '\pm' num2str(pos_Gplus_err)]);
g1.BinWidth = width;  
g2=histogram(center_Gmin,'FaceColor',color(round(end/2),:),'DisplayName',[name{z} 'shift G^-=' num2str(pos_Gmin_av) '\pm' num2str(pos_Gmin_err)]);
g2.BinWidth = width; 
xlabel('Raman shift / cm^{-1}','FontSize',12);
ylabel('counts','FontSize',12);
title('Raman shift G modes','FontSize',12);
set(gca,'Box','on');
xlim([band1Low-10 band1High+10])
legend show

figure(99), hold on;
Iplus_minus=Int_Gplus./Int_Gmin;
binsInt=round((max(Iplus_minus)-min(Iplus_minus))/0.1);
g3=histogram(Iplus_minus,binsInt,'FaceColor',color(round(end/2),:),'DisplayName',name{z}) ;
xlabel('I_{G^+}/I_{G^-}','FontSize',12);
ylabel('counts','FontSize',12);
title('Intensity ratio G modes','FontSize',12);
set(gca,'Box','on');
legend show

figure(98), hold on;
g4=histogram(FWHM_Gmin,'FaceColor',color(round(end/2),:),'DisplayName',[name{z} ' FWHM G^-=' num2str(FWHM_Gmin_av) '\pm' num2str(FWHM_Gmin_err)]);
g4.BinWidth = width_fw;
g5=histogram(FWHM_Gplus,'FaceColor',color(round(end/3),:),'DisplayName',[name{z} ' FWHM G^+=' num2str(FWHM_Gplus_av) '\pm' num2str(FWHM_Gplus_err)]);
g5.BinWidth = width_fw;
xlabel('FWHM / cm^{-1}','FontSize',12);
ylabel('counts','FontSize',12);
title('FWHM G band','FontSize',12);
set(gca,'Box','on');
legend show

else
figure(100), hold on
h1=histogram(center_G,'FaceColor',color(round(end/2),:),'DisplayName',[name{z} '; Shift G=' num2str(pos_G_av) '\pm' num2str(pos_G_err)]);
h1.BinWidth = width;
xlabel('Raman shift / cm^{-1}','FontSize',12);
ylabel('counts','FontSize',12);
title('Raman shift G band','FontSize',12);
legend show;
set(gca,'Box','on');

if lorentz
figure(98), hold on;
g4=histogram(FWHM_G,'FaceColor',color(round(end/2),:),'DisplayName',[name{z} ' FWHM G=' num2str(FWHM_G_av) '\pm' num2str(FWHM_G_err)]);
g4.BinWidth = width_fw;
xlabel('FWHM / cm^{-1}','FontSize',12);
ylabel('counts','FontSize',12);
title('FWHM G band','FontSize',12);
set(gca,'Box','on');
legend show
end

end

% Raman shift D mode
binsShift=round(max(center_D)-min(center_D));
figure(200), hold on;
h2=histogram(center_D,'FaceColor',color(round(end/2),:),'DisplayName',[name{z} '; Shift D=' num2str(pos_D_av) '\pm' num2str(pos_D_err)]);
h2.BinWidth = width;
xlabel('Raman shift / cm^{-1}','FontSize',12);
ylabel('counts','FontSize',12);
title('Raman shift D band','FontSize',12);
set(gca,'Box','on');
legend show

if lorentz
figure(199), hold on;
h3=histogram(FWHM_D,'FaceColor',color(round(end/2),:),'DisplayName',[name{z} '; FWHM D=' num2str(FWHM_D_av) '\pm' num2str(FWHM_D_err)]);
h3.BinWidth = width_fw;
xlabel('FWHM / cm^{-1}','FontSize',12);
ylabel('counts','FontSize',12);
title('FWHM D band','FontSize',12);
set(gca,'Box','on');
legend show
end

% Raman shift 2D mode
binsShift=round(max(center_2D)-min(center_2D));
figure(300), hold on;
h4=histogram(center_2D,'FaceColor',color(round(end/2),:),'DisplayName',[name{z} '; Shift 2D=' num2str(pos_2D_av) '\pm' num2str(pos_2D_err)]);
h4.BinWidth = width;
xlabel('Raman shift / cm^{-1}','FontSize',12);
ylabel('counts','FontSize',12);
title('Raman shift 2D band','FontSize',12);
set(gca,'Box','on');
legend show

if lorentz
figure(299), hold on;
h5=histogram(FWHM_2D,'FaceColor',color(round(end/2),:),'DisplayName',[name{z} '; FWHM 2D=' num2str(FWHM_2D_av) '\pm' num2str(FWHM_2D_err)]);
h5.BinWidth = width_fw;
xlabel('FWHM / cm^{-1}','FontSize',12);
ylabel('counts','FontSize',12);
title('FWHM 2D band','FontSize',12);
set(gca,'Box','on');
legend show
end

% Raman shift RBM
if rbm
binsShift=round((max(PeaksLoc)-min(PeaksLoc))/2);
figure(400), hold on
h4=histogram(PeaksLoc,'FaceColor',color(round(end/2),:),'DisplayName',name{z});
h4.BinWidth = width;
xlabel('Raman shift / cm^{-1}','FontSize',12);
ylabel('counts','FontSize',12);
title('Raman shift RBM modes','FontSize',12);
set(gca,'Box','on');
legend show
end

%% Correlations

if correlations

% Position Band 1 vs I21
if nt  
p1=fit(I21',center_Gplus','poly1');
p1_coeff = coeffvalues(p1);
p1_confint = confint(p1);

figure(1000), hold on;
plot(I21(:),center_Gplus,'o','color',color(round(end/2),:),'MarkerFaceColor',color(round(end/2),:),'DisplayName',name{z});
plot(linspace(min(I21),max(I21)),polyval(p1_coeff,linspace(min(I21),max(I21))),'-','DisplayName',['Linear fit, m=' num2str(p1_coeff(1)) '\pm' num2str((p1_confint(2,1) - p1_confint(1,1))/2)],'color',color(round(end/2),:));
xlabel('I_d/I_g','FontSize',12);
ylabel('Raman shift G^+ band / cm^{-1}','FontSize',12);
title('G^+ band vs I_d/I_g','FontSize',12);
set(gca,'Box','on');
legend show;   

else
    
p1=fit(I21',center_G','poly1');
p1_coeff = coeffvalues(p1);
p1_confint = confint(p1);

figure(1000), hold on;
plot(I21(:),center_G,'o','color',color(round(end/2),:),'MarkerFaceColor',color(round(end/2),:),'DisplayName',name{z});
plot(linspace(min(I21),max(I21)),polyval(p1_coeff,linspace(min(I21),max(I21))),'-','DisplayName',['Linear fit, m=' num2str(p1_coeff(1)) '\pm' num2str((p1_confint(2,1) - p1_confint(1,1))/2)],'color',color(round(end/2),:));
xlabel('I_d/I_g','FontSize',12);
ylabel('Raman shift G band / cm^{-1}','FontSize',12);
title('G band vs I_d/I_g','FontSize',12);
set(gca,'Box','on');
legend show;
end

% Position Band 2 vs I21
p2=fit(I21',center_D','poly1');
p2_coeff = coeffvalues(p2);
p2_confint = confint(p2);

figure(1001), hold on;
plot(I21(:),center_D(:),'o','DisplayName',name{z},'color',color(round(end/2),:),'MarkerFaceColor',color(round(end/2),:))
plot(linspace(min(I21),max(I21)),polyval(p2_coeff,linspace(min(I21),max(I21))),'-','DisplayName',['Linear fit, m=' num2str(p2_coeff(1)) '\pm' num2str((p2_confint(2,1) - p2_confint(1,1))/2)],'color',color(round(end/2),:));
xlabel('I_d/I_g','FontSize',12);
ylabel('Raman shift D band / cm^{-1}','FontSize',12);
title('D band vs I_d/I_g','FontSize',12);
set(gca,'Box','on');
legend show;

% Position Band 3 vs I21
p3=fit(I21',center_2D','poly1');
p3_coeff = coeffvalues(p3);
p3_confint = confint(p3);

figure(1002), hold on;
plot(I21(:),center_2D(:),'o','DisplayName',name{z},'color',color(round(end/2),:),'MarkerFaceColor',color(round(end/2),:))
plot(linspace(min(I21),max(I21)),polyval(p3_coeff,linspace(min(I21),max(I21))),'-','DisplayName',['Linear fit, m=' num2str(p3_coeff(1)) '\pm' num2str((p3_confint(2,1) - p3_confint(1,1))/2)],'color',color(round(end/2),:));
xlabel('I_d/I_g','FontSize',12);
ylabel('Raman shift 2D band / cm^{-1}','FontSize',12);
title('2D band vs I_d/I_g','FontSize',12);
set(gca,'Box','on');
legend show;

% Position Band 3 vs position band 1
if nt
p4=fit(center_Gplus',center_2D','poly1');
p4_coeff = coeffvalues(p4);
p4_confint = confint(p4);

figure(1003), hold on;
plot(center_Gplus(:),center_2D(:),'o','DisplayName',name{z},'color',color(round(end/2),:),'MarkerFaceColor',color(round(end/2),:))
plot(linspace(min(center_Gplus),max(center_Gplus)),polyval(p4_coeff,linspace(min(center_Gplus),max(center_Gplus))),'-','DisplayName',['Linear fit, m=' num2str(p4_coeff(1)) '\pm' num2str((p4_confint(2,1) - p4_confint(1,1))/2)],'color',color(round(end/2),:));
xlabel('Raman Shift G^+ band / cm^{-1}','FontSize',12);
ylabel('Raman shift 2D band / cm^{-1}','FontSize',12);
title('2D vs G^+','FontSize',12);
set(gca,'Box','on');
legend show;
    
else
p4=fit(center_G',center_2D','poly1');
p4_coeff = coeffvalues(p4);
p4_confint = confint(p4);

figure(1003), hold on;
plot(center_G(:),center_2D(:),'o','DisplayName',name{z},'color',color(round(end/2),:),'MarkerFaceColor',color(round(end/2),:))
plot(linspace(min(center_G),max(center_G)),polyval(p4_coeff,linspace(min(center_G),max(center_G))),'-','DisplayName',['Linear fit, m=' num2str(p4_coeff(1)) '\pm' num2str((p4_confint(2,1) - p4_confint(1,1))/2)],'color',color(round(end/2),:));
xlabel('Raman Shift G band / cm^{-1}','FontSize',12);
ylabel('Raman shift 2D band / cm^{-1}','FontSize',12);
title('2D vs G','FontSize',12);
set(gca,'Box','on');
legend show;
end

if lorentz

% Position band 2 vs FWHM band 2
figure(1005), hold on;
plot(center_D(:),FWHM_D(:),'o','DisplayName',name{z},'color',color(round(end/2),:),'MarkerFaceColor',color(round(end/2),:))
xlabel('Raman Shift D band / cm^{-1}','FontSize',12);
ylabel('FWHM D band / cm^{-1}','FontSize',12);
title('FWHM vs Shift D band','FontSize',12);
set(gca,'Box','on');
legend show;

% Position band 3 vs FWHM band 3
figure(1006), hold on;
plot(center_2D(:),FWHM_2D(:),'o','DisplayName',name{z},'color',color(round(end/2),:),'MarkerFaceColor',color(round(end/2),:))
xlabel('Raman Shift 2D band / cm^{-1}','FontSize',12);
ylabel('FWHM 2D band / cm^{-1}','FontSize',12);
title('FWHM vs Shift 2D band','FontSize',12);
set(gca,'Box','on');
legend show;

if nt
% Position band 1 vs FWHM band 1
figure(1004), hold on;
plot(center_Gplus(:),FWHM_Gplus(:),'o','DisplayName',name{z},'color',color(round(end/2),:),'MarkerFaceColor',color(round(end/2),:))
xlabel('Raman Shift G band / cm^{-1}','FontSize',12);
ylabel('FWHM G^+ band / cm^{-1}','FontSize',12);
title('FWHM vs Shift G^+ band','FontSize',12);
set(gca,'Box','on');
legend show;

else
figure(1004), hold on;
plot(center_G(:),FWHM_G(:),'o','DisplayName',name{z},'color',color(round(end/2),:),'MarkerFaceColor',color(round(end/2),:))
xlabel('Raman Shift G band / cm^{-1}','FontSize',12);
ylabel('FWHM G band / cm^{-1}','FontSize',12);
title('FWHM vs Shift G band','FontSize',12);
set(gca,'Box','on');
legend show;

end
end
end

%% MAPS 

if map
raws_real=linspace(0,rawsmum,raws);
columns_real=linspace(0,colmum,col);

Imap=zeros(raws,col);
pos1map=zeros(raws,col);
pos2map=zeros(raws,col);
pos3map=zeros(raws,col);
I_2d_map=zeros(raws,col);

for i=1:raws
   for j=1:col
      Imap(i,j)=I21(col*(i-1)+j);%I21
      pos2map(i,j)=center_D(col*(i-1)+j);
      pos3map(i,j)=center_2D(col*(i-1)+j);
      if nt
      pos1map(i,j)=center_Gplus(col*(i-1)+j);  
      pos1minmap(i,j)=center_Gmin(col*(i-1)+j);   
      else
      pos1map(i,j)=center_G(col*(i-1)+j);
      end
   end
end

% Intensity ratio 1/2
figure(40), hold on
subplot(div(1),div(2),z); hold on
levels=10;
[M,c]=contourf(columns_real,raws_real,Imap,levels);
set(c,'LineColor','none')
xlabel('x / \mum','FontSize',12);
ylabel('y / \mum','FontSize',12);
%caxis([0,1]); % Set custom defined color range
h = colorbar;
title([name(z) 'I_d/I_g']);

% Position band 1
figure(41), hold on
subplot(div(1),div(2),z); hold on
levels=10;
[M1,c1]=contourf(columns_real,raws_real,pos1map,levels);
set(c1,'LineColor','none')
xlabel('x / \mum','FontSize',12);
ylabel('y / \mum','FontSize',12);
%caxis([1583,1602]); % Set custom defined color range
colorbar
if nt
title([name(z) 'G^+ band position']);
else
title([name(z) 'G band position']);
end

% Position band 2
figure(42), hold on
subplot(div(1),div(2),z); hold on
levels=10;
[M2,c2]=contourf(columns_real,raws_real,pos2map,levels);
set(c2,'LineColor','none')
xlabel('x / \mum','FontSize',12);
ylabel('y / \mum','FontSize',12);
%caxis([1339,1356]);
colorbar
title([name(z) 'D band position']);

% Postion band 3
figure(43), hold on
subplot(div(1),div(2),z); hold on
levels=10;
[M3,c3]=contourf(columns_real,raws_real,pos3map,levels);
set(c3,'LineColor','none')
xlabel('x / \mum','FontSize',12);
ylabel('y / \mum','FontSize',12);
%caxis([2670,2700]); % Set custom defined color range
colorbar
set(gca,'ColorScale','log')
title([name(z) '2D band position']);

if nt
% Postion band 1: G-
figure(44), hold on
subplot(div(1),div(2),z); hold on
levels=10;
[M4,c4]=contourf(columns_real,raws_real,pos1minmap,levels);
set(c4,'LineColor','none')
xlabel('x / \mum','FontSize',12);
ylabel('y / \mum','FontSize',12);
%caxis([2670,2700]); % Set custom defined color range
colorbar
set(gca,'ColorScale','log')
title([name(z) 'G^- band position']);
end

end

clearvars Shift Intensity IN

end