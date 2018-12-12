format compact;
close all;
fclose all;
format short;

%% Extended to Sample 30
%% Composed in Level2_RO2REML30_21Nov2018.m
%data = importdata('Level2_RO2REML30Block_9Nov2018.csv',',');

%% Extended to Sample 50
%% Composed in Level2_RO2REML30_21Nov2018.m
data = importdata('Level2_RO2REML50Block_9Nov2018.csv',',');

%% Extended to Sample 100
%% Composed in Level2_RO2REML100_21Nov2018.m
%data = importdata('Level2_RO2REML100Block_9Nov2018.csv',',');

%% Extended to Sample 200
%% Composed in Level2_RO2REML200_21Nov2018.m
%data = importdata('Level2_RO2REML200Block_9Nov2018.csv',',');

%% Extended to Sample 500
%% Composed in Level2_RO2REML500Part1_21Nov2018.csv
%% Composed in Level2_RO2REML500Part2_21Nov2018.csv
%data=[];
%data = [data; importdata('Level2_RO2REML500Part1Block_9Nov2018.csv',',')];
%data = [data; importdata('Level2_RO2REML500Part2Block_9Nov2018.csv',',')];


% D(:,1) = n
% D(:,2) = g
% D(:,3) = ng_reml (extended sample size)
% D(:,4) = Case (factor variability)
% D(:,5) = scen (error variability)
% D(:,6) = ES
% D(:,7) = e_var
% D(:,8) = b00 estimate
% D(:,9) = b01 estimate 
% D(:,10) = b10 estimate 
% D(:,11) = b11 estimate 
% D(:,12) = b00 SE 
% D(:,13) = b01 SE
% D(:,14) = b10 SE
% D(:,15) = b11 SE
% D(:,16) = b00 tstat
% D(:,17) = b00 tstat
% D(:,18) = b01 tstat
% D(:,19) = b10 tstat
% D(:,20) = b11 pvalue
% D(:,21) = b01 pvalue
% D(:,22) = b10 pvalue
% D(:,23) = b11 pvalue
% D(:,24) = r0j var
% D(:,25) = r1j var
% D(:,26) = MSE 

D=data;

stride=1000;
fprintf(1,'i,n,g,reml_ng,factor var,err var,ES,MSE,b00 mean,b01 mean,b10 mean,b11 mean,b00 median,b01 median,b10 median,b11 median,b00 SE,b01 SE,b10 SE,b11 SE,b00 tstat,b01 tstat,b10 tstat,b11 tstat,b00 pvalue,b01 pvalue,b10 pvalue,b11 pvalue,rel b00 bias,rel b01 bias,rel b10 bias,rel b11 bias,b00 power,b01 power,b10 power,b11 power \n')

for i = 0:stride:(size(D,1)-stride)

ES = mean(D(i+(1:stride),6));

%b00
b00_mean=mean(D(i+(1:stride),8));
b00_median=median(D(i+(1:stride),8));
b00_SE=mean(D(i+(1:stride),12));
b00_tstatAVG=b00_mean/b00_SE;
%b00_tstat=mean(D(i:i+stride,16));
b00_pvalue=mean(D(i+(1:stride),20));
b00_bias = (abs(b00_mean - ES))/ES;
b00_power = sum(D(i+(1:stride),20)<.05)/1000;

% b01
b01_mean=mean(D(i+(1:stride),9));
b01_median=median(D(i+(1:stride),9));
b01_SE=mean(D(i+(1:stride),13));
b01_tstatAVG=b01_mean/b01_SE;
b01_tstat=mean(D(i+(1:stride),17));
b01_pvalue=mean(D(i+(1:stride),21));
b01_bias = (abs(b01_mean - ES))/ES;
b01_power = sum(D(i+(1:stride),21)<.05)/1000;

% b10
b10_mean=mean(D(i+(1:stride),10));
b10_median=median(D(i+(1:stride),10));
b10_SE=mean(D(i+(1:stride),14));
b10_tstatAVG=b10_mean/b10_SE;
b10_tstat=mean(D(i+(1:stride),18));
b10_pvalue=mean(D(i+(1:stride),22));
b10_bias = (abs(b10_mean - ES))/ES;
b10_power = sum(D(i+(1:stride),22)<.05)/1000;

% b11
b11_mean=mean(D(i+(1:stride),11));
b11_median=median(D(i+(1:stride),11));
b11_SE=mean(D(i+(1:stride),15));
b11_tstatAVG=b11_mean/b11_SE;
b11_tstat=mean(D(i+(1:stride),19));
b11_pvalue=mean(D(i+(1:stride),23));
b11_bias = (abs(b11_mean - ES))/ES;
b11_power = sum(D(i+(1:stride),23)<.05)/1000;

% MSE
sqrt_mse = mean(D(i+(1:stride),26));

fprintf(1,'%d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f \n',...
    ceil((i+1)/stride),D(i+1,1),D(i+1,2),D(i+1,3),D(i+1,4),D(i+1,5),D(i+1,6),sqrt_mse,...
    b00_mean, b01_mean, b10_mean, b11_mean,...
    b00_median, b01_median, b10_median, b11_median,... 
    b00_SE, b01_SE, b10_SE, b11_SE,...
    b00_tstatAVG, b01_tstatAVG, b10_tstatAVG, b11_tstatAVG,...
    b00_pvalue, b01_pvalue, b10_pvalue, b11_pvalue,... 
    b00_bias, b01_bias, b10_bias, b11_bias,...
    b00_power, b01_power, b10_power, b11_power)
end

