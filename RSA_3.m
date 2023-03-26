%% run RSA between 3 models and neural data (Dissimilarity matrix per ROI)
%% plot 3 figures(3 ROI): 3 mean of correlation in each figure (1 ROI)

clear all
addpath(genpath('/Users/muku/Documents/MATLAB folder/CoSMoMVPA-master/mvpa'));
addpath(genpath('/Users/muku/Documents/MATLAB folder/NifTI_20140122'));

data_filepath=fullfile('/Users/muku/Documents/MATLAB folder/Bracci_analysis/');

%% Load models
%mod = 1431(54*54 conditions) row and 3(number of models) columns
load([data_filepath, 'models/models.mat']);
mod=lower_models_vect;
%correlation between 3 models
r_mod=corr(lower_models_vect);

%% Load neural data
%The RDM file contains DSM of each subject per ROI
load ([data_filepath, 'results/SingleSubj_ROIs_RDM.mat']);
neuralData=RDM.data;

%create ROI variables to calculate the size of ROI
ROI = RDM.ROIs;
nROIs = size(ROI, 2);

%% RUN RSA WITH CORR AND PARTIALCORRI AND COMPARE THE RESULTS (test different correlation coefficients)
for r = 1:nROIs
 
    neural_vect=RDM.data{r};%get all the subjects DSM vector(11 in this case) of a specified ROI

    %correlation
    tempRSA_corr= corr(neural_vect,mod,'Type','Pearson');% correlate with 11 columns(DSM of all subjects and 3 column(3 models)  
    RSA_corr{r}=tempRSA_corr; % store corr results    
    tempRSA_parcorr=partialcorri(neural_vect,mod,'Type','Pearson')%partial correlation, which is more common
    RSA_parcorr{r}=tempRSA_parcorr;% store parcorr results
    
    %fisher transf the RSA results 
    fisher_tempRSA_corr=atanh(tempRSA_parcorr); %tempRSA_parcorr is 11 (subject) ×3(models) on a specified ROI    
    RSA_fisher_corr{r}=fisher_tempRSA_corr;%store results

    %ttest: two-tailed against zero correlation
    [h,p,ci,stats] = ttest(fisher_tempRSA_corr);
    stdevemean = std(fisher_tempRSA_corr)%standard deviation: how dispersed the data is in relation to the mean
    sem = std(fisher_tempRSA_corr)/sqrt(length(fisher_tempRSA_corr))%standard error of the mean: the dispersion of sample means around the population mean.
    stat.h=h;
    stat.p=p;
    stat.std=stdevemean;
    stat.sem=sem;
    stat.stats=stats;
    stat.rois= ROI{r};
         
    %compute noise ceiling
    for s = 1:size(neural_vect,2)
        singleSubj = neural_vect(:,s);
        maskMinus=neural_vect;maskMinus(:,s) = NaN;
        groupMinus=nanmean(maskMinus,2);
        acrossSubj_lower_bound(s,r)=corr(singleSubj,groupMinus);

    stat.noiseCeiling = mean(acrossSubj_lower_bound(:,r));    
    statistics{r}=stat;
    end
    
end

%% Save results
RSA.statistics=statistics;
RSA.results=RSA_fisher_corr;
RSA.results_mean = {mean(RSA_fisher_corr{1});mean(RSA_fisher_corr{2});mean(RSA_fisher_corr{3})}';%average correlation result among all the subjects

%output in the specific place withe the file name(RSA_results)
name_file=fullfile([[data_filepath,'results/'] 'RSA_results']);
save (name_file, 'RSA');

%% plot
%3 figures(3 ROI): each figure contains 3 bar(correlation with 3 models per ROI)
for i = 1:3
subplot(1,3,i)
b = bar(RSA.results_mean{i});
b.FaceColor = 'y'
title(ROI{i})
xlabel(ROI(i))
xticks([1 2 3])
xticklabels({'pixelwise','Shape','Category'})
xtickangle(45)
ylim([-0.02,0.45])
hold on
%plot standard error of the mean
errorbar(RSA.results_mean{i},statistics{1,i}.sem,'.','Color','r');

%plot p-value symbol
xt=get(gca,'XTick');
idx = xt(statistics{1,i}.h==1);
plot(idx,0.01,'*k','Color','r')
    if i == 1
        ylabel('z','FontSize',20)
    end

%plot noise celling
ceiling = statistics{1,i}.noiseCeiling
line([0 length(ROI)+1],[ceiling ceiling],'Color','r','LineWidth',5);
end
hold off


