%% -- compute fMRI RDM
%% one dissimilarity matrix for each subject per ROI
%% plot the average dissimilarity matrix per ROI: In this case, 3 ROI = 3 DMS matrix

clear all

addpath(genpath('/Users/muku/Documents/MATLAB folder/CoSMoMVPA-master/mvpa'));
addpath(genpath('/Users/muku/Documents/MATLAB folder/NifTI_20140122'));

% define subjs and ROIs 
subjs =  {'SUB01','SUB02','SUB03','SUB04','SUB05','SUB09','SUB10','SUB11','SUB13','SUB14','SUB15'};
numSubjs = size (subjs, 2);

ROI =  {'BA17','aFG','pFG'};
resROIs = size (ROI, 2);

%define path to DS folder
ds_filepath=fullfile('/Users/muku/Documents/MATLAB folder/Bracci_analysis/results/'); %DS data

for r= 1:resROIs
    for s = 1:numSubjs
        
        %load DS
        filename=fullfile([ds_filepath, subjs{s}, '_' ROI{r} '_ds']);
        load (filename, 'ds');
        
        % print dataset
        fprintf('Dataset input:\n');
        cosmo_disp(ds);       
       
        %compute avg across runs with cosmo cosmo_fx
        f_ds=cosmo_fx(ds,@(x)mean(x,1),'targets'); %average within the same targets
        
        %compute RDM with cosmo_dissimilarity_matrix_measure (use correlation measure)
        ds_dsm = cosmo_dissimilarity_matrix_measure(f_ds,'metric','correlation','center_data',true); % center data = normarize %pearson correlation
        [unfl,labels,values]=cosmo_unflatten(ds_dsm,1,'set_missing_to',NaN); %remove diagonal part

        %store results dsm (vector);
        temp_dsm=ds_dsm.samples;
        dsm_all(:,s)=temp_dsm;
        
        %store results matric;
        temp_dsm_mat=unfl;
        dsm_mat_all(:,:,s)=temp_dsm_mat;
        
    end    
    dsm_vect{r}=dsm_all;%dsm_vect = row: 1431(54*54), column: 11(number of subjects)
    dsm_unflatten{r}=dsm_mat_all;%dsm_unflatten = DSM matrix per subejcts each ROI 
end

%store results
RDM.data=dsm_vect;
RDM.data_unflatten=dsm_unflatten;
RDM.ROIs=ROI;
RDM.SUB=subjs;

%output RDM in the result folder
name_file=fullfile([ds_filepath, 'SingleSubj_ROIs_RDM']);
save (name_file, 'RDM');

%% plot
for i=1:size(RDM.data_unflatten,2)
    temp_mat=RDM.data_unflatten{i}
    mean_rdm=mean(temp_mat,3); %3 indicates (:,:,here),resulting in average across subjects
    subplot (1,3,i)
    imagesc(mean_rdm)
    axis equal tight
    title(ROI(i),'averaged across subjects')
    xlabel('conditions')
    ylabel('conditions')
    hold on
    set(gcf,'color','w');
    colorbar
end