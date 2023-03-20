clear all

addpath(genpath('/Users/muku/Documents/MATLAB folder/CoSMoMVPA-master/mvpa'));
addpath(genpath('/Users/muku/Documents/MATLAB folder/NifTI_20140122'));

% reset citation list
cosmo_check_external('-tic');

%% Load data
subjs =  {'SUB01','SUB02','SUB03','SUB04','SUB05','SUB09','SUB10','SUB11','SUB13','SUB14','SUB15'};
numSubjs = size (subjs, 2);

naROI =  {'BA17','aFG','pFG'};
nameROIs = size (naROI, 2);

ROI =  {'BA17.nii','aFG.nii','pFG.nii'};
numROIs = size (ROI, 2);

%set paths
study_path=fullfile('/Users/muku/Documents/MATLAB folder/Bracci_analysis/data_1_advancedfmri/'); %fmri data
mask_path=fullfile('/Users/muku/Documents/MATLAB folder/Bracci_analysis/ROI/'); % ROI
results_path=fullfile('/Users/muku/Documents/MATLAB folder/Bracci_analysis/results/');

%beta_index
% a = [];
% n = 1 ;
% while n+53 < 977
%     a = [a,n:n+53] %append
%     n = n + 53 + 7; %add 53 and skip 7
% end
% beta_index = a'

for s = 1:numSubjs
    for r = 1:numROIs
        
        data_path=fullfile([study_path, subjs{s} ]); %sub 
        mask_fn=fullfile([mask_path, ROI{r}]); % ROI           
        
        data_fn=fullfile(data_path,'SPM.mat');
        
        %% TASK: generate COSMO DS structure %%%        
        ds = cosmo_fmri_dataset(data_fn,'mask',mask_fn);

        %% TASK: Create two additional structures for the dimensions of (1) categories and (2) shapes 
        % TASK: define targets and chunks %%%
        ds.sa.targets = repmat((1:54)',16,1);                
        ds.sa.chunks = reshape(repmat((1:16),54,1),[],1);  
        % TASK: define the category dimension %%%
        ds.sa.cat_dim = repmat((1:6)',9*16,1);
        % TASK: define the shape dimension %%%
        ds.sa.shape_dim = reshape(repmat((1:9),6,16),[],1);

        % simple sanity check to ensure all attributes are set properly
        cosmo_check_dataset(ds);
        
        % print dataset
        fprintf('Dataset input:\n');
        cosmo_disp(ds);
        
        % remove constant features
        ds=cosmo_remove_useless_data(ds);
        
        % print dataset
        fprintf('Dataset input:\n');
        cosmo_disp(ds);
        
        name_file=fullfile([results_path, subjs{s}, '_' naROI{r} '_ds']);
        save (name_file, 'ds');
    end
end