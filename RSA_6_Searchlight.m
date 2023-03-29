%% Set data paths
clear all

addpath(genpath('/Users/muku/Documents/MATLAB folder/CoSMoMVPA-master/mvpa'));
addpath(genpath('/Users/muku/Documents/MATLAB folder/NifTI_20140122'));
addpath(genpath('/Users/muku/Documents/MATLAB folder/BrainNetViewer'));

data_filepath=fullfile('/Users/muku/Documents/MATLAB folder/Bracci_analysis/');

study_path='/Users/muku/Documents/MATLAB folder/Bracci_analysis/data_1_advancedfmri/';
output_path='/Users/muku/Documents/MATLAB folder/Bracci_analysis/results/';
%% Load data

subjs =  {'SUB01','SUB02','SUB03','SUB04','SUB05','SUB09','SUB10','SUB11','SUB13','SUB14','SUB15'};
numSubjs = size (subjs, 2);

for s=1:numSubjs
    
    %　Get the data name of each subject = /data_1_advancedfmri/SUB01
    data_path=fullfile([study_path subjs{s}]);
    
    %　set whole-brain mask
    mask_fn=fullfile(data_path, 'mask.img'); % whole brain mask
    
    %　set targets and chunks
    targets=repmat(1:54,1,16);    %targets = number of conditions (1:54(conditions),1(column),16(runs))
    chunks=floor(((1:864)-1)/54)+1; %floor(1.9) = 1, round toward negative integer, so create 54 times of 1
    
    % generate DS / each subject has SPM.mat
    data_fn=fullfile(data_path,'SPM.mat');
    ds = cosmo_fmri_dataset(data_fn,'mask',mask_fn, 'targets',targets,'chunks',chunks); % define targets and chunks
        
    % remove 'useless' (constant and/or non-finite) samples or features
    ds=cosmo_remove_useless_data(ds);
    
    % average across runs so that we have one value for each condition 864
    % to 54 
    ds=cosmo_fx(ds, @(x)mean(x,1), 'targets', 1);
    
    % simple sanity check to ensure all attributes are properly set
    cosmo_check_dataset(ds);
    
    % print dataset
    fprintf('Dataset input:\n');
    cosmo_disp(ds);
    
    % Define feature neighorhoods for each feature (voxel)
    nvoxels_per_searchlight=100;
    
    fprintf('Defining neighborhood for each feature\n');
    
    % set feature neighorhoods and searchilight size
    nbrhood=cosmo_spherical_neighborhood(ds,'count',nvoxels_per_searchlight);
    
    % print neighborhood
    fprintf('Searchlight neighborhood definition:\n');
    cosmo_disp(nbrhood);
    
    % load models
    models =  {'shape','semantic'};
    numMod = size (models, 2);
    
    for m = 1:numMod
        model_fn=fullfile('/Users/muku/Documents/MATLAB folder/Bracci_analysis/models/'); %model data path -- set it
        filename=fullfile([model_fn, models{m}]);
        load (filename, 'dissimilarity');
        target_dsm=dissimilarity;
        
        % run searchlight RSA
        fprintf('Using the following target dsm\n');
        
        t_dsm=target_dsm+target_dsm';% create symetrical matrix to match with the ds
        disp(t_dsm);
        imagesc(t_dsm);
        nsamples=size(ds.samples,1);
        set(gca,'XTick',1:nsamples,'XTickLabel',ds.sa.labels,...
            'YTick',1:nsamples,'YTickLabel',ds.sa.labels)
        
        % set measure 
        measure=@cosmo_target_dsm_corr_measure;
                
        % define your additional args (1)center the data; (2)define your model (to be done)
        measure_args=struct('center_data',true);
        measure_args.target_dsm= t_dsm;
        
        % run searchlight
        % ds = brain_dataset, nbrhood = searchlight_size, measure =
        % correlation, measure_args = define model and center data.
        ds_rsm=cosmo_searchlight(ds,nbrhood,measure,measure_args);
                
        disp('Done Searchlight for model');
        % show results
        cosmo_plot_slices(ds_rsm);
        
        % store results as nii.file
        output_fn=fullfile([output_path, subjs{s},'_' models{m},'_rsm.nii']);
        cosmo_map2fmri(ds_rsm,output_fn);
        
        % save results
        cd(output_path);
        
        % store result as mat.file
        filename_ds_rsm=sprintf([subjs{s},'_' models{m},'_ds_rsm']);
        save(filename_ds_rsm,'ds_rsm');
        
    end
end



%% compute maps for the individual models

study_path2='/Users/muku/Documents/MATLAB folder/Bracci_analysis/results/';

for m=1:numMod
    for s=1:numSubjs      
        data_searchlight=fullfile([study_path2, subjs{s} '_' models{m} '_ds_rsm.mat']); % download mat.file of each subject's each ROI
        ds_cell{s}=cosmo_fmri_dataset(data_searchlight); %ds_cell contains 11 subjects, 1 column for 1 subject       
    end
    
   % compute voxels intersect with cosmo_mask_dim_intersect / only
   % overlap voxel is extracted
   [idx_cell,ds_intersect_cell]=cosmo_mask_dim_intersect(ds_cell,2); %2 indicates each columns fusion
    
    
    %%COMBINE DATASET %%
    for subject_i=1:numSubjs

        % get dataset
        ds=ds_intersect_cell{subject_i}; % ds_intersect_cell = 11×1
        
        % assign chunks
        n_samples=size(ds.samples,1);
        ds.sa.chunks=ones(n_samples,1)*subject_i;
        
        % assign targets: see cosmo_stat on how to assign targets depending on subsequent analysis
        % (one-sample or two-sample t-test, or one-way or repeated-measures ANOVA).
        
        ds.sa.targets=repmat(ones,1); % one sample t-test
              
        % store results
        ds_intersect_cell{subject_i}=ds;
    end
    
    % stack subjects data with cosmo_stack
    ds_all=cosmo_stack(ds_intersect_cell,1,'drop_nonunique'); %1 indicates dim, stack .samples and .sa
    
    % SAVE DATA MAPS FOR DIRECT CONTRAST
    cd(study_path);
    filename=sprintf([models{m},'_groupMaps.mat']);
    save(filename,'ds_all');
     
    %%RUN STATS %%
    % one-sample t-test against 0 (for representational similarity analysis)
    h0_mean=0;
    
    % number of null iterations.
    % values of at least 10,000 are recommended for publication-quality
    niter=10000;
    
    % Set neighborhood for clustering / just one voxel significant is
    % considered as noise, 3 indicates two more than two voxels cluster
    cluster_nbrhood=cosmo_cluster_neighborhood(ds_all,'fmri',3);
    
    % Show a plot with the sorted number of neighbors for each voxel
    n_neighbors_per_feature=cellfun(@numel,cluster_nbrhood.neighbors);
    plot(sort(n_neighbors_per_feature))
    
    % correction for multiple comparison with tfce >> threshold-free cluster enhancement ref: Stephen M. Smith, Thomas E. Nichols (2009)   
    stat_map_z=cosmo_montecarlo_cluster_stat(ds_all, cluster_nbrhood,...
        'niter', niter, 'h0_mean', h0_mean);
    cosmo_plot_slices(stat_map_z)
      
    x = stat_map_z.samples;
    indices = find(abs(x) <1.6449); %threshold z > 1.9600 For a two-tailed test
    x(indices) = 0; % lower than z = 1.6449 is 0, not siginificant
    
    stat_map_z.samples=x;
    
    % Cluster correction - store results
    filename=sprintf([models{m},'_groupRSA_MontecarloCluster.nii']);
    output_fn=fullfile(study_path2,filename);
    % save SPM
    cosmo_map2fmri(stat_map_z,output_fn);
    
    % T-stats uncorrected results
    t_map=cosmo_stat(ds_all,'t');
    cosmo_disp(t_map);
    
    filename=sprintf([models{m},'_groupRSA_uncorrected.nii']);
    output_fn=fullfile(study_path2,filename);
    % save SPM
    cosmo_map2fmri(t_map,output_fn);
    
end


%%Contrast maps %%
% 11 (number of subject) × 48805 (number of voxels) in each model (shape
% and semantic)
groupdata_shape=fullfile([study_path2, 'shape_groupMaps.mat']);
groupdata_semantic=fullfile([study_path2, 'semantic_groupMaps.mat']);

ds_shape=cosmo_fmri_dataset(groupdata_shape);
ds_semantic=cosmo_fmri_dataset(groupdata_semantic);

%not interested in negative correlation, transformed negative value into 0
x = ds_shape.samples;
xindices = find(x < -0.00);
x(xindices) = 0;
ds_shape.samples=x;

q = ds_semantic.samples;
qindices = find(q < -0.00);
q(qindices) = 0;
ds_semantic.samples=q;

% shape vs category
data_test=ds_semantic;
data_test.samples=[];
data_test.samples = ds_semantic.samples - ds_shape.samples;

     
%%RUN STATS %%
% one-sample t-test against 0 (for representational similarity analysis)
h0_mean=0;

% number of null iterations.
% values of at least 10,000 are recommended for publication-quality
niter=1000;

% Set neighborhood for clustering
cluster_nbrhood=cosmo_cluster_neighborhood(data_test,'fmri',3);

% Show a plot with the sorted number of neighbors for each voxel
n_neighbors_per_feature=cellfun(@numel,cluster_nbrhood.neighbors);
plot(sort(n_neighbors_per_feature))

%correction for multiple comparison with tfce >> threshold-free cluster enhancement ref: Stephen M. Smith, Thomas E. Nichols (2009)
stat_map_z=cosmo_montecarlo_cluster_stat(data_test, cluster_nbrhood,...
'niter', niter, 'h0_mean', h0_mean);
cosmo_plot_slices(stat_map_z)

x = stat_map_z.samples;
indices = find(abs(x) <1.6449); %threshold z > 1.9600 For a two-tailed test
x(indices) = 0;
stat_map_z.samples=x;

%Cluster correction - store results
filename='category_vs_shape_RSA_MontecarloCluster.nii';
output_fn=fullfile(study_path2,filename);
% save SPM
cosmo_map2fmri(stat_map_z,output_fn);
