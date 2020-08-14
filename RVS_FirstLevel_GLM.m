close all; clear all; clc

%% Preparation
% type in the project direction
project_dir = 'RVS_EEG';

subj_ID     = {'01', '02', '03', '04', '09', '10', '11', '12', '13', '14', ...
               '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', ...
               '25', '26', '27', '28'};
           
analysis_dir = fullfile(project_dir, 'processed_data_1');
cd(analysis_dir);          
cond_file = {'rehearsal'};      
filter = {'^scondition_..\.nii$', '^scondition_....\.nii$'};
spm('defaults', 'EEG');

%%
for sj = 1:length(subj_ID)    
    images = {};
    for j = 1:2
        subj_dir = fullfile(analysis_dir, ['sj_' subj_ID{sj}], cond_file{1});        
        % find all images
        images = [images; cellstr(spm_select('ExtFPList', [subj_dir '/images'], filter{j}))];
    end           
        tmp_corr = regexp(images, ['(?<!\d)(\d{2})(?=.nii)'], 'match');
        conds_corr= str2double([tmp_corr{:}])';
        conds_corr(:,2) = conds_corr/10;
        conds_corr(:,1) = 1;
        cons_corr =num2str(conds_corr);
        tmp_error  = regexp(images, ['(?<!\d)(\d{4})(?=.nii)'], 'match');
        conds_error= str2double([tmp_error{:}])';
        cons_error =num2str(conds_error);
        conds_error(:,2) = (conds_error-1000)/10;
        conds_error(:,1) = 2;
        conditions = [conds_corr' conds_error'];
        
        stats_dir = fullfile(subj_dir, 'stats');
        if exist(stats_dir, 'dir') ~= 7
            mkdir(stats_dir)
        end
        %%
        %************************************************************************
        matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(stats_dir);
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject.scans = images;%cellstr(images);
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject.conds = conditions';
        
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name = 'frequency';
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept = 1;               % 0 = iid
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 0;           % 1 = unequal
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;              % Grand mean scaling
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;
        
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name = 'performance';
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept = 1;               % 0 = iid
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 0;           % 1 = unequal
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;              % Grand mean scaling
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;
        
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.inter.fnums = [1 2];
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;        
        spm_jobman('run', matlabbatch)
        clear matlabbatch
        
        %% ========================================================================
        %               MODEL ESTIMATION
        % =========================================================================
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(stats_dir, 'SPM.mat')};
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;        
        display('MODEL ESTIMATION');
        spm_jobman('run', matlabbatch)
        clear matlabbatch
    
        %% ========================================================================
        %                        COMPUTE CONTRASTS
        %==========================================================================
        matlabbatch{1}.spm.stats.con.spmmat = {fullfile(stats_dir, 'SPM.mat')};
        matlabbatch{1}.spm.stats.con.delete = 1;        % delete existing contrasts
        
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'paramatric_pos_correct';
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights =  [-0.75 -0.25 0.25 0.75];
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'paramatric_pos_incorrect';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights =  [0 0 0 0 -0.75 -0.25 0.25 0.75 ];
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'paramatric_pos';
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights =  [-0.75 -0.25 0.25 0.75 -0.75 -0.25 0.25 0.75 ];
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
        
        matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'correct > incorrect';
        matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights =  [1 1 1 1 -1 -1 -1 -1];
        matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
        
        matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'incorrect > correct';
        matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights =  [-1 -1 -1 -1 1 1 1 1];
        matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
        
        matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'pos_correct';
        matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights =  [1 1 1 1 0 0 0 0];
        matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
        
        matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = 'pos_incorrect';
        matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights =  [0 0 0 0 1 1 1 1];
        matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
        
        % F contrast
        matlabbatch{1}.spm.stats.con.consess{8}.fcon.name = 'omnibus';
        matlabbatch{1}.spm.stats.con.consess{8}.fcon.weights =  eye(8);
        matlabbatch{1}.spm.stats.con.consess{8}.fcon.sessrep = 'none';
        
        display('Compute contrasts');
        spm_jobman('run', matlabbatch)
        clear matlabbatch D conditions images subj_dir stats_dir   
end

%% Plotting
hdr = spm_vol('spmT_0001.nii');
[Y, XYZ] = spm_read_vols(hdr);
load('chan.mat');

% plot the t-contrast in time frequency domain
h = figure(2);
h = imagesc(Y(:,:,47));
axis xy

set(gca,'XTick',[1 6 11 16 21 26 31])
set(gca,'XTicklabel',{'-0.25', '0', '0.25', '', '', '1.00', '1.25'})
set(gca,'YTick',[1.5 4 6.5 9 11.5 14 16.5 19 21.5])
set(gca,'YTicklabel',{'5', '10', '15', '20', '25', '30', '35', '40', '45'})

% plot for the all stimlus
set(gca,'XTick',[0 1.25 2.5 3.75 5 6.25 7.5 8.75 10 11.25])
set(gca,'XTicklabel',{'0', '1.25', '2.5', '3.75', '5', '6.25', '7.5', '8.75', '10', '11.25'})
line([0,0], [0, 48], 'Color', 'black', 'LineStyle', '--', 'Linewidth', 2)
line([1.25,1.25], [0, 48], 'Color', 'black', 'LineStyle', '--', 'Linewidth', 2)
line([2.5,2.5], [0, 48], 'Color', 'black', 'LineStyle', '--', 'Linewidth', 2)
line([3.75,3.75], [0, 48], 'Color', 'black', 'LineStyle', '--', 'Linewidth', 2)
line([5,5], [0, 48], 'Color', 'black', 'LineStyle', '--', 'Linewidth', 2)
line([6.25,6.25], [0, 48], 'Color', 'black', 'LineStyle', '--', 'Linewidth', 2)
line([7.5,7.5], [0, 48], 'Color', 'black', 'LineStyle', '--', 'Linewidth', 2)
line([8.75,8.75], [0, 48], 'Color', 'black', 'LineStyle', '--', 'Linewidth', 2)
line([10,10], [0, 48], 'Color', 'black', 'LineStyle', '--', 'Linewidth', 2)
line([11.25,11.25], [0, 48], 'Color', 'black', 'LineStyle', '--', 'Linewidth', 2)

box off
line([0,0], [0, 48], 'Color', 'black', 'LineStyle', '--', 'Linewidth', 2)
line([11,11], [0, 48], 'Color', 'black', 'LineStyle', '--', 'Linewidth', 2)
line([31,31], [0, 48], 'Color', 'black', 'LineStyle', '--', 'Linewidth', 2)
colorbar

%% Check the averaged electrodes 
imagesc(D.time, D.frequencies, squeeze(mean(D(D.indchannel({'C2', 'C4', 'C6', 'CP2' 'CP6', 'CP4'}),:,:),1)));
axis xy
colorbar

%% Plot all the channels 
mask = '^mrtf.*\.mat';
D = spm_eeg_load('plotting');

hdr = spm_vol('spmT_0001.nii');
[Y, XYZ] = spm_read_vols(hdr);
Y_tmp = permute(Y,[3 1 2]);
for i = 1:8
    Y_vision(:,:,:,i) = Y_tmp;
end
D(:,:,:,:) = Y_vision;
clear Y_tmp Y_vision Y

%% 
data = fttimelock(D);
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.colorbar = 'yes';
cfg.baseline = 'yes';
cfg.colormap = jet;
cfg.baselinetype = 'relchange';
ft_multiplotTFR(cfg, data);
