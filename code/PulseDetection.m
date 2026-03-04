% Pulse detection Script from BOLD DATA 
% TUE008 NIMG Study - Corinna Schulz - 2025 

% Newbold paper cast spike detection
clc 
clear all 

% Set filepath for preprocessed DATA from conn 
path_prepro = '/mnt/ghrelin/TUE008_NIMG/Study/derivatives/conn/finalsampleCONN/data/BIDS/dataset/'; 

% Set output path 
Pulse_output = '/mnt/TUE_data/TUE008/Data/TUE008_NIMG/Study/derivatives/conn/finalsampleCONN/Pulse_Detection/'; 

% Load atlas file for ROIs 
atlas_TUE008 = niftiread('/mnt/ghrelin/TUE008_NIMG/Atlas/rrHOex_atlas_TUE008_NIMG.nii'); 
% Get number of regions 
regions = (unique(atlas_TUE008(~isnan(atlas_TUE008(:)))));
% Set SEED Region Hypothalamus, set atlas Index (look up label file) 
region_hypoL = 148; 
region_hypoR = 147; 
region_NaCL = 105; 
region_NaCR = 104;
region_VTA_L = 143; 
region_VTA_R  = 144; 
region_Cau_L = 94; 
region_Cau_R = 95; 
region_Put_L = 96; 
region_Put_R  = 97; 
region_Ins_L = 3; 
region_Ins_R  = 4; 

regions_to_do = [region_hypoL,region_hypoR,region_NaCL, region_NaCR, ...
                region_VTA_L, region_VTA_R, ...
                region_Cau_L, region_Cau_R, region_Put_L, region_Put_R, ...
                region_Ins_L, region_Ins_R];

% Set number of phases
% BASELINE, BOLUS, INFUSION 1, INFUSION 2, IMT
phases = 10; % attention, conn labelled all phases (from 2 sessions) from 1-10 with IMT as run-01 and run-06 
subjects = 26; 


% Get condition list 
matching = readtable('/mnt/TUE_general/Projects/TUE008_DFG_Ghrelin/10_analysis/CONN_STUDY_ID_MATCHING.xlsx');
cond_ghrelin = readtable('/mnt/TUE_general/Projects/TUE008_DFG_Ghrelin/03_Administration/0.2 NIMG/4.randomization/TUE008_NIMG_conditions_Final.xlsx');


% Init FC matrix to save longformat  
Pulse_all = [];

%% (1) For each subject get Hypothalamus ROI extracted
for sub = 1:(subjects)
    fprintf("Subject %d\n", sub);

    % For each subject calculate ROIs for all 5 Phases (and two sessions)
    % BASELINE, BOLUS, INFUSION 1, INFUSION 2, IMT

    % Initialise condition counter
    cond = 0;

    for pha = 1:phases % Conn labelled Phases for Sessions conseq. 1-10 
        fprintf("Phase %d\n", pha);

        cond = cond + 1;
        % Get filename for Voxel Matrix
        if sub < 10 && pha < 10
            mat_file_name = strcat('sub-000', char(num2str(sub)),'/func/dssub-000' ,char(num2str(sub)) , '_run-0' ,char(num2str(pha)), '_bold.nii');
        elseif sub >= 10 && pha < 10
            mat_file_name = strcat('sub-00', char(num2str(sub)),'/func/dssub-00' ,char(num2str(sub)) , '_run-0' ,char(num2str(pha)), '_bold.nii');
        elseif sub < 10 && pha == 10
            mat_file_name = strcat('sub-000', char(num2str(sub)),'/func/dssub-000' ,char(num2str(sub)) , '_run-' ,char(num2str(pha)), '_bold.nii');
        elseif sub >= 10 && pha == 10
            mat_file_name = strcat('sub-00', char(num2str(sub)),'/func/dssub-00' ,char(num2str(sub)) , '_run-' ,char(num2str(pha)), '_bold.nii');

        end


        % Load the denoised and smoothed nii
        nii_denoised = niftiread(fullfile(path_prepro, mat_file_name ));
        % size(nii_denoised) Dimensions: 91, 109, 91, 426

        nii_denoised(nii_denoised == 0) = NaN;


        % Get ROI Value for Hypothalamus
        % Get ROI Value for ROI

        % regions = regions(2:end); %discard unique 0 value, not a ROI

        % Get number of volumes/ timepoints
        volumes = size(nii_denoised,4);

      
        % Loop through every Atlas region and get timeseries for that region
        % by looping thorugh the Volumes and averaging those voxels activity
        for t = 1:volumes
            run = nii_denoised(:,:,:,t);
            ROI_timeseries_HR(t) = nanmean(run(atlas_TUE008 == 147));
            ROI_timeseries_HL(t) = nanmean(run(atlas_TUE008 == 148));
            ROI_timeseries_hypo(t) = (ROI_timeseries_HR(t)+ ROI_timeseries_HL(t))/2; 
            
            ROI_timeseries_Ins_R(t) = nanmean(run(atlas_TUE008 == 3));
            ROI_timeseries_INS_L(t) = nanmean(run(atlas_TUE008 == 4));
            ROI_timeseries_INS(t) = (ROI_timeseries_Ins_R(t)+ ROI_timeseries_INS_L(t))/2; 
            
            ROI_timeseries_CA_R(t) = nanmean(run(atlas_TUE008 == 94));
            ROI_timeseries_CA_L(t) = nanmean(run(atlas_TUE008 == 95));
            ROI_timeseries_Caudate(t) = (ROI_timeseries_CA_R(t)+ ROI_timeseries_CA_L(t))/2; 
            
                       ROI_timeseries_put_R(t) = nanmean(run(atlas_TUE008 == 96));
            ROI_timeseries_put_L(t) = nanmean(run(atlas_TUE008 == 97));
            ROI_timeseries_Putamen(t) = (ROI_timeseries_put_R(t)+ ROI_timeseries_put_L(t))/2; 
            
            
            ROI_timeseries_NaCR(t) = nanmean(run(atlas_TUE008 == 104));
            ROI_timeseries_NaCL(t) = nanmean(run(atlas_TUE008 == 105));
            ROI_timeseries_NAcc(t) = (ROI_timeseries_NaCR(t)+ ROI_timeseries_NaCL(t))/2; 
            
             ROI_timeseries_VTA_L(t) = nanmean(run(atlas_TUE008 == 143));
            ROI_timeseries_VTA_R(t) = nanmean(run(atlas_TUE008 == 144));
            ROI_timeseries_VTA(t) = (ROI_timeseries_VTA_L(t)+ ROI_timeseries_VTA_R(t))/2; 
                 
                         ROI_timeseries_NTS(t) = nanmean(run(atlas_TUE008 == 153));

            
        end
        

        % Save all ROI Timeseries per Phase, Roi, Timepoint
        All_Phases_Hypo_timeseries{pha} = ROI_timeseries_hypo(:);
        clear ROI_timeseries_hypo
        All_Phases_NAcc_timeseries{pha} = ROI_timeseries_NAcc(:);
        clear ROI_timeseries_NAcc
        All_Phases_VTA_timeseries{pha} = ROI_timeseries_VTA(:);
        clear ROI_timeseries_VTA
        All_Phases_Insula_timeseries{pha} = ROI_timeseries_INS(:);
        clear ROI_timeseries_INS
        All_Phases_Caudate_timeseries{pha} = ROI_timeseries_Caudate(:);
        clear ROI_timeseries_Caudate
        All_Phases_Putamen_timeseries{pha} = ROI_timeseries_Putamen(:);
        clear ROI_timeseries_Putamen
        All_Phases_NTS_timeseries{pha} = ROI_timeseries_NTS(:); 

    end

    All_Subs_Phases_HYPO_timeseries{sub} = All_Phases_Hypo_timeseries;
    All_Subs_Phases_NAcc_timeseries{sub} = All_Phases_NAcc_timeseries;
    All_Subs_Phases_VTA_timeseries{sub} = All_Phases_VTA_timeseries;

    All_Subs_Phases_INS_timeseries{sub} = All_Phases_Insula_timeseries;
    All_Subs_Phases_CAUD_timeseries{sub} = All_Phases_Caudate_timeseries;
    All_Subs_Phases_PUT_timeseries{sub} = All_Phases_Putamen_timeseries;

        All_Subs_Phases_NTS_timeseries{sub} = All_Phases_NTS_timeseries;

end


save("All_Subs_Phases_timeseries.mat","All_Subs_Phases_NTS_timeseries", "All_Subs_Phases_VTA_timeseries", "All_Subs_Phases_HYPO_timeseries", "All_Subs_Phases_NAcc_timeseries","All_Subs_Phases_PUT_timeseries","All_Subs_Phases_CAUD_timeseries","All_Subs_Phases_INS_timeseries");


% Now plot example 

%% --- plot phases 1:5 vs 6:10 for one example subject ---
sub_idx = 16;   

Matching_ID = matching.STUDY_ID(matching.CONN_ID(:) == sub_idx)
Ghrelin = cond_ghrelin.Ghrelin(cond_ghrelin.ID(:) == Matching_ID)
gh_flag = Ghrelin(1); 

phase_ts = All_Subs_Phases_VTA_timeseries{sub_idx};  
%phase_ts = All_Phases_Hypo_timeseries; 
phase_ts = All_Subs_Phases_INS_timeseries{sub_idx};  


% Colors
cG = [0.11 0.62 0.12];   % Ghrelin
cP = [0.20 0.45 0.74];   % Placebo

% Helper: mean-normalize to 0
normalize_to_mean = @(x) (x - mean(x,'omitnan'));

% Peak marking toggle
mark_peaks = true;

figure('Color','w','Name',sprintf('Subject %02d — Ghrelin vs Placebo (mean-normalized)', sub_idx));
tl = tiledlayout(5,1,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('Subject %02d — Ghrelin (green) vs Placebo (blue) — mean-normalized, overlaid per phase pair', sub_idx));

for k = 1:5
    a = phase_ts{k};       % phase 1..5 (session A)
    b = phase_ts{k+5};     % phase 6..10 (session B)

    nexttile; hold on;

    if isempty(a) || isempty(b)
        axis off; title(sprintf('Phase %d vs %d (no data)', k, k+5));
        continue;
    end

    % Mean-normalize
    aN = normalize_to_mean(a(:));
    bN = normalize_to_mean(b(:));

    % OPTIONAL: also equalize variance if desired
    % aN = aN / std(aN,'omitnan'); bN = bN / std(bN,'omitnan');

    % Resample both to the same length (shorter of the two)
    nA = numel(aN); nB = numel(bN); nC = min(nA,nB);
    xa = linspace(1,nA,nC);
    xb = linspace(1,nB,nC);
    aR = interp1(1:nA, aN, xa, 'linear');
    bR = interp1(1:nB, bN, xb, 'linear');
    x  = 1:nC;  % volumes; replace with time if you have TR, e.g., x=(0:nC-1)*TR;

    % Decide which trace is Ghrelin vs Placebo for this subject
    if gh_flag == 1
        % Phases 1–5 = Ghrelin; 6–10 = Placebo
        pG = plot(x, aR, 'LineWidth', 1.5, 'Color', cG);    % Ghrelin
        pP = plot(x, bR, 'LineWidth', 1.5, 'Color', cP);    % Placebo
        condLabel = sprintf('P%d (Ghrelin) vs P%d (Placebo)', k, k+5);
    else
        % Phases 1–5 = Placebo; 6–10 = Ghrelin
        pG = plot(x, bR, 'LineWidth', 1.5, 'Color', cG);    % Ghrelin
        pP = plot(x, aR, 'LineWidth', 1.5, 'Color', cP);    % Placebo
        condLabel = sprintf('P%d (Placebo) vs P%d (Ghrelin)', k, k+5);
    end

    yline(0, ':', 'Color', [0.7 0.7 0.7]);
    xlabel('Time (volumes)'); ylabel('BOLD (mean-normalized)');
    title(condLabel);

    if mark_peaks
        % Mark peaks for each condition trace with robust prominence
        promG = max(1.25*mad(pG.YData,1), 0);  % add a floor if needed
        promP = max(1.25*mad(pP.YData,1), 0);
        [pksG,locG] = findpeaks(pG.YData, 'MinPeakProminence', promG, ...
                                'MinPeakDistance', max(1,round(nC*0.01)));
        [pksP,locP] = findpeaks(pP.YData, 'MinPeakProminence', promP, ...
                                'MinPeakDistance', max(1,round(nC*0.01)));
        plot(x(locG), pksG, 'o', 'MarkerFaceColor', cG, 'MarkerEdgeColor', 'w', 'MarkerSize', 4);
        plot(x(locP), pksP, 's', 'MarkerFaceColor', cP, 'MarkerEdgeColor', 'w', 'MarkerSize', 4);
    end

    box on; hold off;
end

% Legend (put it on the first tile)
axes(tl.Children(end)); % first tile axes
legend({'Ghrelin','Placebo'}, 'Location','best');


% PULSE DETECTION 
% -------------------------------------------------------------
%% ----- Parameters for pulse detection -----
TR = 2.8;                    % <-- set your TR in seconds
%refractory_sec = 5;         % minimum time between pulses
% minPeakDist = max(1, round(refractory_sec / TR));

AMP_THR = 0.25;           % signal change 
prom_factor = 0;          % scale on MAD for robustness
prom_floor  = 0.00;          % absolute floor on prominence 

%% ----- Overlay and pulse detection (1↔6 ... 5↔10), colored by condition -----
sub_idx   =6;  % choose subject
phase_ts  = S.All_Subs_Phases_NAcc_timeseries{sub_idx};   % cell {1..10}
%phase_ts  = All_Subs_Phases_HYPO_timeseries{sub_idx};   % cell {1..10}

% Get condition flag: 1 => phases 1-5 are Ghrelin; 0 => phases 1-5 are Placebo
gh_flag = cond_ghrelin{sub_idx, 'Ghrelin'};

% --- Colors ---
color_Placebo = [0.0, 0.0, 0.545];        % darkblue
color_Ghrelin = [0.72, 0.53, 0.04];       % darkgoldenrod
color_Pulse   = [0.8, 0, 0];              % red pulse marker

normalize_to_mean = @(x) (x - mean(x, 'omitnan'));

% --- Collect outputs ---
pulse_list_G = cell(5,1); pulse_count_G = zeros(5,1);
pulse_list_P = cell(5,1); pulse_count_P = zeros(5,1);

figure('Color','w','Name',sprintf('Subject %02d — Pulse detection example (NAcc)', sub_idx));
tl = tiledlayout(5,1,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('Subject %02d — Mean-normalized overlays with pulse detection', sub_idx));

for k = 1:5
    a = phase_ts{k};       % Phase 1–5 (Session A)
    b = phase_ts{k+5};     % Phase 6–10 (Session B)
    nexttile; hold on;

    if isempty(a) || isempty(b)
        axis off; title(sprintf('Phase %d vs %d (no data)', k, k+5)); continue;
    end

    % --- Normalize each run (mean-center) ---
    aN = normalize_to_mean(a(:));
    bN = normalize_to_mean(b(:));
    if k == 1
        IDP = 'Task'
    elseif k == 2
        IDP = 'Baseline'
    else
        IDP = 'Task Free'
    end
    % --- Determine Ghrelin vs Placebo assignment ---
    if gh_flag == 1
        sigG = aN; sigP = bN;
        lbl = strcat('Ghrelin vs. Saline - ', IDP);
    else
        sigG = bN; sigP = aN;
        lbl = strcat('Saline vs. Ghrelin - ', IDP);
    end

    % --- Time vector in TR units (no interpolation!) ---
    nTR_G = numel(sigG);
    nTR_P = numel(sigP);
    tG = (0:nTR_G-1);    % TR indices for Ghrelin
    tP = (0:nTR_P-1);    % TR indices for Placebo

    % --- Plot traces ---
    pG = plot(tG, sigG, 'LineWidth', 2.2, 'Color', color_Ghrelin);
    pP = plot(tP, sigP, 'LineWidth', 2.2, 'Color', color_Placebo);
    yline(0, ':', 'Color', [0.7 0.7 0.7]);
    xlabel('Time (TRs)'); ylabel('BOLD (a.u.)');
    title(lbl); box on;

    % --- Pulse detection (no interpolation) ---
    [pksG, locG] = findpeaks(sigG, 'MinPeakHeight', AMP_THR);
    [pksP, locP] = findpeaks(sigP, 'MinPeakHeight', AMP_THR);

    % --- Mark detected pulses ---
    plot(locG-1, pksG, 'x', 'Color', color_Pulse, 'MarkerSize', 8, 'LineWidth', 2.5);
    plot(locP-1, pksP, 'x', 'Color', color_Pulse, 'MarkerSize', 8, 'LineWidth', 2.5);

    % --- Save results ---
    pulse_list_G{k} = locG(:);
    pulse_list_P{k} = locP(:);
    pulse_count_G(k) = numel(locG);
    pulse_count_P(k) = numel(locP);

    legend([pG pP], {'Ghrelin', 'Placebo'}, 'Location','best');
    hold off;
end

    


% ----- Save figure in fullscreen -----
set(gcf, 'Units', 'normalized', 'Position', [0 0 0.5 0.5]);  % fullscreen figure
saveas(gcf, fullfile('/mnt/TUE_data/TUE008/Data/TUE008_NIMG/Study/derivatives/conn/finalsampleCONN/Pulse_Detection/', sprintf('Sub%02d_PulseDetection_ExampleNAcc_Thr34.png', sub_idx)));

disp(table((1:5).', pulse_count_G, pulse_count_P, ...
     'VariableNames', {'Pair(1-6..5-10)','GhrelinPulses','PlaceboPulses'}));

 %% PULSE DETECTION FOR ALL SUBJS AND ROIS 

 % ---------------- Detection params ----------------
 AMP_THR       = 0.25;    % %PSC amplitude threshold
MAD_K         = 0;   % if >0, use MinPeakProminence = MAD_K * MAD(trace)
MIN_DIST_FRAC = 0;   % min peak distance as fraction of run length

own_detection = 0; % 0 if use paper code 

% ---------------- Load time series ----------------
S = load("All_Subs_Phases_timeseries.mat"); 
fieldnames(S)

% Wrap into a struct for iteration
ROIs = struct( ...
   'name', {'Hypothalamus','NAcc','VTA','Insula','Caudate','Putamen'}, ...
   'var',  {'All_Subs_Phases_HYPO_timeseries','All_Subs_Phases_NAcc_timeseries','All_Subs_Phases_VTA_timeseries', ...
            'All_Subs_Phases_INS_timeseries','All_Subs_Phases_CAUD_timeseries','All_Subs_Phases_PUT_timeseries'} );


% ---------------- Helpers ----------------
psc = @(x) 100 * (x(:) ./ median(x(:),'omitnan') - 1);  % median-based PSC

% ---------------- Accumulator for result table ----------------
rows = {};
row_i = 0;

% ---------------- Main: loop subjects & ROIs ----------------
for sub = 1:26
    gh_flag = logical(cond_ghrelin{sub, 3});  % 1 => phases 1–5 = Ghrelin
    sesA = 1:5; 
    sesB = 6:10; %later decide ghre or plac

    for r = 1:numel(ROIs)
        ts_cell = S.(ROIs(r).var){sub};  % expected: 1x10 or 10x1 cell of vectors
        assert(numel(ts_cell) >= phases, 'Subject %d %s has <10 phases', sub, ROIs(r).name);

        pulses_per_phase = zeros(phases,1);
        vols_per_phase   = zeros(phases,1);
        pulse_idx_per_phase = cell(phases,1);   

        for pha = 1:phases
            x = ts_cell{pha};
            if isempty(x)
                pulses_per_phase(pha) = 0; 
                vols_per_phase(pha)   = 0; 
                pulse_idx_per_phase{pha} = []; 
                continue;
            end
            x = x(:);
            vols_per_phase(pha) = numel(x);

            % Convert to PSC and detect peaks
            x_psc = psc(x);

            % Detect Pulses with certain threshold 
            [~,locs] = findpeaks(x_psc, 'MinPeakHeight', AMP_THR);

            
            pulses_per_phase(pha)   = numel(locs);
            pulse_idx_per_phase{pha} = locs(:)';   % save WHEN pulses occur (volume indices)
        end

        % Map phases to conditions (per subject)
        if gh_flag
            gh_pulses = pulses_per_phase(sesA);
            pl_pulses = pulses_per_phase(sesB);
            gh_idx    = pulse_idx_per_phase(sesA);  
            pl_idx    = pulse_idx_per_phase(sesB);  
        else
            gh_pulses = pulses_per_phase(sesB);
            pl_pulses = pulses_per_phase(sesA);
            gh_idx    = pulse_idx_per_phase(sesB);  
            pl_idx    = pulse_idx_per_phase(sesA);
        end

        % Append rows: Ghrelin, Placebo
        row_i = row_i + 1;
        rows(row_i,:) = {sub, ROIs(r).name, 'Ghrelin', gh_pulses', gh_idx};   
        row_i = row_i + 1;
        rows(row_i,:) = {sub, ROIs(r).name, 'Placebo', pl_pulses', pl_idx};  
    end
end

PulseTable = cell2table(rows, 'VariableNames', ...
    {'Subject','ROI','Condition','PulseCount','PulseIdx'});


% ---------------- Save ----------------
out_csv = fullfile(Pulse_output, 'PulseCounts_byCondition_ROIs_fromSavedTS.csv');
PulseTableForCsv = removevars(PulseTable, 'PulseIdx');   
writetable(PulseTableForCsv, out_csv);
save(fullfile(Pulse_output, 'PulseCounts_withLocs.mat'), 'PulseTable'); 



%% --- Overlaid histograms per phase (Ghrelin vs Saline), pooled across subjects ---
cP = [0.0, 0.0, 0.545];        % darkblue
cG = [0.72, 0.53, 0.04];       % darkgoldenrod
phases_in_session = 5;          % your PulseCount vectors are 1x5 per row

roi_plots = ["Hypothalamus","NAcc","VTA","Insula","Caudate","Putamen"]; 

for r = 1:size(roi_plots,2)
    roi_to_plot = roi_plots(r) ;   

    % Filters
    selROI = strcmp(PulseTable.ROI, roi_to_plot);
    isG    = strcmp(PulseTable.Condition, 'Ghrelin');
    isP    = strcmp(PulseTable.Condition, 'Placebo');

    % Helper to extract kth element from each cell row (skips rows with shorter vectors)
    getK = @(cells,k) cellfun(@(v) (k<=numel(v)) * v(k) + (k>numel(v)) * NaN, cells);

    fig = figure('Color','w','Name',sprintf('Pulse count histograms — %s', roi_to_plot));
    set(fig,'Units','normalized','OuterPosition',[0 0 1 1]);
    t = tiledlayout(3,2,'TileSpacing','compact','Padding','compact'); % 6 tiles; we’ll use first 5
    title(t, sprintf('ROI: %s — Pulse count histograms', roi_to_plot), ...
      'FontWeight','bold');

    for k = 1:phases_in_session
        nexttile; hold on;

        gvals = PulseTable.PulseCount(selROI & isG,k);
        pvals =  PulseTable.PulseCount(selROI & isP,k);

        % Common integer bin edges for fair comparison
        pooled = [gvals; pvals];
        edges = (floor(min(pooled))-1.5) : 1 : (ceil(max(pooled))+1.5);

       % histogram(gvals, (edges), 'FaceColor', cG, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        h1 = histfit(gvals, 12,'poisson');
        h1(1).FaceColor = cG;
        h1(2).Color =  cG;
        set(h1(1), 'FaceAlpha', 0.5);
        hold on
        h2 = histfit(pvals, 20,'poisson');
        h2(1).FaceColor = cP;
        h2(2).Color =  cP;
        set(h2(1), 'FaceAlpha', 0.5);
        %histogram(pvals, edges, 'FaceColor', cP, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        if k == 1
            IDP = 'Task'
        elseif k == 2
            IDP = 'Baseline'
        elseif k ==3 
            IDP = 'Task Free 1'
        elseif k == 4 
            IDP = 'Task Free 2' 
        else
            IDP = 'Task Free 3'
        end
         title(strcat('Phase - ', IDP));
        xlabel('Pulse count'); ylabel('Subjects');
        legend({'Ghrelin','Ghrelin fit','Saline','Saline fit'}, 'Location','best');
        box on; hold off;
    end
    
    set(gcf, 'Units', 'normalized', 'Position', [0 0 1 1]);  % fullscreen figure
    saveas(gcf, fullfile('/mnt/TUE_data/TUE008/Data/TUE008_NIMG/Study/derivatives/conn/finalsampleCONN/Pulse_Detection/', sprintf('PulseHist_028_%s.png', char(roi_to_plot))));


end 
close all 

%% test 

rows = [];
for i = 1:height(PulseTable)
    pc = PulseTable.PulseCount(i,1:5);            % 1x5 counts
    n  = numel(pc);
    rows = [rows; table( ...
        repmat(PulseTable.Subject(i), n, 1), ...
        repmat(string(PulseTable.ROI(i)), n, 1), ...
        repmat(string(PulseTable.Condition(i)), n, 1), ...
        (1:n)', pc(:), ...
        'VariableNames', {'Subject','ROI','Condition','Phase','PulseCount'})];
end

PulseLong = rows;
PulseLong.Subject   = categorical(PulseLong.Subject);
PulseLong.ROI       = categorical(PulseLong.ROI);
PulseLong.Condition = categorical(PulseLong.Condition, {'Placebo','Ghrelin'}); % set baseline if desired
PulseLong.Phase     = categorical(PulseLong.Phase);   % treat phases as discrete


%% ---- GLME: Poisson (counts) with log link ----


% Remove baseline for now 

% T = PulseLong(PulseLong.Phase ~= '1', :);
% T = PulseLong(PulseLong.Phase == '1' |PulseLong.Phase == '2' , :);

T = PulseLong;


pnum = double(T.Phase);   % already numeric

% Map to 3 levels: 0 (Phase 2), 1 (Phase 1), 2 (Phases 3–5)
BaseCond = NaN(size(pnum));
BaseCond(pnum == 2) = 0;
BaseCond(pnum == 1) = 2;
BaseCond(ismember(pnum, [3 4 5])) = 1;

% Add to table (numeric or categorical—pick one)
T.BaseCond = BaseCond;  % numeric 0/1/2
T.BaseCond = categorical(T.BaseCond, [0 1 2], {'Baseline','TaskFree','Task'});
T.BaseCond = reordercats(T.BaseCond, {'Baseline','TaskFree','Task'});

PulseLong_NAcc = T(T.ROI == 'NAcc', :);
PulseLong_HYPO = T(T.ROI == 'Hypothalamus', :);
PulseLong_VTA = T(T.ROI == 'VTA', :);
PulseLong_INS = T(T.ROI == 'Insula', :);
PulseLong_C = T(T.ROI == 'Caudate', :);
PulseLong_P = T(T.ROI == 'Putamen', :);

PulseLong_NAcc.Phase     = removecats(categorical(PulseLong_NAcc.Phase));  % drops the now-missing '2'
PulseLong_VTA.Phase     = removecats(categorical(PulseLong_VTA.Phase));  % drops the now-missing '2'
PulseLong_HYPO.Phase     = removecats(categorical(PulseLong_HYPO.Phase));  % drops the now-missing '2'
PulseLong_INS.Phase     = removecats(categorical(PulseLong_INS.Phase));  % drops the now-missing '2'
PulseLong_C.Phase     = removecats(categorical(PulseLong_C.Phase));  % drops the now-missing '2'
PulseLong_P.Phase     = removecats(categorical(PulseLong_P.Phase));  % drops the now-missing '2'

T = PulseLong_HYPO; 
T = PulseLong_VTA; 
T = PulseLong_NAcc; 
T = PulseLong_INS; 
T = PulseLong_C; 
T = PulseLong_P; 

ph = categories(T.Phase);

%ph = setdiff(ph, '3', 'stable');           % all others, keep original order
%T.Phase = reordercats(T.Phase, ['3'; ph]); % '3' becomes the baseline
T.BaseCond     = (categorical(PulseLong_P.BaseCond)); 


glme = fitglme(T, ...
   'PulseCount ~ Condition* BaseCond + (1 |Subject)', ...
   'Distribution','Poisson','Link','log','FitMethod','Laplace');

disp(glme)


T = PulseLong_NAcc; 
writetable(T,fullfile(strcat(Pulse_output, 'PulseCount_NAcc_Thresh28.xlsx')))  ; 

T_all = vertcat(PulseLong_HYPO, PulseLong_VTA, PulseLong_NAcc, PulseLong_INS, PulseLong_C,PulseLong_P); 
writetable(T_all,fullfile(strcat(Pulse_output, 'PulseCount_allROIs_Thresh28.xlsx')))  ; 

%% 

%% -- Verteilung Diff 

VT = PulseLong(PulseLong.Phase ~= '1' &  PulseLong.Phase ~= '2', :);

VT = PulseLong(PulseLong.Phase == '4',:);



PulseLong_NAcc = VT(VT.ROI == 'NAcc', :);
PulseLong_HYPO = VT(VT.ROI == 'Hypothalamus', :);
PulseLong_VTA = VT(VT.ROI == 'VTA', :);
PulseLong_Put = VT(VT.ROI == 'Putamen', :);

Dist_G = PulseLong_Put(PulseLong_Put.Condition == 'Ghrelin',:)
Dist_P = PulseLong_Put(PulseLong_Put.Condition == 'Placebo',:)

  [H,P] = kstest2(Dist_G.PulseCount, Dist_P.PulseCount)


  Dist_G = PulseLong_HYPO(PulseLong_HYPO.Condition == 'Ghrelin',:)
Dist_P = PulseLong_HYPO(PulseLong_HYPO.Condition == 'Placebo',:)

  [H,P] = kstest2(Dist_G.PulseCount, Dist_P.PulseCount)


  
Dist_G = PulseLong_VTA(PulseLong_VTA.Condition == 'Ghrelin',:)
Dist_P = PulseLong_VTA(PulseLong_VTA.Condition == 'Placebo',:)

  [H,P] = kstest2(Dist_G.PulseCount, Dist_P.PulseCount)


  VT = PulseLong_NAcc(PulseLong_NAcc.BaseCond ~= '0',:);

Dist_G = VT(VT.Condition == 'Ghrelin',:)
Dist_P = VT(VT.Condition == 'Placebo',:)

  [H,P] = kstest2(Dist_G.PulseCount, Dist_P.PulseCount)
  
  
%% ================== Co-occurrence analysis across ROI pairs ==================

% ---- Parameters ----
LAG_MAX_VOL = 2;   % "few volumes later" maximum lag (inclusive)

cooc_mat_path = fullfile(Pulse_output, 'PulseCounts_withLocs.mat');
tmp = load(cooc_mat_path);
PulseTable = tmp.PulseTable;

% ---- Pairs to compare (ROI1 -> ROI2) ----
Pairs = {
    'Hypothalamus','NAcc';
    'NAcc','VTA'; 
    'NAcc','Caudate';
    'Hypothalamus','Insula';
    'NAcc','Putamen';
    'Hypothalamus','Putamen';
    'Hypothalamus','Caudate'
};

subjects = unique(PulseTable.Subject);
conds    = {'Ghrelin','Placebo'};

% ---- Results accumulators ----
rows_sum   = {};   % summary (CSV-friendly)
rows_detail = {};  % detailed (kept in MAT)

for s = reshape(subjects,1,[])
    for c = 1:numel(conds)
        cond = conds{c};

        for p = 1:size(Pairs,1)
            ROI1 = Pairs{p,1};
            ROI2 = Pairs{p,2};

            % fetch rows for this subject/condition/ROI
            i1 = find(PulseTable.Subject==s & strcmp(PulseTable.Condition,cond) & strcmp(PulseTable.ROI,ROI1), 1);
            i2 = find(PulseTable.Subject==s & strcmp(PulseTable.Condition,cond) & strcmp(PulseTable.ROI,ROI2), 1);
            if isempty(i1) || isempty(i2)
                % no data → add NA rows for transparency (optional: skip)
                for ph = 1:5
                    rows_sum(end+1,:) = {s, cond, [ROI1 '-' ROI2], ph, LAG_MAX_VOL, NaN, NaN, NaN, NaN, NaN}; 
                    rows_detail(end+1,:) = {s, cond, [ROI1 '-' ROI2], ph, LAG_MAX_VOL, [], [], [], []};
                end
                
            end

            ts1 = PulseTable.PulseIdx{i1};   % 1x5 cell: indices per phase within condition
            ts2 = PulseTable.PulseIdx{i2};   % 1x5 cell

            nph = min(5, min(numel(ts1), numel(ts2)));
            for ph = 1:nph
                v1 = ts1{ph}; if isempty(v1), v1 = []; end
                v2 = ts2{ph}; if isempty(v2), v2 = []; end
                v1 = sort(v1(:)');  % ensure row vector & sorted
                v2 = sort(v2(:)');

                % match ROI1→ROI2 at same time or later within LAG_MAX_VOL
                [pairs12, lags12] = match_forward(v1, v2, LAG_MAX_VOL);
                n1 = numel(v1); n2 = numel(v2);
                nm = size(pairs12,1);
                nsame = sum(lags12==0);
                nlater = sum(lags12>0);

                % summary row (CSV-safe)
                rows_sum(end+1,:) = {s, cond, [ROI1 '-' ROI2], ph, LAG_MAX_VOL, n1, n2, nm, nsame, nlater}; %#ok<AGROW>

                % detail row (kept in MAT)
                second_times = pairs12(:,2)';  % indices in ROI2 where matches occurred
                rows_detail(end+1,:) = {s, cond, [ROI1 '-' ROI2], ph, LAG_MAX_VOL, pairs12, second_times, lags12(:)', v1}; %#ok<AGROW>
            end
        end
    end
end

% ---- Build tables ----
CoocSummary = cell2table(rows_sum, 'VariableNames', ...
    {'Subject','Condition','Pair','Phase','MaxLagVol', ...
     'ROI1_Peaks','ROI2_Peaks','Matches','SameTime','Later'});

CoocDetail = cell2table(rows_detail, 'VariableNames', ...
    {'Subject','Condition','Pair','Phase','MaxLagVol', ...
     'MatchedPairs_ROI1_ROI2','MatchedSecondIdx','MatchedLags', ...
     'ROI1_AllIdx'});

% STATS 
% Categorical variables with explicit reference level
CoocSummary.Subject   = categorical(CoocSummary.Subject);
CoocSummary.Pair      = categorical(CoocSummary.Pair);
CoocSummary.Phase     = categorical(CoocSummary.Phase);  % phases 1..5
CoocSummary.Condition = categorical(CoocSummary.Condition, {'Placebo','Ghrelin'}); % Placebo = ref

% Drop rows with zero "opportunities" (cannot compute log offset)

pairSel = ["Hypothalamus-NAcc","NAcc-VTA","Hypothalamus-Insula","Hypothalamus-Caudate","Hypothalamus-Putamen","NAcc-Caudate","NAcc-Putamen"]


%pstr = string(CoocSummary.Phase);
%pnum = str2double(pstr);                        
%CoocSummary.BaseCond = double(pnum ~= 2);

pnum = double(CoocSummary.Phase);   % already numeric

% Map to 3 levels: 0 (Phase 2), 1 (Phase 1), 2 (Phases 3–5)
BaseCond = NaN(size(pnum));
BaseCond(pnum == 2) = 0;
BaseCond(pnum == 1) = 2;
BaseCond(ismember(pnum, [3 4 5])) = 1;

% Add to table (numeric or categorical—pick one)
CoocSummary.BaseCond = BaseCond;  % numeric 0/1/2
CoocSummary.BaseCond = categorical(CoocSummary.BaseCond, [0 1 2], {'Baseline','TaskFree','Task'});
CoocSummary.BaseCond = reordercats(CoocSummary.BaseCond, {'Baseline','TaskFree','Task'});


for i= 1:size(pairSel,2)
    T = CoocSummary(CoocSummary.ROI1_Peaks > 0, :);
    idxPair = strcmpi(cellstr(T.Pair), pairSel(i));
    T = T(idxPair, :);
    
    pairSel(i)
    
    % T = T(T.Phase ~= '2', :);
    %T = T(T.Phase ~= '1', :);

    %T.Phase     = removecats(categorical(T.Phase));  % drops the now-missing '2'
    %T.BaseCond     = (categorical(T.BaseCond));  % drops the now-missing '2'

    %ph = categories(T.Phase);
    %ph = setdiff(ph, '3', 'stable');           % all others, keep original order
    %T.Phase = reordercats(T.Phase, ['3'; ph]); % '3' becomes the baseline


    % Poisson GLME for co-occurrence rate: Matches per ROI1 peak
    glme = fitglme(T, ...
       'Matches ~ Condition * BaseCond + (1 |Subject)', ...
       'Distribution','Poisson','Link','log','FitMethod','Laplace');

    disp(glme)
    
    writetable(T,fullfile(strcat(Pulse_output, 'PulseCooccurance_Clean_28_' ,pairSel{i}, '.xlsx')))  ; 
end 



%% Plotting 

% Pairs to plot (same list you used)
Pairs = {
    'Hypothalamus','NAcc';
    'Hypothalamus','Putamen';
};

subjects = unique(PulseTable.Subject);

% Colors
%cG = [0.10 0.60 0.20];   % Ghrelin
%cP = [0.20 0.45 0.85];   % Placebo
cM = [0.85 0.33 0.10];   % Matched markers
cLineG = cG;             % connector line (Ghrelin)
cLineP = cP;             % connector line (Placebo)
gridCol = [0.75 0.75 0.75];
baseTS = All_Subs_Phases_HYPO_timeseries;

for p = 1:size(Pairs,1)
    ROI1 = Pairs{p,1};
    ROI2 = Pairs{p,2};
    pairName = [ROI1 '-' ROI2];

    for s = reshape(subjects,1,[])
        % ---- subject-specific phase lengths & offsets ----
        ts_cell_len = baseTS{s};                  % {1..10} cell
        phase_len   = cellfun(@(v) numel(v), ts_cell_len(:));  % 10x1
        offsets = [0; cumsum(phase_len(1:9))];    % length 10, 0-based offsets
        totalLen = offsets(end) + phase_len(10);

        % ---- mapping phases→condition for this subject ----
        gh_flag = logical(cond_ghrelin{s,3});     % 1: phases 1–5 = Ghrelin
        if gh_flag
            gh_phases = 1:5; pl_phases = 6:10;
        else
            gh_phases = 6:10; pl_phases = 1:5;
        end

        % ---- fetch rows for ROI1 & ROI2 (both conditions) ----
        i1G = find(PulseTable.Subject==s & strcmp(PulseTable.Condition,'Ghrelin') & strcmp(PulseTable.ROI,ROI1), 1);
        i1P = find(PulseTable.Subject==s & strcmp(PulseTable.Condition,'Placebo') & strcmp(PulseTable.ROI,ROI1), 1);
        i2G = find(PulseTable.Subject==s & strcmp(PulseTable.Condition,'Ghrelin') & strcmp(PulseTable.ROI,ROI2), 1);
        i2P = find(PulseTable.Subject==s & strcmp(PulseTable.Condition,'Placebo') & strcmp(PulseTable.ROI,ROI2), 1);


        % ---- convert per-phase local indices → global x (concatenate phases) ----
        % helper to expand a 1x5 cell (within-condition) into global indices
        expand_global = @(cells, abs_phases) ...
            sort(cell2mat(arrayfun(@(k) offsets(abs_phases(k)) + cells{k}(:), 1:numel(abs_phases), 'uni', 0)));

        % --- ROI1 (Ghrelin) ---
        tmp = arrayfun(@(kk) offsets(gh_phases(kk)) + double(PulseTable.PulseIdx{i1G}{kk}(:)), ...
                       1:numel(gh_phases), 'uni', 0);
        g1 = sort(vertcat(tmp{:}));

        % --- ROI1 (Placebo) ---
        tmp = arrayfun(@(kk) offsets(pl_phases(kk)) + double(PulseTable.PulseIdx{i1P}{kk}(:)), ...
                       1:numel(pl_phases), 'uni', 0);
        p1 = sort(vertcat(tmp{:}));

        % --- ROI2 (Ghrelin) ---
        tmp = arrayfun(@(kk) offsets(gh_phases(kk)) + double(PulseTable.PulseIdx{i2G}{kk}(:)), ...
                       1:numel(gh_phases), 'uni', 0);
        g2 = sort(vertcat(tmp{:}));

        % --- ROI2 (Placebo) ---
        tmp = arrayfun(@(kk) offsets(pl_phases(kk)) + double(PulseTable.PulseIdx{i2P}{kk}(:)), ...
                       1:numel(pl_phases), 'uni', 0);
        p2 = sort(vertcat(tmp{:}));

        % ---- match co-occurrences within each condition (across its 5 phases) ----
        [pairsG, lagsG] = match_forward(g1, g2, LAG_MAX_VOL);
        [pairsP, lagsP] = match_forward(p1, p2, LAG_MAX_VOL);

        % split matched/unmatched
        m1G = pairsG(:,1)'; m2G = pairsG(:,2)'; u1G = setdiff(g1, m1G,'stable'); u2G = setdiff(g2, m2G,'stable');
        m1P = pairsP(:,1)'; m2P = pairsP(:,2)'; u1P = setdiff(p1, m1P,'stable'); u2P = setdiff(p2, m2P,'stable');

        % ---- figure (single axes; all phases appended; both conditions colored) ----
        fig = figure('Color','w','Name', sprintf('S%02d — %s — All phases appended', s, pairName));
        set(fig,'Units','normalized','OuterPosition',[0 0 1 1])
        
        ax = axes(fig); hold(ax,'on');

        % background: light bands for each absolute phase (optional)
        yl = [0.5 2.5];
        for ph = 1:10
            if mod(ph,2)==0
                patch([offsets(ph) offsets(ph)+phase_len(ph) offsets(ph)+phase_len(ph) offsets(ph)], ...
                      [yl(1) yl(1) yl(2) yl(2)], [0.95 0.95 0.95], 'EdgeColor','none','Parent',ax);
            end
        end
        uistack(findobj(ax,'Type','patch'),'bottom');

        % Ghrelin pulses (unmatched faint, matched highlighted)
        if ~isempty(u1G), scatter(u1G, repmat(2,1,numel(u1G)), 18, cG, 'filled', 'MarkerFaceAlpha',0.25, 'MarkerEdgeAlpha',0.25); end
        if ~isempty(u2G), scatter(u2G, repmat(1,1,numel(u2G)), 18, cG, 'filled', 'MarkerFaceAlpha',0.25, 'MarkerEdgeAlpha',0.25); end
        if ~isempty(m1G), scatter(m1G, repmat(2,1,numel(m1G)), 34, cG, 'filled'); end
        if ~isempty(m2G), scatter(m2G, repmat(1,1,numel(m2G)), 34, cG, 'filled'); end
        for k = 1:numel(m1G)
            plot([m1G(k) m2G(k)], [2 1], '--', 'Color', cLineG, 'LineWidth', 2);
        end

        % Placebo pulses
        if ~isempty(u1P), scatter(u1P, repmat(2,1,numel(u1P)), 18, cP, 'filled', 'MarkerFaceAlpha',0.25, 'MarkerEdgeAlpha',0.25); end
        if ~isempty(u2P), scatter(u2P, repmat(1,1,numel(u2P)), 18, cP, 'filled', 'MarkerFaceAlpha',0.25, 'MarkerEdgeAlpha',0.25); end
        if ~isempty(m1P), scatter(m1P, repmat(2,1,numel(m1P)), 34, cP, 'filled'); end
        if ~isempty(m2P), scatter(m2P, repmat(1,1,numel(m2P)), 34, cP, 'filled'); end
        for k = 1:numel(m1P)
            plot([m1P(k) m2P(k)], [2 1], '--', 'Color', cLineP, 'LineWidth', 2);
        end

        % condition blocks (optional labels)
        xG = [offsets(gh_phases(1)), offsets(gh_phases(end))+phase_len(gh_phases(end))];
        xP = [offsets(pl_phases(1)), offsets(pl_phases(end))+phase_len(pl_phases(end))];
        text(mean(xG), 2.35, 'Ghrelin', 'HorizontalAlignment','center', 'Color', cG, 'FontWeight','bold');
        text(mean(xP), 2.35, 'Placebo', 'HorizontalAlignment','center', 'Color', cP, 'FontWeight','bold');

        % axes cosmetics
        xlim([0, totalLen+1]);
        yticks([1 2]); yticklabels({ROI2, ROI1});
        ylim([0.5 2.5]);
        xlabel('Volume (all phases appended)'); ylabel('ROI');
        title(sprintf('Example Subject %02d — %s — Max lag = %d vols  |  G: (matches %d), P: (matches %d)', ...
               s, pairName, LAG_MAX_VOL,  size(pairsG,1), size(pairsP,1)));

        % phase boundaries
        for ph = 2:10, xline(offsets(ph), ':', 'Color', gridCol); end
        box on; hold off;

        % save
        safePair = strrep(pairName,'-','_');
        out_png = fullfile(Pulse_output, sprintf('Example_Cooc_Appended_S%02d_%s_AllPhases.png', s, safePair));
        exportgraphics(fig, out_png, 'Resolution', 300, 'BackgroundColor','white');
        fprintf('Saved: %s\n', out_png);
        
        close all 
    end
end


%% 

np = numel(pairSel);

% Colors
cG = [0.10 0.60 0.20];  % Ghrelin
cP = [0.20 0.45 0.85];  % Placebo

% Bins for discrete lags (0..LAG_MAX_VOL)
edges = (-0.5) : 1 : (LAG_MAX_VOL + 0.5);

% Big figure
fig = figure('Color','w','Name','Co-occurrence lag densities (all subjects)');
set(fig,'Units','normalized','OuterPosition',[0 0 1 1]); 

ncols = 2; nrows = ceil(np / ncols);
tl = tiledlayout(fig, nrows, ncols, 'TileSpacing','compact', 'Padding','compact');
title(tl, sprintf('Co-occurrence Lag Densities (0..%d vols) — All Subjects', LAG_MAX_VOL), 'FontWeight','bold');

for i = 1:np
    pairName = pairSel(i);
    nexttile; hold on;

    % Pull lags across all subjects & phases for this pair & condition
    idxG = strcmpi(CoocDetail.Condition, 'Ghrelin')  & strcmpi(CoocDetail.Pair, pairName);
    idxP = strcmpi(CoocDetail.Condition, 'Placebo') & strcmpi(CoocDetail.Pair, pairName);

        % ---- collect lags for this pair/condition (robust to row/col + empties) ----
        lagsG = [];
        if any(idxG)
            Cg = CoocDetail.MatchedLags(idxG);                 % cell array
            Cg = Cg(~cellfun('isempty', Cg));                  % drop empties (optional)
            Cg = cellfun(@(v) v(:), Cg, 'uni', 0);             % force column vectors
            if ~isempty(Cg), lagsG = vertcat(Cg{:}); end       % stack into one column
        end

        lagsP = [];
        if any(idxP)
            Cp = CoocDetail.MatchedLags(idxP);
            Cp = Cp(~cellfun('isempty', Cp));
            Cp = cellfun(@(v) v(:), Cp, 'uni', 0);
            if ~isempty(Cp), lagsP = vertcat(Cp{:}); end
        end
    % Histograms normalized to probability
    if ~isempty(lagsG)
        histogram(lagsG, edges, 'Normalization','probability', ...
                  'FaceColor', cG, 'EdgeColor', 'none', 'FaceAlpha', 0.45);
    end
    if ~isempty(lagsP)
        histogram(lagsP, edges, 'Normalization','probability', ...
                  'FaceColor', cP, 'EdgeColor', 'none', 'FaceAlpha', 0.45);
    end

    % Rug marks (optional): tiny ticks at each observed lag
    if ~isempty(lagsG), scatter(lagsG + (rand(size(lagsG))-0.5)*0.04, -0.01*ones(size(lagsG)), 6, cG, 'filled'); end
    if ~isempty(lagsP), scatter(lagsP + (rand(size(lagsP))-0.5)*0.04, -0.02*ones(size(lagsP)), 6, cP, 'filled'); end

    % Cosmetics
    xlim([-0.5, LAG_MAX_VOL+0.5]); xticks(0:LAG_MAX_VOL);
    ylim([0, max(0.001, ylim(gca)*[0;1])]); % keep y>=0
    xlabel('Lag (volumes)'); ylabel('Probability');
    nG = numel(lagsG); nP = numel(lagsP);
    nsameG = sum(lagsG==0); nsameP = sum(lagsP==0);
    box on; hold off;
end

% Put a single legend in the last tile (or create a layout legend)
lg = legend(tl.Children(1), {'Ghrelin','Placebo'}, 'Location','southoutside', 'Orientation','horizontal');
lg.Layout.Tile = 'south';

% Save
out_png = fullfile(Pulse_output, sprintf('Cooc_LagDensities_AllSubjects_%dvol.png', LAG_MAX_VOL));
exportgraphics(fig, out_png, 'Resolution', 300, 'BackgroundColor','white');
fprintf('Saved: %s\n', out_png);

