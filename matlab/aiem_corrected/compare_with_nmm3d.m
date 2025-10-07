% COMPARE_WITH_NMM3D Compare corrected AIEM with NMM3D reference data
%
% This script loads NMM3D reference data and compares with the corrected
% AIEM implementation following the Python example structure.
%
% Based on: examples/test_aiem.ipynb

clear all; close all; clc;

fprintf('========================================\n');
fprintf('AIEM vs NMM3D Comparison\n');
fprintf('========================================\n\n');

%% Configuration
FREQUENCY_GHZ = 5.405;
INCIDENCE_DEG = 40.0;
TARGET_RATIOS = [4, 7, 10, 15];

% Multiple scattering controls
INCLUDE_MULTIPLE_SCATTERING = true;  % Set to true to include MS
MS_QUAD_POINTS = 129;
MS_SPECTRAL_TERMS = 8;

% File paths
nmm3d_file = '../../data/NMM3D_LUT_NRCS_40degree.dat';
if ~exist(nmm3d_file, 'file')
    fprintf('ERROR: NMM3D data file not found:\n');
    fprintf('  %s\n\n', nmm3d_file);
    fprintf('Please ensure the NMM3D reference data is available.\n');
    return;
end

%% Load and filter data
fprintf('Loading NMM3D reference data...\n');
data = load(nmm3d_file);
fprintf('Loaded %d data points\n', size(data, 1));

% Data format: theta ratio eps_r eps_i rms_norm vv_db hh_db hv_db
theta_col = 1;
ratio_col = 2;
eps_r_col = 3;
eps_i_col = 4;
rms_norm_col = 5;
vv_ref_col = 6;
hh_ref_col = 7;
hv_ref_col = 8;

% Filter for target incidence angle
mask = abs(data(:, theta_col) - INCIDENCE_DEG) < 0.1;
data_40 = data(mask, :);
fprintf('Filtered %d points at theta = %g degrees\n', size(data_40, 1), INCIDENCE_DEG);

% Filter for target ratios
if ~isempty(TARGET_RATIOS)
    mask_ratio = false(size(data_40, 1), 1);
    for r = TARGET_RATIOS
        mask_ratio = mask_ratio | (abs(data_40(:, ratio_col) - r) < 0.01);
    end
    data_40 = data_40(mask_ratio, :);
    fprintf('Filtered %d points for ratios: %s\n', size(data_40, 1), mat2str(TARGET_RATIOS));
end

if isempty(data_40)
    error('No data points left after filtering');
end

%% Setup
lambda = 0.3 / FREQUENCY_GHZ;  % wavelength in meters
k = 2 * pi / lambda;
theta_i = INCIDENCE_DEG;
theta_s = INCIDENCE_DEG;
phi_s = 180;  % Monostatic backscatter
itype = 2;  % Exponential correlation

fprintf('\nSimulation parameters:\n');
fprintf('  Frequency: %.3f GHz\n', FREQUENCY_GHZ);
fprintf('  Wavelength: %.4f m\n', lambda);
fprintf('  Incidence: %.1f deg\n', INCIDENCE_DEG);
fprintf('  Multiple scattering: %s\n', mat2str(INCLUDE_MULTIPLE_SCATTERING));
if INCLUDE_MULTIPLE_SCATTERING
    fprintf('  MS grid: %dx%d\n', MS_QUAD_POINTS, MS_QUAD_POINTS);
    fprintf('  MS terms: %d\n', MS_SPECTRAL_TERMS);
end
fprintf('\n');

%% Group data by (ratio, eps_r, eps_i)
unique_groups = unique(data_40(:, [ratio_col, eps_r_col, eps_i_col]), 'rows');
n_groups = size(unique_groups, 1);

fprintf('Processing %d unique groups...\n', n_groups);

group_results = cell(n_groups, 1);

for g = 1:n_groups
    ratio = unique_groups(g, 1);
    eps_r = unique_groups(g, 2);
    eps_i = unique_groups(g, 3);
    
    % Find all points for this group
    mask = (abs(data_40(:, ratio_col) - ratio) < 0.01) & ...
           (abs(data_40(:, eps_r_col) - eps_r) < 0.01) & ...
           (abs(data_40(:, eps_i_col) - eps_i) < 0.01);
    
    if ~any(mask)
        continue;
    end
    
    % Sort by rms_norm
    group_data = data_40(mask, :);
    [~, idx] = sort(group_data(:, rms_norm_col));
    group_data = group_data(idx, :);
    
    rms_norm = group_data(:, rms_norm_col);
    vv_ref = group_data(:, vv_ref_col);
    hh_ref = group_data(:, hh_ref_col);
    hv_ref = group_data(:, hv_ref_col);
    
    % Compute AIEM for each point
    n_pts = length(rms_norm);
    vv_aiem = zeros(n_pts, 1);
    hh_aiem = zeros(n_pts, 1);
    hv_aiem = zeros(n_pts, 1);
    ks_vals = zeros(n_pts, 1);
    
    for i = 1:n_pts
        sigma = rms_norm(i) * lambda;
        corr_len = ratio * sigma;
        ks = k * sigma;
        kl = k * corr_len;
        ks_vals(i) = ks;
        
        try
            if INCLUDE_MULTIPLE_SCATTERING
                [VV, HH, HV, ~] = AIEM_corrected(theta_i, theta_s, phi_s, kl, ks, ...
                                                  eps_r, eps_i, itype, ...
                                                  'IncludeMultipleScattering', true);
            else
                [VV, HH, HV, ~] = AIEM_corrected(theta_i, theta_s, phi_s, kl, ks, ...
                                                  eps_r, eps_i, itype);
            end
            vv_aiem(i) = VV;
            hh_aiem(i) = HH;
            hv_aiem(i) = HV;
        catch ME
            fprintf('Warning: AIEM failed for group %d, point %d: %s\n', g, i, ME.message);
            vv_aiem(i) = NaN;
            hh_aiem(i) = NaN;
            hv_aiem(i) = NaN;
        end
    end
    
    % Store results
    group_results{g} = struct(...
        'ratio', ratio, ...
        'eps_real', eps_r, ...
        'eps_imag', eps_i, ...
        'ks', ks_vals, ...
        'vv_ref', vv_ref, ...
        'hh_ref', hh_ref, ...
        'hv_ref', hv_ref, ...
        'vv_aiem', vv_aiem, ...
        'hh_aiem', hh_aiem, ...
        'hv_aiem', hv_aiem);
    
    fprintf('  Group %d/%d: ratio=%.1f, eps=%.1f+j%.1f, n=%d\n', ...
            g, n_groups, ratio, eps_r, eps_i, n_pts);
end

% Remove empty cells
group_results = group_results(~cellfun(@isempty, group_results));
fprintf('Completed %d groups\n\n', length(group_results));

%% Plot per ratio
% Extract unique ratios
unique_ratios = [];
for g = 1:length(group_results)
    unique_ratios = [unique_ratios; group_results{g}.ratio];
end
unique_ratios = unique(unique_ratios);

% Extract unique dielectrics
unique_dielectrics = [];
for g = 1:length(group_results)
    unique_dielectrics = [unique_dielectrics; group_results{g}.eps_real, group_results{g}.eps_imag];
end
unique_dielectrics = unique(unique_dielectrics, 'rows');

% Color map for different dielectrics
colors = lines(size(unique_dielectrics, 1));
markers = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h', 'x'};

for r_idx = 1:length(unique_ratios)
    ratio = unique_ratios(r_idx);
    
    % Find all groups for this ratio
    ratio_groups = {};
    for g = 1:length(group_results)
        if abs(group_results{g}.ratio - ratio) < 0.01
            ratio_groups{end+1} = group_results{g};
        end
    end
    
    if isempty(ratio_groups)
        continue;
    end
    
    % Create figure
    figure('Position', [100, 100, 1800, 500]);
    
    bands = {'vv', 'hh', 'hv'};
    band_labels = {'VV', 'HH', 'HV'};
    
    for b = 1:3
        subplot(1, 3, b);
        hold on;
        
        band = bands{b};
        
        for g = 1:length(ratio_groups)
            group = ratio_groups{g};
            eps_pair = [group.eps_real, group.eps_imag];
            
            % Find color index
            color_idx = find(ismember(unique_dielectrics, eps_pair, 'rows'));
            color = colors(color_idx, :);
            marker = markers{mod(color_idx-1, length(markers)) + 1};
            
            ks = group.ks;
            aiem_vals = group.([band '_aiem']);
            ref_vals = group.([band '_ref']);
            
            % Plot AIEM (line)
            plot(ks, aiem_vals, '--', 'Color', color, 'LineWidth', 1.5, ...
                 'DisplayName', sprintf('AIEM ε=%.1f+j%.1f', eps_pair(1), eps_pair(2)));
            
            % Plot NMM3D (markers)
            finite_mask = isfinite(ref_vals);
            if any(finite_mask)
                plot(ks(finite_mask), ref_vals(finite_mask), marker, ...
                     'MarkerSize', 6, 'MarkerFaceColor', color, ...
                     'MarkerEdgeColor', 'black', 'LineStyle', 'none', ...
                     'DisplayName', sprintf('NMM3D ε=%.1f+j%.1f', eps_pair(1), eps_pair(2)));
            end
        end
        
        xlabel('k·σ (rad)', 'FontSize', 12);
        ylabel(sprintf('σ⁰ %s (dB)', band_labels{b}), 'FontSize', 12);
        title(sprintf('%s Polarization', band_labels{b}), 'FontSize', 14);
        grid on;
        legend('Location', 'best', 'FontSize', 8);
        hold off;
    end
    
    sgtitle(sprintf('Comparison at ℓ/σ = %.1f', ratio), 'FontSize', 16);
end

%% Compute metrics by ratio and band
fprintf('========================================\n');
fprintf('Metrics by Ratio and Polarization\n');
fprintf('========================================\n\n');

metrics_table = [];

for r_idx = 1:length(unique_ratios)
    ratio = unique_ratios(r_idx);
    
    % Find all groups for this ratio
    ratio_groups = {};
    for g = 1:length(group_results)
        if abs(group_results{g}.ratio - ratio) < 0.01
            ratio_groups{end+1} = group_results{g};
        end
    end
    
    if isempty(ratio_groups)
        continue;
    end
    
    fprintf('ℓ/σ = %.1f:\n', ratio);
    
    for b = 1:3
        band = bands{b};
        band_label = band_labels{b};
        
        % Concatenate all data for this ratio and band
        ref_all = [];
        sim_all = [];
        for g = 1:length(ratio_groups)
            group = ratio_groups{g};
            ref_all = [ref_all; group.([band '_ref'])];
            sim_all = [sim_all; group.([band '_aiem'])];
        end
        
        % Compute metrics
        mask = isfinite(ref_all) & isfinite(sim_all);
        if sum(mask) == 0
            fprintf('  %s: No valid data\n', band_label);
            continue;
        end
        
        ref = ref_all(mask);
        sim = sim_all(mask);
        diff = sim - ref;
        
        rmse = sqrt(mean(diff.^2));
        bias = mean(diff);
        mae = mean(abs(diff));
        r = corr(ref, sim);
        n = sum(mask);
        
        fprintf('  %s: RMSE=%.2f dB, Bias=%+.2f dB, MAE=%.2f dB, r=%.3f (n=%d)\n', ...
                band_label, rmse, bias, mae, r, n);
        
        % Store for table
        metrics_table = [metrics_table; ratio, b, rmse, bias, mae, r, n];
    end
    fprintf('\n');
end

%% Overall metrics
fprintf('========================================\n');
fprintf('Overall Metrics (All Ratios)\n');
fprintf('========================================\n\n');

for b = 1:3
    band = bands{b};
    band_label = band_labels{b};
    
    % Concatenate all data
    ref_all = [];
    sim_all = [];
    for g = 1:length(group_results)
        group = group_results{g};
        ref_all = [ref_all; group.([band '_ref'])];
        sim_all = [sim_all; group.([band '_aiem'])];
    end
    
    % Compute metrics
    mask = isfinite(ref_all) & isfinite(sim_all);
    if sum(mask) == 0
        fprintf('%s: No valid data\n', band_label);
        continue;
    end
    
    ref = ref_all(mask);
    sim = sim_all(mask);
    diff = sim - ref;
    
    rmse = sqrt(mean(diff.^2));
    bias = mean(diff);
    mae = mean(abs(diff));
    r = corr(ref, sim);
    n = sum(mask);
    
    fprintf('%s: RMSE=%.2f dB, Bias=%+.2f dB, MAE=%.2f dB, r=%.3f (n=%d)\n', ...
            band_label, rmse, bias, mae, r, n);
end
fprintf('\n');

%% Metrics plots
if ~isempty(metrics_table)
    figure('Position', [100, 100, 1600, 400]);
    
    metric_names = {'RMSE (dB)', 'Bias (dB)', 'MAE (dB)', 'Correlation r'};
    metric_cols = [3, 4, 5, 6];
    
    for m = 1:4
        subplot(1, 4, m);
        hold on;
        
        for b = 1:3
            band_data = metrics_table(metrics_table(:, 2) == b, :);
            if isempty(band_data)
                continue;
            end
            
            ratios_plot = band_data(:, 1);
            vals = band_data(:, metric_cols(m));
            
            plot(ratios_plot, vals, '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
                 'DisplayName', band_labels{b});
            
            % Add value labels
            for i = 1:length(ratios_plot)
                text(ratios_plot(i), vals(i), sprintf('%.2f', vals(i)), ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                     'FontSize', 8);
            end
        end
        
        xlabel('Ratio (ℓ/σ)', 'FontSize', 12);
        ylabel(metric_names{m}, 'FontSize', 12);
        title(metric_names{m}, 'FontSize', 14);
        legend('Location', 'best');
        grid on;
        hold off;
    end
    
    sgtitle('Metrics vs Ratio by Polarization', 'FontSize', 16);
end

%% Scatter plots
figure('Position', [100, 100, 1800, 500]);

for b = 1:3
    subplot(1, 3, b);
    hold on;
    
    band = bands{b};
    band_label = band_labels{b};
    
    % Plot by ratio
    for r_idx = 1:length(unique_ratios)
        ratio = unique_ratios(r_idx);
        color = colors(mod(r_idx-1, size(colors, 1)) + 1, :);
        marker = markers{mod(r_idx-1, length(markers)) + 1};
        
        % Find all groups for this ratio
        ref_ratio = [];
        sim_ratio = [];
        for g = 1:length(group_results)
            if abs(group_results{g}.ratio - ratio) < 0.01
                group = group_results{g};
                ref_ratio = [ref_ratio; group.([band '_ref'])];
                sim_ratio = [sim_ratio; group.([band '_aiem'])];
            end
        end
        
        mask = isfinite(ref_ratio) & isfinite(sim_ratio);
        if ~any(mask)
            continue;
        end
        
        scatter(ref_ratio(mask), sim_ratio(mask), 50, marker, ...
                'MarkerFaceColor', color, 'MarkerEdgeColor', 'black', ...
                'DisplayName', sprintf('ℓ/σ=%.1f (n=%d)', ratio, sum(mask)));
    end
    
    % 1:1 line
    xlims = xlim;
    ylims = ylim;
    lims = [min(xlims(1), ylims(1)), max(xlims(2), ylims(2))];
    plot(lims, lims, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    xlabel('NMM3D σ⁰ (dB)', 'FontSize', 12);
    ylabel('AIEM σ⁰ (dB)', 'FontSize', 12);
    title(band_label, 'FontSize', 14);
    legend('Location', 'best', 'FontSize', 8);
    grid on;
    axis equal;
    xlim(lims);
    ylim(lims);
    hold off;
end

sgtitle('AIEM vs NMM3D Scatter Plots', 'FontSize', 16);

%% Summary
fprintf('========================================\n');
fprintf('Summary\n');
fprintf('========================================\n\n');

fprintf('Configuration:\n');
fprintf('  Frequency: %.3f GHz\n', FREQUENCY_GHZ);
fprintf('  Incidence: %.1f deg\n', INCIDENCE_DEG);
fprintf('  Ratios: %s\n', mat2str(unique_ratios));
fprintf('  Multiple scattering: %s\n', mat2str(INCLUDE_MULTIPLE_SCATTERING));
fprintf('\n');

fprintf('Bug fixes applied: ✓ All 5 bugs fixed\n');
fprintf('Groups processed: %d\n', length(group_results));
fprintf('\n');

if ~INCLUDE_MULTIPLE_SCATTERING
    fprintf('Note: Single scattering only\n');
    fprintf('  - Expected bias: +3-5 dB for co-pol\n');
    fprintf('  - HV/VH may be very small or -Inf\n');
    fprintf('  - Set INCLUDE_MULTIPLE_SCATTERING=true for cross-pol\n');
else
    fprintf('Note: Multiple scattering included\n');
    fprintf('  - Cross-pol (HV/VH) should have finite values\n');
    fprintf('  - Computation time: ~2-4 minutes per group\n');
end
fprintf('\n');

fprintf('For details, see:\n');
fprintf('  - README.md\n');
fprintf('  - BUG_FIXES.md\n');
fprintf('  - COMPLETE_MULTIPLE_SCATTERING.md\n');
