function sim_default()
% Minimal, no-input driver: generate V1 inputs and run 3 simulations.

% --- Paths
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
dataDir  = fullfile(repoRoot, 'data');
utilsDir = fullfile(repoRoot, 'src', 'utils');
addpath(utilsDir);

% --- Config
simulation_type = 'default';
weight_name     = 'default';
W_fname = fullfile(dataDir, 'weightSpatRecTrans_sigmaRX_0d10_sigmaRR_0d20_Pts_0d15_0d15_tuning_th_0d60.mat');

% --- Generate V1 firing rates
theta1 = 180;  % degrees
theta2 =  90;  % degrees
c1 = 0.5;
c2 = 0.5;

out = utils_simulation.gen_V1_FiringRates(theta1/180*pi, theta2/180*pi, c1, c2); % expects fields fr_left, fr_right, fr_both

% --- Save three V1 inputs (1=left-only, 2=right-only, 3=both)
V1_files = cell(1,3);
for i = 1:3
    base = strrep(sprintf('V1fr_ori_%d_%d_con_%.1f_%.1f_%d', theta1, theta2, c1, c2, i),'.','d');
    V1_files{i} = fullfile(dataDir, [base '.mat']);
end

% Save as variable 'fr' to match downstream expectations
fr = out.fr_left;  save(V1_files{1}, 'fr');  % case 1
fr = out.fr_right; save(V1_files{2}, 'fr');  % case 2
fr = out.fr_both;  save(V1_files{3}, 'fr');  % case 3

% --- Durations for each case
Ts = [2e3, 2e3, 2e3];  % ms

% --- Run simulations
for i = 1%:3
    fr_fname = V1_files{i};
    T        = Ts(i);
    [~, v1_tag] = fileparts(fr_fname);

    % Build output filename (replace dots with 'd' so it's filesystem-friendly)
    simname  = sprintf('Sim_%s_Weight_%s_V1_%s_T_%.1f', simulation_type, weight_name, v1_tag, T/1e3);
    simname  = strrep(simname, '.', 'd');
    outfile  = fullfile(dataDir, [simname '.mat']);

    param_changes = { ...
        'filename', outfile; ...
        'T',        T; ...
        'fr_fname', fr_fname; ...
        'W_fname',  W_fname ...
    };

    spiking_simulation.spiking_sim(simulation_type, [], param_changes);
end
end
