function results = get_Fisher_information(data_folder, file_prefix, file_ind1, Nfile, Ntr, ds, batch_ID, opts)
% GET_FISHER_INFORMATION
% Syntax:
%   results = get_Fisher_information(data_folder, file_prefix, file_ind1, Nfile, Ntr, ds, batch_ID)
%   results = get_Fisher_information(..., opts)
% Description:
%   Computes naive and bias-corrected Fisher Information (FI) between two
%   stimulus conditions using trial-wise binary responses X1 (cond1) and
%   X2 (cond2). For a set of neuron subset sizes, it samples neuron IDs
%   from 'ind1' and aggregates X1/X2 across multiple files, then computes
%   df = (mean(X2) - mean(X1)) / ds, C = (cov(X1) + cov(X2))/2, and
%   FI_naive = df * pinv(C) * df'. A small-sample bias-corrected estimate
%   FI_BC is also returned.
% Inputs:
%   data_folder : string, folder containing input/output MAT files (with trailing filesep or not).
%   file_prefix : string, common prefix for per-file data, e.g. 'Default_..._'.
%   file_ind1   : string, MAT file containing vector 'ind1' of neuron indices to sample from.
%   Nfile       : positive integer, number of cond1/cond2 files to load (IDs 1..Nfile).
%   Ntr         : positive integer, total #trials across all files (if [], inferred).
%   ds          : positive scalar, stimulus step size used for FI.
%   batch_ID    : scalar or string used in file names after an underscore.
%   opts        : struct of optional parameters:
%                 - Nneuron: vector of subset sizes (default [1 2 4 8 16 31 62 125 250 500 1000 2000 4000])
%                 - save_prefix: default 'FI_'
%                 - use_parfor: logical (default false)
%                 - verbose: logical (default true)
% Outputs (struct array, one per subset size):
%   results(rr) with fields: Nid, ds, Nfile, Ntr, FInaive, FI_BC, neuron_ids, out_file
% Example:
%   results = get_Fisher_information(d, p, 'ind1.mat', 4, [], 0.1, 7);
%   % loads cond1_1..4 & cond2_1..4 and saves FI_* files per subset size
%
% Notes (robustness & performance):
%   - Validates presence of variables 'X1','X2' in files, and 'ind1' in file_ind1.
%   - Infers trials-per-file from the first 'cond1_1_*' file if Ntr is [].
%   - Converts to double only when needed for covariance to save memory.
%   - You can set opts.use_parfor=true to parallelize over subset sizes.

  arguments
    data_folder (1,1) string
    file_prefix (1,1) string
    file_ind1   (1,1) string
    Nfile (1,1) double {mustBeInteger, mustBePositive}
    Ntr   double {mustBeNonnegative} = []
    ds    (1,1) double {mustBePositive}
    batch_ID
    opts.Nneuron double = [1,2,4,8,16,31,62,125,250,500,1000,2000,4000]
    opts.save_prefix (1,1) string = "FI_"
    opts.use_parfor (1,1) logical = false
    opts.verbose (1,1) logical = true
  end

  if ~endsWith(data_folder, filesep)
      data_folder = data_folder + filesep;
  end

  % Load neuron pool
  S = load(data_folder + file_ind1, 'ind1');
  if ~isfield(S, 'ind1')
      error('%s: %s must contain variable ''ind1''.', mfilename, file_ind1);
  end
  ind1 = S.ind1(:);  % column vector
  if isempty(ind1)
      error('%s: ind1 in %s is empty.', mfilename, file_ind1);
  end

  % Infer trials-per-file from the first cond1 file if needed
  firstFile = sprintf('%scond1_%d_%s.mat', file_prefix, 1, string(batch_ID));
  F = load(data_folder + firstFile, 'X1');
  if ~isfield(F, 'X1')
      error('%s: %s missing variable X1.', mfilename, firstFile);
  end
  trials_per_file = size(F.X1, 1);
  if isempty(Ntr)
      Ntr = trials_per_file * Nfile;
  else
      if Ntr ~= trials_per_file * Nfile
          error('%s: Ntr (%d) != trials_per_file*Nfile (%d).', mfilename, Ntr, trials_per_file*Nfile);
      end
  end

  Nsizes = numel(opts.Nneuron);
  results = repmat(struct('Nid',[],'ds',ds,'Nfile',Nfile,'Ntr',Ntr,'FInaive',[],'FI_BC',[],'neuron_ids',[],'out_file',''), Nsizes, 1);

  if opts.use_parfor
      parfor rr = 1:Nsizes
          results(rr) = local_compute_one_size(rr, opts.Nneuron(rr), ind1, data_folder, file_prefix, batch_ID, Nfile, trials_per_file, Ntr, ds, opts.save_prefix, opts.verbose);
      end
  else
      for rr = 1:Nsizes
          results(rr) = local_compute_one_size(rr, opts.Nneuron(rr), ind1, data_folder, file_prefix, batch_ID, Nfile, trials_per_file, Ntr, ds, opts.save_prefix, opts.verbose);
      end
  end
end

function out = local_compute_one_size(rr, Nneur, ind1, data_folder, file_prefix, batch_ID, Nfile, trials_per_file, Ntr, ds, save_prefix, verbose)
  if Nneur > numel(ind1)
      error('Requested Nneuron=%d exceeds available %d in ind1.', Nneur, numel(ind1));
  end

  % Sample neuron IDs (fallback to randperm if randsample is unavailable)
  if exist('randsample','file') == 2 %#ok<EXIST>
      Nid = randsample(ind1, Nneur);
  else
      Nid = ind1(randperm(numel(ind1), Nneur));
  end

  if verbose, fprintf('[%d] Loading %d files for N=%d...
', rr, Nfile, Nneur); tic; end

  % Preallocate trial-by-neuron matrices (int8 to match source, reduces memory)
  X1 = zeros(Ntr, Nneur, 'int8');
  X2 = zeros(Ntr, Nneur, 'int8');

  rowStart = 1;
  for ID = 1:Nfile
      f1 = sprintf('%scond1_%d_%s.mat', file_prefix, ID, string(batch_ID));
      f2 = sprintf('%scond2_%d_%s.mat', file_prefix, ID, string(batch_ID));

      D1 = load(data_folder + f1, 'X1');
      D2 = load(data_folder + f2, 'X2');
      if ~isfield(D1,'X1'), error('%s missing X1.', f1); end
      if ~isfield(D2,'X2'), error('%s missing X2.', f2); end

      if size(D1.X1,2) < max(Nid) || size(D2.X2,2) < max(Nid)
          error('File %s/%s has fewer columns than max(ind1).', f1, f2);
      end

      rows = size(D1.X1,1);
      if rows ~= trials_per_file || rows ~= size(D2.X2,1)
          error('Trial count mismatch between %s and %s or expected %d.', f1, f2, trials_per_file);
      end

      X1(rowStart:rowStart+rows-1, :) = D1.X1(:, Nid);
      X2(rowStart:rowStart+rows-1, :) = D2.X2(:, Nid);
      rowStart = rowStart + rows;
  end

  if verbose, fprintf('[%d] Finished loading. Computing FI...
', rr); end

  % Means and covariance (convert to double for covariance)
  fm1 = mean(X1, 1);
  fm2 = mean(X2, 1);
  COV = (cov(double(X1)) + cov(double(X2))) / 2;
  df  = (fm2 - fm1) / ds;           % 1 x N

  % Naive FI
  FINaive = df * pinv(COV) * df';   % scalar

  % Bias-corrected FI (see Kanitscheider et al., 2015;/original code)
  N = Nneur;
  FI_BC = FINaive * (2*Ntr - N - 3) / (2*Ntr - 2) - 2*N/(Ntr * ds^2);

  out = struct('Nid',Nid,'ds',ds,'Nfile',Nfile,'Ntr',Ntr,'FInaive',FINaive,'FI_BC',FI_BC,
               'neuron_ids',Nid,'out_file','');

  % Save
  outname = sprintf('%s%s%d_%s.mat', save_prefix, file_prefix, Nneur, string(batch_ID));
  save(data_folder + outname, 'Nid','ds','Nfile','FInaive','FI_BC','Ntr');
  out.out_file = data_folder + outname;

  if verbose, t=toc; fprintf('[%d] Done N=%d in %.2fs. Saved to %s
', rr, Nneur, t, out.out_file); end
end
