function compile_all
% Compile all EIF normalization MEX files into sim_script/mex
% Folder layout (relative to this script):
%   sim_script/
%     compile_all.m
%     mex/               (created if missing)
%   src/
%     common/            (EIF_common.c, EIF_common.h)
%     variants/          (multiple EIF_normalization_*.c)

clc;

thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));          % up from sim_script → repo root
simScriptDir = fullfile(repoRoot, 'sim_script');
mexOutDir    = fullfile(simScriptDir, 'mex');
srcDir       = fullfile(repoRoot, 'src');
codesDir     = fullfile(srcDir, 'simulation/C_source_codes');

if ~exist(mexOutDir, 'dir'); mkdir(mexOutDir); end

% Detect platform suffix (just for pretty printing)
arch = mexext; % e.g., 'mexmaca64', 'mexa64', 'mexw64'

fprintf('Building with ''%s''.\n\n', mex.getCompilerConfigurations('C','Selected').ShortName);

%% 1) Compile the standalone default (NO common linkage)
defaultSrc = fullfile(codesDir, 'EIF_normalization_Default.c');
assert(exist(defaultSrc, 'file')==2, 'Missing %s', defaultSrc);

fprintf('Compiling standalone default:\n');
try
    mex('-silent', '-outdir', mexOutDir, defaultSrc);
    tgt = fullfile(mexOutDir, ['EIF_normalization_Default.', arch]);
    fprintf('  ✅ EIF_normalization_Default -> %s\n\n', tgt);
catch ME
    fprintf(2, '  ❌ Failed: EIF_normalization_Default.c\n    %s\n\n', ME.message);
end

%% 2) Compile all other variants (link WITH common)
S = dir(fullfile(codesDir, 'EIF_normalization_*.c'));
% exclude Default from the variant loop
S = S(~strcmpi({S.name}, 'EIF_normalization_Default.c'));

commonC = fullfile(codesDir, 'EIF_common.c');
commonH = fullfile(codesDir, 'EIF_common.h');
assert(exist(commonC,'file')==2 && exist(commonH,'file')==2, ...
    'Expected common sources not found in %s', codesDir);

if ~isempty(S)
    fprintf('Compiling variants (linked with EIF_common.c):\n');
end

ok = 0; fail = 0;
for k = 1:numel(S)
    src = fullfile(S(k).folder, S(k).name);
    [~, base] = fileparts(src);

    try
        % Note: add include path to find EIF_common.h
        mex('-silent', '-outdir', mexOutDir, ['-I', codesDir], src, commonC);
        tgt = fullfile(mexOutDir, [base, '.', arch]);
        fprintf('  ✅ %s -> %s\n', base, tgt);
        ok = ok + 1;
    catch ME
        fprintf(2, '  ❌ Failed: %s\n    %s\n', S(k).name, ME.message);
        fail = fail + 1;
    end
end

%% Summary + directory check
fprintf('\nBuild summary\n  Passed: %d\n  Failed: %d\n', ok+1, fail); % +1 for Default

% Show what actually landed in the output folder so it matches Finder
listing = dir(fullfile(mexOutDir, ['EIF_normalization_*.', arch]));
if isempty(listing)
    fprintf(2, '  (No MEX files found in %s)\n', mexOutDir);
else
    fprintf('  Files in %s:\n', mexOutDir);
    for i = 1:numel(listing)
        fprintf('   - %s\n', listing(i).name);
    end
end

end
