function saveFigure(varargin)
%HF_VIEWSAVE  Save the current figure to disk and optionally open it in a viewer.
%
% Syntax (name–value options; recommended)
%   HF_viewsave('-name','Fig01','-format','pdf','-view',1, ...
%               '-path','./result','-res',400, ...
%               '-label','Figure 1','-colorspace','', '-save',1)
%
% Legacy syntax (kept for backward compatibility)
%   HF_viewsave(projectName, baseName, viewFlag, format, tag, res)
%     projectName : string used to build default path
%     baseName    : file name (without extension)
%     viewFlag    : 0/1 whether to open after saving
%     format      : e.g., 'pdf','epsc2','epsc','jpg','tif',...
%     tag         : explicit extension (if omitted, inferred from format)
%     res         : DPI (numeric)
%
% Description
%   Saves the current figure (gcf) to disk in the requested format, with
%   sensible defaults and cross-platform viewer opening. Supports adding a
%   bold label beneath the figure.
%
% Inputs (name–value options; all optional unless noted)
%   -name        : (required) base file name, no extension (char/string)
%   -path        : output directory (char/string). Empty→current folder.
%   -view        : logical {0/1}. If 1, attempt to open after saving.
%   -format      : graphics device for PRINT (default 'pdf').
%   -tag         : explicit file extension override (default = inferred).
%   -res         : DPI resolution (default 400; set automatically if needed).
%   -label       : char/string. If nonempty, add bold label below the figure.
%   -colorspace  : char/string passed to PRINT (e.g., '-cmyk' or '-rgb'). '' = none.
%   -save        : logical {0/1}. If 0, skip printing; still view if -view=1.
%
% Outputs
%   None. Side effects: writes files to disk; may open a viewer.
%
% Example
%   HF_viewsave('-name','SuppFig_Summary','-path','result','-format','pdf','-view',1);
%
% Notes
%   - EPS→PDF conversion on UNIX uses 'pstopdf' if available, otherwise 'ps2pdf'.
%   - macOS opens with Preview; Windows uses the default file association;
%     Linux uses xdg-open (falls back to evince/eog).
%
% See also: PRINT, OPEN

% =========================
% Constants (no magic numbers)
% =========================
DEFAULT_FORMAT       = 'pdf';
DEFAULT_RES_DPI      = 400;
DEFAULT_VIEW         = 0;
DEFAULT_COLORSPACE   = '';
DEFAULT_SAVE_FLAG    = 1;
MIN_BOTTOM_MARGIN_IN = 0.5;   % inches: ensure room for label
LABEL_FONT_SIZE_PT   = 14;
LABEL_FONT_WEIGHT    = 'bold';

% =========================
% Parse & normalize options
% =========================
P = struct();
if useNewStyle(varargin)            % name–value style
    P = parseNameValue(varargin{:});
else                                % legacy positional style
    P = parseLegacy(varargin{:});
end


% Fill defaults & validate
P = applyDefaults(P, DEFAULT_FORMAT, DEFAULT_RES_DPI, DEFAULT_VIEW, ...
    DEFAULT_COLORSPACE, DEFAULT_SAVE_FLAG);
validateOptions(P);

% Normalize output folder & build final path
if isempty(P.path)
    outdir = '';
else
    outdir = ensureTrailingSep(P.path);
    ensureFolderExists(outdir);
end
outfileBase = [outdir, P.name];  % full path without extension

% =========================
% Renderer-dependent tweaks
% =========================
rend = get(gcf,'Renderer');
if ~strcmpi(rend,'painters')
    if isempty(P.res) || ~isfinite(P.res)
        P.res = DEFAULT_RES_DPI;
    end
end

% =========================
% Optional figure label
% =========================
if ~isempty(P.label)
    % Ensure some bottom margin to place the label nicely
    pos = get(gcf,'PaperPosition');  % [left bottom width height], inches
    if pos(2) < MIN_BOTTOM_MARGIN_IN
        pos(2) = MIN_BOTTOM_MARGIN_IN;
        set(gcf,'PaperPosition',pos);
    end

    ax = axes('OuterPosition',[0 0 1 1], 'Visible','off'); 
    text(0.5, -0.5, P.label, ...
        'Units','normalized', ...
        'FontSize', LABEL_FONT_SIZE_PT, ...
        'FontWeight', LABEL_FONT_WEIGHT, ...
        'HorizontalAlignment','center', ...
        'Parent', ax); 
end

% =========================
% Determine extension/tag
% =========================
fmtLower   = lower(string(P.format));
tag        = inferTag(fmtLower, P.tag);
outfileTag = [outfileBase, '.', tag];

% =========================
% Save/print
% =========================
if P.save
    if exist(outfileTag, 'file'), delete(outfileTag); end

    args = {gcf, ['-d', char(fmtLower)]};
    if ~isempty(P.colorspace)
        args{end+1} = P.colorspace; 
    end
    if ~isempty(P.res) && isfinite(P.res)
        args{end+1} = ['-r', num2str(P.res)];
    end
    args{end+1} = outfileTag;

    try
        % print(args{:});
        print(gcf,['-d',P.format],[P.colorspace],['-r',num2str(P.res)],[outfileBase,'.',tag]);
    catch ME
        error('Failed to print figure to "%s": %s', outfileTag, ME.message);
    end

    % If EPS, convert to PDF to match original intent for viewing convenience
    if strcmpi(tag, 'eps')
        convertEpsToPdf(outfileBase);
    end
end

% =========================
% View (open) if requested
% =========================
if P.view
    targetPath = pickBestToOpen(outfileBase, tag); % prefer PDF if present
    openInViewer(targetPath);
end

end

% ======================================================================
% Helpers
% ======================================================================

function tf = useNewStyle(vs)
% Detect name–value style by checking if the 3rd argument looks like an option flag.
    tf = ~isempty(vs) ...
         && numel(vs) >= 3 ...
         && (ischar(vs{3}) || (isstring(vs{3}) && isscalar(vs{3})));
end

function P = parseNameValue(varargin)
% Parse name–value pairs; accepts leading '-' in names but does not require it.
    if mod(nargin,2) ~= 0
        error('HF_viewsave:ArgPairs', ...
            'Name–value arguments must come in pairs.');
    end
    P = struct();
    for k = 1:2:numel(varargin)
        key = varargin{k};
        val = varargin{k+1};
        if ~(ischar(key) || (isstring(key) && isscalar(key)))
            error('HF_viewsave:BadOptionName', ...
                'Option name at position %d must be a char or string.', k);
        end
        key = char(key);
        if ~isempty(key) && key(1) == '-'
            key = key(1+1:end); % drop leading dash
        end
        P.(lower(key)) = val;
    end
end

function P = parseLegacy(varargin)
% Legacy positional interface: (projectName, baseName, viewFlag, format, tag, res)
    if nargin < 3
        error('HF_viewsave:LegacyArgs', ...
            'Legacy form requires at least 3 inputs: projectName, baseName, viewFlag.');
    end
    projectName = varargin{1};
    baseName    = varargin{2};
    viewFlag    = varargin{3};

    if ~(ischar(baseName) || (isstring(baseName) && isscalar(baseName)))
        error('HF_viewsave:BaseNameType', 'baseName must be a char or string.');
    end

    % Try to reproduce the original path build. If setgetDirs exists, use it.
    outPath = '';
    try
        if exist('setgetDirs', 'file') == 2
            try
                D = setgetDirs();
                outPath = fullfile(D.Root, 'project', char(projectName), 'figures');
            catch
                outPath = fullfile(pwd, char(projectName), 'figures');
            end
        else
            outPath = fullfile(pwd, char(projectName), 'figures');
        end
    catch
        outPath = fullfile(pwd, char(projectName), 'figures');
    end

    P = struct();
    P.name   = char(baseName);
    P.view   = logical(viewFlag);
    P.path   = outPath;

    if nargin >= 4 && ~isempty(varargin{4}), P.format = varargin{4}; end
    if nargin >= 5 && ~isempty(varargin{5}), P.tag    = varargin{5}; end
    if nargin >= 6 && ~isempty(varargin{6}), P.res    = varargin{6}; end
end

function P = applyDefaults(P, defFmt, defDpi, defView, defColor, defSave)
% Apply defaults when fields are missing.
    if ~isfield(P,'name'),       P.name = '';      end
    if ~isfield(P,'view'),       P.view = defView; end
    if ~isfield(P,'format'),     P.format = defFmt; end
    if ~isfield(P,'tag'),        P.tag = [];       end
    if ~isfield(P,'res'),        P.res = defDpi;   end
    if ~isfield(P,'path'),       P.path = [];      end
    if ~isfield(P,'label'),      P.label = [];     end
    if ~isfield(P,'colorspace'), P.colorspace = defColor;  end
    if ~isfield(P,'save'),       P.save = defSave; end
end

function validateOptions(P)
% Validate option types and values with explicit error messages.
    if ~(ischar(P.name) || (isstring(P.name) && isscalar(P.name))) || isempty(P.name)
        error('HF_viewsave:MissingName', ...
            'You must provide a nonempty -name (base file name without extension).');
    end
    if ~(isscalar(P.view) && (islogical(P.view) || isnumeric(P.view)))
        error('HF_viewsave:ViewType', '-view must be 0 or 1.');
    end
    if ~isempty(P.res)
        if ~(isnumeric(P.res) && isfinite(P.res) && P.res > 0)
            error('HF_viewsave:ResType', '-res must be a positive finite numeric value.');
        end
    end
    if ~(isempty(P.path) || ischar(P.path) || (isstring(P.path) && isscalar(P.path)))
        error('HF_viewsave:PathType', '-path must be empty or a char/string.');
    end
    if ~(ischar(P.format) || (isstring(P.format) && isscalar(P.format)))
        error('HF_viewsave:FormatType', '-format must be a char/string.');
    end
    if ~(isempty(P.label) || ischar(P.label) || (isstring(P.label) && isscalar(P.label)))
        error('HF_viewsave:LabelType', '-label must be empty or a char/string.');
    end
    if ~(isscalar(P.save) && (islogical(P.save) || isnumeric(P.save)))
        error('HF_viewsave:SaveType', '-save must be 0 or 1.');
    end
end

function out = ensureTrailingSep(p)
% Ensure trailing file separator.
    p = char(p);
    if isempty(p), out = p; return; end
    if p(end) ~= filesep
        out = [p, filesep];
    else
        out = p;
    end
end

function ensureFolderExists(d)
% Create directory if it does not exist.
    if ~isempty(d) && exist(d,'dir') ~= 7
        ok = mkdir(d);
        if ~ok
            error('HF_viewsave:MkdirFailed', ...
                'Failed to create output directory: %s', d);
        end
    end
end

function tag = inferTag(fmtLower, explicitTag)
% Map MATLAB device names to common extensions.
    if ~isempty(explicitTag)
        tag = char(lower(string(explicitTag)));
        return
    end
    fmtLower = char(fmtLower);
    switch fmtLower
        case {'ps','psc','psc2'}
            tag = 'ps';
        case {'eps','epsc2','epsc'}
            tag = 'eps';
        case {'meta'}
            tag = 'emf';
        otherwise
            tag = fmtLower; % e.g., 'pdf','jpg','tif','png',...
    end
end

function convertEpsToPdf(outfileBase)
% Convert OUTFILEBASE.eps → OUTFILEBASE.pdf using available tools.
    epsPath = [outfileBase, '.eps'];
    pdfPath = [outfileBase, '.pdf'];
    if ~exist(epsPath,'file')
        return
    end
    if exist(pdfPath,'file')
        delete(pdfPath);
    end

    if ispc
        error('HF_viewsave:DistillerMissing', ...
            ['EPS→PDF conversion on Windows requires Acrobat Distiller.\n' ...
             'Please convert manually or install Ghostscript and run ps2pdf.']);
    else
        [st1,~] = system(sprintf('pstopdf "%s" -o "%s"', epsPath, pdfPath));
        if st1 ~= 0
            [st2,~] = system(sprintf('ps2pdf "%s" "%s"', epsPath, pdfPath));
            if st2 ~= 0
                warning('HF_viewsave:NoConverter', ...
                    'Could not convert EPS to PDF (no pstopdf/ps2pdf). EPS left as-is.');
            end
        end
    end
end

function target = pickBestToOpen(base, tag)
% Prefer PDF if it exists; else open the tagged file.
    pdfPath = [base, '.pdf'];
    tagged   = [base, '.', tag];
    if exist(pdfPath,'file')
        target = pdfPath;
    else
        target = tagged;
    end
end

function openInViewer(fullpath)
% Cross-platform open in default viewer.
    if ~exist(fullpath,'file')
        warning('HF_viewsave:OpenMissing', 'Cannot open viewer: file not found: %s', fullpath);
        return
    end

    if ismac
        system(sprintf('open -a Preview "%s"', fullpath));
    elseif ispc
        try
            winopen(fullpath);
        catch
            system(sprintf('start "" "%s"', fullpath));
        end
    else
        st = system(sprintf('xdg-open "%s" >/dev/null 2>&1 &', fullpath));
        if st ~= 0
            [~,~,ext] = fileparts(fullpath);
            if strcmpi(ext,'.pdf')
                system(sprintf('evince "%s" >/dev/null 2>&1 &', fullpath));
            else
                system(sprintf('eog "%s" >/dev/null 2>&1 &', fullpath));
            end
        end
    end
end
