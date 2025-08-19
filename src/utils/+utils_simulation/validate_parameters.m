function validate_parameters(params, options)
% VALIDATE_PARAMETERS - Check parameter consistency and requirements
%
% H1 Line: Validate simulation parameters and options
%
% Syntax:
%   validate_parameters(params, options)
%
% Description:
%   Checks that all required parameters are present and consistent.
%   Throws informative errors if validation fails.
%
% Inputs:
%   params  - Parameter structure
%   options - Options structure
%
% Example:
%   validate_parameters(params, options);

    % Check required fields for saving
    if options.save && ~isfield(params, 'filename')
        error('SpikingSim:MissingFilename', ...
              'Options.save is true but no filename specified');
    end
    
    % Check correlation computation requirements
    if options.CompCorr && ~isfield(params, 'Nc')
        error('SpikingSim:MissingNc', ...
              'Options.CompCorr is true but Nc not specified');
    end
    
    % Check file loading requirements
    if options.loadfr1 && ~isfield(params, 'fr_fname')
        error('SpikingSim:MissingFrFile', ...
              'Options.loadfr1 is true but fr_fname not specified');
    end
    
    if options.useWfile && ~isfield(params, 'W_fname')
        error('SpikingSim:MissingWFile', ...
              'Options.useWfile is true but W_fname not specified');
    end
    
    % Validate network size
    if params.N > 200000
        error('SpikingSim:NetworkTooLarge', ...
              'Network size N=%d exceeds maximum (200000)', params.N);
    end
    
    % Validate numerical parameters
    if params.dt <= 0
        error('SpikingSim:InvalidTimeStep', ...
              'Time step dt must be positive');
    end
    
    if params.T <= 0
        error('SpikingSim:InvalidDuration', ...
              'Simulation duration T must be positive');
    end
    
end
