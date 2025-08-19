function varargout = raster2D_ani(s0, t0, t1, Ne1)
%RASTER2D_ANI Animate a 2‑D raster from spike times.
%   RASTER2D_ANI(s0, t0, t1, Ne1) animates spikes binned in 2 ms steps
%   between times t0 and t1 (inclusive). If s0 has two rows, it is
%   interpreted as [t; linearIndex]; if it has three rows, it is
%   interpreted as [t; x; y]. Ne1 is the number of neurons per dimension.
%
%   - For s0 with 2 rows, indices are mapped to (x,y) on an Ne1-by-Ne1 grid:
%       x = ceil(idx / Ne1),  y = mod(idx-1, Ne1) + 1.
%   - The very first frame initializes axes to show only the first half of
%     the linear indices (idx < Ne1^2/2), with x-limits [0, Ne1/2] and
%     y-limits [0, Ne1]. Subsequent frames plot all valid spikes.
%
%   A = RASTER2D_ANI(...) also returns the movie frame struct array A,
%   where each element is compatible with VIDEOFRAME/GETFRAME.
%
%   Inputs
%     s0  : 2×Ns or 3×Ns spike array (milliseconds recommended)
%     t0  : start time (same units as s0(1,:))
%     t1  : end time   (same units as s0(1,:), t1 ≥ t0)
%     Ne1 : # of neurons per dimension (positive integer)
%
%   Output (optional)
%     A   : 1×Nframes struct('cdata',..., 'colormap', ...)
%
%   Notes
%     • Bin size is fixed at 2 (same units as t0/t1), matching original code.
%     • Axis ticks are kept as in your script; adjust inside the code if needed.

    % -----------------------
    % Parameters (kept as in original)
    % -----------------------
    dta = 2;                       % Bin size for raster animation
    timea = t0:dta:t1;             % Frame times (bin right edges)
    numframes = numel(timea);

    % -----------------------
    % Basic validation
    % -----------------------
    if ~(isnumeric(s0) && ismatrix(s0))
        error('s0 must be a numeric 2×N or 3×N matrix.');
    end
    if size(s0,1)~=2 && size(s0,1)~=3
        error('s0 must have 2 or 3 rows: [t; idx] or [t; x; y].');
    end
    if ~isscalar(t0) || ~isscalar(t1) || ~isfinite(t0) || ~isfinite(t1) || t1 < t0
        error('t0 and t1 must be finite scalars with t1 >= t0.');
    end
    if ~isscalar(Ne1) || ~isfinite(Ne1) || Ne1 <= 0 || Ne1 ~= floor(Ne1)
        error('Ne1 must be a positive integer.');
    end

    % -----------------------
    % Figure setup
    % -----------------------
    fig1 = figure('Color','w','Position',[350 100 550 450]);
    ax = axes('Parent',fig1);
    hold(ax, 'on');

    % Preallocate movie frames
    A(1:numframes) = struct('cdata', [], 'colormap', []);

    % -----------------------
    % First frame: initialize plot & axes (preserve original behavior)
    % -----------------------
    i = 1;
    Is = find( s0(1,:) <= timea(i) & s0(1,:) > timea(i)-dta & ...
               s0(2,:) > 0 & s0(2,:) < Ne1^2/2 );  %#ok<EFIND>

    if size(s0,1) == 2
        if ~isempty(Is)
            x = ceil(s0(2,Is)/Ne1);
            y = mod(s0(2,Is)-1, Ne1) + 1;
        else
            x = 0; y = 0; % empty placeholder (kept from original)
        end
    else
        x = s0(2,Is);
        y = s0(3,Is);
        if isempty(Is), x = 0; y = 0; end
    end

    h = plot(ax, x, y, 'k.', 'MarkerSize', 4);
    % Fix axes as in original first frame
    axis(ax, [0, Ne1/2, 0, Ne1]);
    xlabel(ax, 'cell index 1', 'FontSize', 18);
    ylabel(ax, 'cell index 2', 'FontSize', 18);
    set(ax, 'FontSize', 15, 'XTick', [0 100 200], 'YTick', [0 100 200]);

    title(ax, sprintf('t=%d msec', round(timea(i)) - t0));
    A(i) = getframe(fig1);

    % -----------------------
    % Remaining frames
    % -----------------------
    for i = 2:numframes
        % Bin spikes in (timea(i)-dta, timea(i)] and s0(2,:)>0
        Is = find( s0(1,:) <= timea(i) & s0(1,:) > timea(i)-dta & s0(2,:) > 0 ); %#ok<EFIND>

        if size(s0,1) == 2
            if ~isempty(Is)
                x = ceil(s0(2,Is)/Ne1);
                y = mod(s0(2,Is)-1, Ne1) + 1;
            else
                x = []; y = [];
            end
        else
            x = s0(2,Is);
            y = s0(3,Is);
        end

        % Update data (faster than refreshdata)
        set(h, 'XData', x, 'YData', y);
        drawnow;

        title(ax, sprintf('t=%d msec', round(timea(i)) - t0));
        A(i) = getframe(fig1);
    end

    % -----------------------
    % Optional output
    % -----------------------
    if nargout == 1
        varargout{1} = A;
    end
end
