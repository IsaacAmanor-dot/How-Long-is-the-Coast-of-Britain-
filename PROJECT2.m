%% AMS 595 / DCS 525 — Project 2  % Task 1: Mandelbrot Fractal — it = fractal(c)  % Description of purpose
clear; clc; close all;  % Clear workspace, command window, and close figures

%% Configuration  % Parameters for Task 1 visualization and escape test
MAX_ITER = 100;  % Maximum number of iterations
BAILOUT   = 2.0; % Divergence threshold for |z|

xlim_viz  = [-2.5, 1.0];  % Real-axis limits for visualization
ylim_viz  = [-1.5, 1.5];  % Imag-axis limits for visualization
nx = 800;  % Number of points along real axis (resolution)
ny = 600;  % Number of points along imaginary axis (resolution)

%% Self-Tests  % Quick checks for correctness of fractal()
fprintf('Running Task 1 sanity checks...\n');  % Status message
it0 = fractal(0 + 0i, MAX_ITER, BAILOUT);  % Test with c = 0 (inside set)
assert(it0 == 0, 'Expected it=0 for c=0 (did not diverge).');  % Verify
it_div = fractal(1 + 1i, MAX_ITER, BAILOUT);  % Test with c = 1+1i (diverges)
assert(it_div > 0 && it_div <= MAX_ITER, 'Expected divergence for c=1+1i.');  % Verify
it2 = fractal(2 + 0i, MAX_ITER, BAILOUT);  % Test with c = 2 (diverges)
assert(it2 > 0, 'Expected divergence for c=2.');  % Verify
fprintf('Sanity checks passed. Examples:\n');  % Print message
fprintf('  c = 0         → it = %d (0 = no divergence within %d iters)\n', it0, MAX_ITER);  % Example output
fprintf('  c = 1 + 1i    → it = %d (diverges)\n', it_div);  % Example output
fprintf('  c = 2 + 0i    → it = %d (diverges)\n', it2);  % Example output

%% Optional Visualization  % Render escape-time image for the Mandelbrot set
do_viz = true;  % Toggle visualization
if do_viz  % If visualization is enabled
    fprintf('Generating optional Mandelbrot visualization (%dx%d)...\n', ny, nx);  % Status message
    x = linspace(xlim_viz(1), xlim_viz(2), nx);  % Real-axis grid
    y = linspace(ylim_viz(1), ylim_viz(2), ny);  % Imag-axis grid
    it_img = zeros(ny, nx, 'uint16');  % Preallocate escape-time image
    for iy = 1:ny  % Loop over rows
        yy = y(iy);  % Current imaginary coordinate
        c_row = x + 1i*yy;  % Row of complex numbers
        it_row = arrayfun(@(cc) fractal(cc, MAX_ITER, BAILOUT), c_row);  % Compute escape times for row
        it_img(iy, :) = uint16(it_row);  % Store in image
    end  % End loop over rows
    figure('Color','w');  % Create new figure
    imagesc(x, y, it_img);  % Plot escape-time image
    axis image xy;  % Equal aspect ratio, y increasing upward
    colormap(hot(256)); colorbar;  % Apply colormap and show colorbar
    title(sprintf('Mandelbrot escape-time (MAX\\_ITER=%d, BAILOUT=%.1f)', MAX_ITER, BAILOUT), 'Interpreter','tex');  % Title
    xlabel('Re(c)'); ylabel('Im(c)');  % Axis labels
end  % End optional visualization
drawnow;  % ensure the figure renders before the next section runs
saveas(gcf, 'T1.png');  % or 'T1.pdf' if you prefer vector output


%% Task 2: auto-bracket + bisection  % Finds boundary with per-x bracketing
% NEW
clc;  % keep console clean but do not kill Task 1 figure  % Housekeeping
MAX_ITER = 100;  % Iterations for fractal
BAILOUT  = 2.0;  % Escape radius
xmin = -2.0; xmax = 1.0;  % X range
nx = 1600;  % Number of x-samples (raise to improve coverage)
x_grid = linspace(xmin, xmax, nx);  % X grid
tol = 1e-6;  % Bisection tolerance on y
max_bisect_iters = 60;  % Max iterations for bisection
ymax_scan = 2.0;  % Absolute cap for vertical scans
dy_init = 0.02;  % Starting step for bracket search
grow = 1.5;  % Geometric growth factor for step size
y_top = nan(size(x_grid));  % Preallocate top boundary
y_bot = nan(size(x_grid));  % Preallocate bottom boundary
hit_top = false(size(x_grid));  % Hit mask for top
hit_bot = false(size(x_grid));  % Hit mask for bottom
for i = 1:nx  % Loop over x
    x = x_grid(i);  % Current x
    fn = indicator_fn_at_x(x, MAX_ITER, BAILOUT);  % Indicator at this x
    [s,e] = find_bracket(fn, 0.0, +1, dy_init, grow, ymax_scan);  % Try to bracket above the axis
    if ~isempty(s)  % If bracket found
        y_top(i) = bisection(fn, s, e, tol, max_bisect_iters);  % Bisect to boundary
        hit_top(i) = true;  % Mark success
    end  % End bracket-above branch
    [s,e] = find_bracket(fn, 0.0, -1, dy_init, grow, ymax_scan);  % Try to bracket below the axis
    if ~isempty(s)  % If bracket found
        y_bot(i) = bisection(fn, s, e, tol, max_bisect_iters);  % Bisect to boundary
        hit_bot(i) = true;  % Mark success
    end  % End bracket-below branch
end  % End loop over x
valid_xt = x_grid(hit_top);  % X where top boundary found
valid_yt = y_top(hit_top);  % Y for top boundary
valid_xb = x_grid(hit_bot);  % X where bottom boundary found
valid_yb = y_bot(hit_bot);  % Y for bottom boundary
fprintf('Top boundary: %d / %d   Bottom boundary: %d / %d\n', numel(valid_xt), nx, numel(valid_xb), nx);  % Coverage report
fprintf('Unique x with any boundary: %d\n', numel(unique([valid_xt, valid_xb])));  % Combined coverage
figure('Color','w');  % New figure
plot(valid_xt, valid_yt, '.', 'MarkerSize', 5); hold on;  % Plot top
plot(valid_xb, valid_yb, '.', 'MarkerSize', 5); hold off;  % Plot bottom
xlabel('x'); ylabel('y'); title('Mandelbrot boundary via robust bisection'); grid on; legend({'Top','Bottom'}, 'Location','best');  % Labels and legend

% Clean and save boundary data from Task 2  % Utility snippet
valid_x = unique([valid_xt(:); valid_xb(:)]);  % Merge top/bottom x and make unique
[valid_x, idx] = sort(valid_x);  % Sort x
yy_top = nan(size(valid_x));  % Preallocate y (top)
yy_bot = nan(size(valid_x));  % Preallocate y (bottom)
[tf,loc] = ismember(valid_xt, valid_x); yy_top(loc(tf)) = valid_yt(tf);  % Map top y to merged x
[tf,loc] = ismember(valid_xb, valid_x); yy_bot(loc(tf)) = valid_yb(tf);  % Map bottom y to merged x
save('boundary_raw.mat','valid_x','yy_top','yy_bot');  % Persist for Task 3/4

%% Task 3: Polynomial fit (order 15) of Mandelbrot top boundary  % Fits y(x)=poly15 over trimmed boundary
clc;  % keep console clean but do not kill Task 2 figure 
have_ws = evalin('base','exist(''valid_xt'',''var'') && exist(''valid_yt'',''var'')');  % Check if data in workspace
if have_ws  % Prefer fresh data from Task 2
    x_in = evalin('base','valid_xt(:)');  % Column x (top)
    y_in = evalin('base','valid_yt(:)');  % Column y (top)
else  % Fall back to saved file
    S = load('boundary_raw.mat','valid_x','yy_top');  % Load saved data
    x_in = S.valid_x(:);  % X from file
    y_in = S.yy_top(:);  % Top y from file
end  % End data sourcing branch
mask = ~isnan(x_in) & ~isnan(y_in);  % Remove NaNs
x = x_in(mask);  % Clean x
y = y_in(mask);  % Clean y
[x, ord] = sort(x); y = y(ord);  % Sort by x
eps_flat = 1e-3;  % Threshold to consider y ~ 0 (flat tails)
keep = abs(y) > eps_flat;  % Keep only curved boundary
xk = x(keep);  % Trimmed x
yk = y(keep);  % Trimmed y
q = quantile(xk, [0.01 0.99]);  % Drop extreme 1% on each side
m2 = xk >= q(1) & xk <= q(2);  % Mask central region
xk = xk(m2); yk = yk(m2);  % Apply mask
order = 15;  % Required polynomial order
p = polyfit(xk, yk, order);  % Coefficients of best-fit polynomial (highest power first)
xs = linspace(min(xk), max(xk), 2000);  % Dense x for smooth curve
ys = polyval(p, xs);  % Evaluate polynomial on dense x
fprintf('Fitted order-%d polynomial to %d points (from %d raw top points).\n', order, numel(xk), numel(x));  % Fit report
fprintf('Fit domain: [%.6f, %.6f]\n', xs(1), xs(end));  % Domain of validity
figure('Color','w');  % New figure
plot(xk, yk, '.', 'MarkerSize', 6); hold on;  % Plot data
plot(xs, ys, '-', 'LineWidth', 1.5); hold off;  % Plot fit
xlabel('x'); ylabel('y'); title('Order-15 polynomial fit to Mandelbrot top boundary'); grid on; legend({'Boundary data','Polynomial fit'}, 'Location','best');  % Labels
save('polyfit_top_order15.mat','p','xs','ys','xk','yk');  % Save coefficients, sample curve, and training data

%% Task 4 driver: compute polynomial boundary length  % Uses p,xk from Task 3
clear; clc;  % Housekeeping
S = load('polyfit_top_order15.mat','p','xk');  % Load fit and training x-range
p = S.p;  % Polynomial coefficients (order 15)
s = min(S.xk);  % Left bound of trusted fit
e = max(S.xk);  % Right bound of trusted fit
L = poly_len(p, s, e);  % Curve length of the fitted boundary
fprintf('Polynomial boundary length on [%.6f, %.6f]:  L = %.10f\n', s, e, L);  % Report length with high precision
xs = linspace(s, e, 2000);  % Dense x-grid over fit domain
dp = polyder(p);  % Derivative coefficients
figure('Color','w'); plot(xs, sqrt(1 + (polyval(dp, xs)).^2), 'LineWidth', 1.2); grid on; xlabel('x'); ylabel('sqrt(1 + (y''(x))^2)'); title('Arc-length integrand over fit domain');  % Plot integrand

%%  Local functions (single definitions; no duplicates) 

function it = fractal(c, MAX_ITER, BAILOUT)  % Escape iteration count; 0 if no escape within MAX_ITER
    if nargin < 2 || isempty(MAX_ITER), MAX_ITER = 100; end  % Default iterations
    if nargin < 3 || isempty(BAILOUT),  BAILOUT   = 2.0; end  % Default bailout
    z = complex(0.0, 0.0);  % Initialize z = 0
    it = 0;  % Default output: no divergence
    for k = 1:MAX_ITER  % Iterate up to max iterations
        z = z*z + c;  % Mandelbrot update
        if abs(z) > BAILOUT  % Check divergence
            it = k;  % Record divergence iteration
            return;  % Exit function early
        end  % End bailout check
    end  % End iteration loop
end  % End fractal

function fn = indicator_fn_at_x(x, MAX_ITER, BAILOUT)  % Indicator along vertical line at fixed x
    fn = @(y) (fractal(x + 1i*y, MAX_ITER, BAILOUT) > 0)*2 - 1;  % +1 outside (diverges), -1 inside (bounded within MAX_ITER)
end  % End indicator_fn_at_x

function [s,e] = find_bracket(fn, y0, dirn, dy, grow, ymax)  % Find [s,e] with sign change starting at y0 in direction dirn
    si = fn(y0);  % Sign at start
    if si == 0, s = y0; e = y0; return; end  % Exact boundary (rare)
    y = y0;  % Current y
    step = dy*dirn;  % Initial step with direction
    s = []; e = [];  % Default (empty means not found)
    while abs(y) <= ymax  % Scan until cap
        y_next = y + step;  % Next probe location
        sj = fn(y_next);  % Sign at next
        if sj == 0  % Exact boundary hit
            s = min(y,y_next); e = max(y,y_next); return;  % Return tight bracket
        end  % End exact-boundary check
        if si*sj < 0  % Sign change detected
            s = min(y,y_next); e = max(y,y_next); return;  % Return bracket
        end  % End sign-change check
        y = y_next;  % Advance
        step = step*grow;  % Grow step geometrically
    end  % End scan loop
end  % End find_bracket

function m = bisection(fn_f, s, e, tol, max_iters)  % Bisection with guaranteed sign change on [s,e]
    fs = fn_f(s); fe = fn_f(e);  % Endpoint signs
    if fs*fe >= 0, m = NaN; return; end  % No sign change
    for k = 1:max_iters  % Iterate
        m = 0.5*(s + e);  % Midpoint
        fm = fn_f(m);  % Mid sign
        if abs(e - s) < tol || fm == 0, return; end  % Stop when small or exact
        if fs*fm < 0  % Root in [s,m]
            e = m; fe = fm;  % Tighten right
        else  % Otherwise root in [m,e]
            s = m; fs = fm;  % Tighten left
        end  % End interval selection
    end  % End bisection loop
end  % End bisection

function L = poly_len(p, s, e)  % Curve length of polynomial y=polyval(p,x) on [s,e]
    if e < s, [s,e] = deal(e,s); end  % Ensure ascending limits
    dp = polyder(p);  % Derivative polynomial coefficients
    ds = @(x) sqrt(1 + (polyval(dp, x)).^2);  % Arc-length integrand sqrt(1 + (dy/dx)^2)
    L = integral(ds, s, e, 'RelTol',1e-8, 'AbsTol',1e-12);  % Numerical integration
end  % End poly_len
