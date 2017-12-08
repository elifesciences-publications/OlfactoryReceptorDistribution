function genenv_example(fname, n, beta)
% GENENV_EXAMPLE Generate GIF showing random environment generation.

pause_frames = 10;
cov_frames = 10;
ln_mu = 0;
ln_sigma = 0.8326;
clim = [-2 2];

f = figure;
f.Units = 'pixels';
f.Position = [f.Position(1:2) 100 100];
set(f, 'color', 'white');

% axis tight manual; % this ensures that getframe() returns a consistent size
cmap0 = divergent([0.21 0.17 0.53], [0.98 0.40 0.17], 256);
colormap(cmap0);

axis equal;
axis off;
% ax = gca;

first_frame = true;

% write the animation for generating the correlation matrix
M0 = randcorr(n, beta, 'callback', @callback, 'shuffle', false);

% write some pause frames
for i = 1:pause_frames
    clf;
    imagesc(M0, clim);
    drawnow;
    
    write_frame(f);
end

% bring in the variances
variances = lognrnd(ln_mu, ln_sigma, n, 1);
stdevs = sqrt(variances);

M = diag(stdevs)*M0*diag(stdevs);
for i = 1:cov_frames
    crt_t = (i-1)/(cov_frames-1);
    crt_M = (1 - crt_t)*M0 + crt_t*M;
    clf;
    imagesc(crt_M, clim);
    drawnow;
    
    write_frame(f);
end

% write some pause frames
for i = 1:pause_frames
    clf;
    imagesc(M, clim);
    drawnow;
    
    write_frame(f);
end

    function callback(S)
        clf;
        imagesc(S, clim);
        drawnow;
        
        write_frame(f);
    end

    function [imind, cm] = capture_frame(f)
        % capture image from figure
        frame = getframe(f);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
    end

    function write_frame(f)
        % write image from figure to GIF file
        hold on;
        plot([0.5 n+0.5 n+0.5 0.5 0.5], [0.5 0.5 n+0.5 n+0.5 0.5], 'k');
        axis equal;
        axis off;
        [imind, cm] = capture_frame(f);
        if first_frame
            imwrite(imind, cm, fname, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
            first_frame = false;
        else
            imwrite(imind, cm, fname, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
        end        
    end
  
end