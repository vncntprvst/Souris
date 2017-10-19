function YDeN = WT_anidenoise(Y, ispara, iter, dt, kappa, opt)
%Yblur = anidenoise(Y) denoise using anisotropic diffusion
%    denoise image/sequence
%    Y is input image/sequence
%    Yblur is output image/sequence
%    Jinghao Lu 05/17/2016

    h = tic;
    disp(['Begin anisotropic diffusion denoising'])
    nframes = size(Y, 3);
    YDeN = zeros(size(Y));
    %%% current init; need revising %%%
    if nargin < 2 || isempty(ispara)
        ispara = 0;
    end
    
    if nargin < 3 || isempty(iter)
        iter = 10; %%% original 10 %%%
    end
    
    if nargin < 4 || isempty(dt)
        dt = 0.05;
    end
    
    if nargin < 5 || isempty(kappa)
        kappa = 0.5;
    end
    
    if nargin < 6 || isempty(opt)
        opt = 1;
    end
    
    if nframes == 1
        YDeN = WT_anisodiff_JL(Y, iter, dt, kappa, opt);   
    else
        if isempty(gcp('nocreate'))
            parpool(feature('numCores'));
        end
        parfor i = 1: nframes
            I = Y(:, :, i);
            YDeN(:, :, i) = WT_anisodiff_JL(I, iter, dt, kappa, opt);
%             disp(num2str(i))
        end
    end
    
    %%% norm after remove %%%
    YDeN = WT_normalize(YDeN);
    time = toc(h);
    disp(['Done anisotropic diffusion denoising: ', num2str(time)])
end