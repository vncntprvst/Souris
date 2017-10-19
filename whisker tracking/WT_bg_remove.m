function Ydebg = WT_bg_remove(Y, sz)
%deb = bg_remove(tt) remove background
%    Remove background use difference between imopen and original sequence
%    Y is input sequence; blurred is preferred
%    Ydebg is output bg removed sequence
%    use imopen with strel 'disk' with 15
%    updated version: 'disk' with 9; and can use anisotropic diffusion
%    Jinghao Lu 01/16/2016

    h = tic;
    disp(['Begin bg remove'])
    if isempty(gcp('nocreate'))
        parpool(feature('numCores'));
    end
    if nargin < 2
        sz = 9;
    end
    nframes = size(Y, 3);
    Ydebg = zeros(size(Y));
    k = strel('disk', sz);
    parfor i = 1: nframes
        I = Y(:, :, i);
        bg = imopen(I, k);
        tmp = I - bg;
        tmp(1: sz, 1: sz) = tmp(sz + 1, sz + 1);
        tmp(1: sz, end - sz: end) = tmp(sz + 1, end - sz - 1);
        tmp(end - sz: end, 1: sz) = tmp(end - sz - 1, end - sz - 1);
        tmp(end - sz: end, end - sz: end) = tmp(end - sz - 1, end - sz - 1);        
        Ydebg(:, :, i) = tmp;
        if mod(i, 100) == 0
            disp(['Done #', num2str(i), '/', num2str(nframes)])
        end
    end    
    
    %%% norm after remove %%%
    Ydebg = WT_normalize(Ydebg);
    time = toc(h);
    disp(['Done bg remove, time: ', num2str(time)])
end