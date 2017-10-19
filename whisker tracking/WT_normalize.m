function norm_frame = WT_normalize(frame_in, dim)
%norm_frame = normalize normalize intensity to [0, 1]
%    frame_in is input frame or video
%    dim is dimension flag, optional
%    Jinghao Lu 01/16/2016
%     tic
    if nargin < 2
        dim = 4; %%% treat the tensor as a whole %%%
    end
    if dim == 4
%         disp('Begin normalize')
        mn = nanmin(nanmin(nanmin(frame_in)));
        norm_frame = frame_in - mn;
        clear frame_in
        mx = nanmax(nanmax(nanmax(norm_frame)));
        norm_frame = norm_frame / mx;
    else
        norm_frame = bsxfun(@minus, frame_in, nanmin(frame_in, [], dim));
        norm_frame = bsxfun(@rdivide, norm_frame, (nanmax(norm_frame, [], dim)));
%         norm_frame = bsxfun(@rdivide, norm_frame, (repmat(max(norm_frame, [], dim), 1, size(frame_in, dim))));
    end
%     time = toc;
%     disp(['Done normalize, time: ', num2str(time)])
end