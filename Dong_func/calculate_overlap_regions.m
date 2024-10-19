function IDX = calculate_overlap_regions(Xiter, k)
    % this function calculate the overlap regions between the activation regions
    % Xiter: the activation regions
    % k: the size of the kernel
    % IDX: the overlap region, size(IDX) = size(Xiter)

    [h, w, num_kernels] = size(Xiter);
    IDX = false(size(Xiter));
    Mask = false(size(Xiter));
    % create square mask around every activation region with the size of its kernel by convolving the activation region with a square kernel
    for n = 1:num_kernels
        Mask(:,:,n) = conv2(Xiter(:,:,n), ones(k(n,1), k(n,2)), 'same') > 0;
    end
    % calculate the overlap region by and operation between target mask and other masks
    for n = 1:num_kernels
        for m = 1:num_kernels
            if n ~= m
                IDX(:,:,n) = IDX(:,:,n) | (Mask(:,:,n) & Mask(:,:,m));
            end
        end
    end
end