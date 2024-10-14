function A1 = initialize_kernels(Y, num_kernels, kernel_size, kerneltype)
    A1 = cell(1, num_kernels);

    switch lower(kerneltype)
        case 'random'
            for n = 1:num_kernels
                A1{n} = proj2oblique(randn(kernel_size(n,:)));
            end

        case 'selected'
            figure;
            imagesc(Y);
            colorbar;
            axis square;
            title('Select initial kernel regions');

            [img_height, img_width] = size(Y);

            for n = 1:num_kernels
                % Create a draggable rectangle for selection
                h = imrect(gca, [1 1 min(kernel_size(n,2), img_width) min(kernel_size(n,1), img_height)]);
                
                % Set constraints to keep the rectangle within the image bounds
                setPositionConstraintFcn(h, @(pos) constrainPosition(pos, img_width, img_height, kernel_size(n,:)));
                
                position = wait(h);

                % Extract the selected region
                x1 = round(position(1));
                y1 = round(position(2));
                x2 = min(x1 + kernel_size(n,2) - 1, img_width);
                y2 = min(y1 + kernel_size(n,1) - 1, img_height);

                selected_kernel = Y(y1:y2, x1:x2);

                % Project onto the oblique manifold
                A1{n} = proj2oblique(selected_kernel);
            end
            close(gcf);

        otherwise
            error('Invalid kernel initialization type. Choose "random" or "selected".');
    end
end

function new_pos = constrainPosition(pos, img_width, img_height, kernel_size)
    % Constrain the position of the rectangle to stay within the image bounds
    new_pos = pos;
    new_pos(1) = max(1, min(pos(1), img_width - kernel_size(2) + 1));
    new_pos(2) = max(1, min(pos(2), img_height - kernel_size(1) + 1));
    new_pos(3) = kernel_size(2);
    new_pos(4) = kernel_size(1);
end
