function [xi, ti, fi] = line_slice_intersect(x1, t1, x2, t2, a, d)
%LINE_SLICE_INTERSECT Find intersections of a line with a slice in nD+1 space.
%
%   Syntax:
%      [xi, ti, fi] = line_slice_intersect(x1, t1, x2, t2, a, d)
%
%   Inputs:
%      x1 - Space coordinates of the first points. 
%      t1 - Time coordinates of the first points.
%      x2 - Space coordinates of the second points.
%      t2 - Time coordinates of the second points.
%      a  - Scalar value at which the d-th entry of x is set.
%      d  - Index of the dimension to slice on.
%
%   Outputs:
%      xi - Space coordinates of valid points of intersection.
%      ti - Time coordinates of valid points of intersection.
%      fi - Logical array indicating invalid intersections.
%
%   Author: Yi Cai
%   Date: 2024-12-11
%   Version: 1.1

    if ~isscalar(t1) || ~isscalar(t2)
        error('t1 and t2 must be scalar.');
    end

    input_format = detect_input_format(x1);

    % Normalize inputs to column-wise format for computation
    x1 = normalize_input(x1);
    x2 = normalize_input(x2);
    t1 = t1(:);
    t2 = t2(:);
    
    % Compute s for all points
    h = x2(:, d) - x1(:, d); % Change along dimension d
    valid_mask = abs(h) > 1e-10;    % Non-parallel lines
    s = (a - x1(:, d)) ./ h;    % Compute s (NaN for invalid)
    
    % Compute intersection points and times
    xi = x1 + (x2 - x1) .* s; % Intersect points
    ti = t1 + (t2 - t1) .* s; % Intersect times
    
    % Check valid time ranges
    t_min = min(t1, t2);
    t_max = max(t1, t2);
    time_mask = (ti >= t_min) & (ti <= t_max);
    
    % Combine masks
    valid_mask = valid_mask & time_mask;
    fi = ~valid_mask; % Logical array of invalid intersections
    
    % Extract valid intersections
    xi = xi(valid_mask, :);
    ti = ti(valid_mask);
    
    % Format outputs to match the input format
    xi = format_output(xi, input_format);
end

function input_format = detect_input_format(x)
%DETECT_INPUT_FORMAT Detect the format of the input x.
    if iscell(x)
        input_format = 'cell';
    elseif ismatrix(x) && size(x, 2) > 1
        input_format = 'matrix';
    else
        input_format = 'vector';
    end
end

function x_normalized = normalize_input(x)
%NORMALIZE_INPUT Normalize input x to column-wise format.
    if iscell(x)
        x_normalized = cell2mat(x);
    elseif ismatrix(x) && size(x, 2) > 1
        x_normalized = x;
    else
        x_normalized = x(:);
    end
end

function x_formatted = format_output(x, input_format)
%FORMAT_OUTPUT Convert x back to the original input format.
    if strcmp(input_format, 'cell')
        if isvector(x)
            x_formatted = num2cell(x(:));
        else
            x_formatted = num2cell(x, 1);
        end
    elseif strcmp(input_format, 'matrix')
        x_formatted = x;
    else
        x_formatted = x(:)';
    end
end
