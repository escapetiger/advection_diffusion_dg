function error_table = convergence(errors, resolutions, error_descriptions, error_names, desired_rates)
%CONVERGENCE Calculate convergence orders and generate an error table.
%
% Inputs:
%   errors            - Cell matrix where each row corresponds to a type of error
%                       (e.g., 'L1', 'L2', 'Linf'), and each column corresponds to
%                       a refinement level.
%   resolutions       - Array of scalar resolutions for convergence calculations
%                       (e.g., [16, 32, 64]).
%   error_descriptions - Cell array of strings describing each resolution
%                        (e.g., {'16x16', '32x32', '64x64'}).
%   error_names       - Cell array of strings for naming each error type
%                       (e.g., {'L1', 'L2', 'Linf'}). If empty, defaults to
%                       'Error_Type_1', 'Error_Type_2', etc.
%   desired_rates     - (Optional) Array of desired convergence rates for each error type.
%
% Output:
%   error_table - Table of errors and convergence orders for all types of errors.

% Validate input dimensions
[num_types, num_levels] = size(errors);
if num_levels ~= length(resolutions)
    error('Number of refinement levels must match the number of columns in errors.');
end
if nargin < 3 || isempty(error_descriptions)
    % Default to string versions of scalar resolutions
    error_descriptions = arrayfun(@(r) sprintf('%dx', r), resolutions, 'UniformOutput', false);
end
if length(error_descriptions) ~= num_levels
    error('Number of resolution descriptions must match the number of levels in resolutions.');
end
if nargin < 4 || isempty(error_names)
    % Default to generic error names if not provided
    error_names = arrayfun(@(i) sprintf('Error_Type_%d', i), 1:num_types, 'UniformOutput', false);
elseif length(error_names) ~= num_types
    error('Number of error names must match the number of error types (rows in errors).');
end
if nargin < 5 || isempty(desired_rates)
    desired_rates = zeros(1, num_types); % Default: no desired rates
elseif length(desired_rates) ~= num_types
    error('Number of desired rates must match the number of error types (rows in errors).');
end

% Initialize arrays for errors and convergence orders
error_values = zeros(num_types, num_levels);
error_orders = zeros(num_types, num_levels);

% Extract errors from the cell matrix
for i = 1:num_types
    for j = 1:num_levels
        error_values(i, j) = errors{i, j};
    end
end

% Calculate convergence orders
for i = 1:num_types
    for j = 2:num_levels
        h_ratio = resolutions(j-1) / resolutions(j); % Refinement ratio
        if error_values(i, j) > 0 && error_values(i, j-1) > 0
            % Compute convergence order
            error_orders(i, j) = log(error_values(i, j-1) / error_values(i, j)) / log(h_ratio);
        else
            % Leave as zero if invalid
            error_orders(i, j) = 0;
        end
    end
end

% Prepare table data
table_data = cell(num_levels, 2 * num_types + 1);
table_data(:, 1) = error_descriptions(:); % Resolution descriptions
for i = 1:num_types
    table_data(:, 2*i) = num2cell(error_values(i, :).');
    table_data(:, 2*i+1) = num2cell(error_orders(i, :).');
end

% Create table
var_names = {'Resolution'};
for i = 1:num_types
    var_names{end+1} = sprintf('%s_Error', error_names{i});
    var_names{end+1} = sprintf('%s_Order', error_names{i});
end
error_table = cell2table(table_data, 'VariableNames', var_names);

% % Plot errors on a log-log scale
% figure;
% hold on;
% colors = lines(num_types); % Generate distinct colors for each error type
% 
% % Loop through each type of error
% for i = 1:num_types
%     % Plot error values
%     loglog(resolutions, error_values(i, :), '-o', 'Color', colors(i, :), ...
%         'DisplayName', sprintf('%s Error', error_names{i}), ...
%         'LineWidth', 2);
% 
%     % Plot desired convergence rate (dashed line)
%     if desired_rates(i) > 0
%         % Reference line: scaling as h^(desired_rate)
%         ref_slope = error_values(i, 1) * (resolutions / resolutions(1)).^(desired_rates(i));
%         loglog(resolutions, ref_slope, '--', 'Color', colors(i, :), ...
%             'DisplayName', sprintf('%s Desired Rate (%d)', error_names{i}, desired_rates(i)), ...
%             'LineWidth', 2);
%     end
% end
% 
% % Set log-log scale
% set(gca, 'XScale', 'log', 'YScale', 'log');
% xlabel('Resolution');
% ylabel('Error');
% legend('show', 'Location', 'Best');
% grid on;
% title('Error Convergence (Log-Log Scale)');
% hold off;
end
