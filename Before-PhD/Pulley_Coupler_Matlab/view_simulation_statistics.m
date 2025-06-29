addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'Matlab Fcn Lib')));\n% Define the root directory
root_dir = 'C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\PG_Coupler_Jan_9_2025\log';

% Get all subdirectories
folders = dir(fullfile(root_dir, 'Ring*'));
folders = folders([folders.isdir]);

% Initialize arrays to store data
body_widths = [];
radii = [];
momentum_diffs = struct();

% Loop through each subdirectory
for i = 1:length(folders)
    folder_name = folders(i).name;
    folder_path = fullfile(root_dir, folder_name);
    
    % Extract body waveguide width and radius from the folder name
    tokens = regexp(folder_name, 'Ring([\d\.]+)um_([\d\.]+)x2\.00um', 'tokens');
    if isempty(tokens)
        continue;
    end
    radius = str2double(tokens{1}{1});
    body_width = str2double(tokens{1}{2});

    % Locate the @1064.0nm_*.mat file
    mat_file = dir(fullfile(folder_path, '@532.0nm_*.mat'));
    if isempty(mat_file)
        continue;
    end
    
    % Load the results struct
    mat_data = load(fullfile(folder_path, mat_file.name));
    results = mat_data.results;
    
    % Store body width and radius
    body_widths(end + 1) = body_width;
    radii(end + 1) = radius;

    % Extract body_TE value (scalar)
    body_TE = results.body_TE;
    
    % Process each bus width key (w_X)
    field_names = fieldnames(results);
    for j = 1:length(field_names)
        field = field_names{j};
        if startsWith(field, 'w')
            if isfield(results.(field), 'bus_TE')
                bus_TE = results.(field).bus_TE;
                min_diff = min(abs(bus_TE - body_TE));
                
                % Store results using the original field name (with underscores)
                if ~isfield(momentum_diffs, field)
                    momentum_diffs.(field) = [];
                end
                momentum_diffs.(field) = [momentum_diffs.(field); body_width, radius, min_diff];
            end
        end
    end
end

% Plot the results for each w_X with correctly formatted labels
figure;
hold on;
fields = fieldnames(momentum_diffs);
colors = lines(length(fields));
legend_entries = cell(length(fields), 1);

for k = 1:length(fields)
    data = momentum_diffs.(fields{k});
    scatter3(data(:,1), data(:,2), data(:,3), 50, colors(k,:), 'filled');
    
    % Convert field name (w_2_5 to 2.5) for labeling
    numeric_label = strrep(fields{k}(2:end), '_', '.');
    legend_entries{k} = ['\beta - \beta_{Body TE} for w = ', numeric_label];
end

% Label the plot
xlabel('Body Waveguide Width (\mum)');
ylabel('Radius (\mum)');
zlabel('Minimum Momentum Mismatch');
title('Momentum Mismatch Analysis');
legend(legend_entries, 'Location', 'best');
grid on;
view(135, 30);
hold off;


%%
% 
% Initialize an array to store the overall minimum momentum mismatch for each (x, y)
unique_data = [];

% Loop through all stored momentum differences to find the minimum for each (x, y)
for k = 1:length(fields)
    data = momentum_diffs.(fields{k});
    
    % Concatenate all (x, y, momentum) values into a single matrix
    unique_data = [unique_data; data];
end

% Find unique (x, y) pairs and their minimum z values
[x_unique, ~, idx] = unique(unique_data(:,1:2), 'rows');
min_z = accumarray(idx, unique_data(:,3), [], @min);

% Combine the results into a single matrix
final_results = [x_unique, min_z];

% Plot the overall minimum momentum mismatch for each (x, y)
figure;
scatter3(final_results(:,1), final_results(:,2), final_results(:,3), 50, 'filled');
hold on;
fill3([min(final_results(:,1)) max(final_results(:,1)) max(final_results(:,1)) min(final_results(:,1))], ...
      [min(final_results(:,2)) min(final_results(:,2)) max(final_results(:,2)) max(final_results(:,2))], ...
      [0 0 0 0], 'k', 'FaceAlpha', 0.3);
hold off;
% Label the plot
xlabel('Body Waveguide Width (\mum)');
ylabel('Radius (\mum)');
zlabel('Minimum Momentum Mismatch');
title('Overall Minimum Momentum Mismatch');
zlim([-2 20]);
grid on;
view(135, 30);
colorbar;
