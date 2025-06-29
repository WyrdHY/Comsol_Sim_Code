split_functions("C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Code\Matlab Fcn Lib\Fcn_Lib.m");
function split_functions(fcn_lib_path)
    % Open the file for reading
    fid = fopen(fcn_lib_path, 'r');
    if fid == -1
        error('Could not open the file: %s', fcn_lib_path);
    end

    file_contents = fread(fid, '*char')';
    fclose(fid);

    % Regular expression to match MATLAB function definitions including multi-line signatures
    func_pattern = '(?<=\n|^)\s*function\s+([\s\S]+?)\n([\s\S]*?)(?=\n\s*function\s+|$)';

    % Find all matches
    funcs = regexp(file_contents, func_pattern, 'tokens');

    if isempty(funcs)
        error('No functions found in the file.');
    end
    

    
    % Get the full path of the currently executing script
    if isdeployed
        script_path = ctfroot;  % For deployed applications
    else
        script_path = matlab.desktop.editor.getActiveFilename;
    end
    
    % Extract the directory where the script is located
    script_folder = fileparts(script_path);
    
    % Define the target folder within the same directory as Partition.m
    target_dir = fullfile(script_folder, 'FCN_Lib');
    
    % Create the folder if it does not exist
    if ~exist(target_dir, 'dir')
        mkdir(target_dir);
        disp(['Folder created: ', target_dir]);
    else
        disp(['Folder already exists: ', target_dir]);
    end

    % Process each function found
    for i = 1:length(funcs)
        func_signature = strtrim(funcs{i}{1});
        func_body = strtrim(funcs{i}{2});

        % Extract function name by considering different possible signatures
        func_name_match = regexp(func_signature, '^\s*\[?.*?\]?\s*=\s*([a-zA-Z_]\w*)', 'tokens', 'once');
        
        if isempty(func_name_match)
            func_name_match = regexp(func_signature, '^\s*([a-zA-Z_]\w*)', 'tokens', 'once');
        end
        
        if isempty(func_name_match)
            warning('Could not determine function name for:\n%s\nSkipping...', func_signature);
            continue;
        end
        
        % Use the function name found in the signature
        func_name = strtrim(func_name_match{1});

        % Write each function to a new .m file in the Fcn_lib folder
        new_file_name = fullfile(target_dir, [func_name, '.m']);
        fid_new = fopen(new_file_name, 'w');
        if fid_new == -1
            warning('Could not create file: %s', new_file_name);
            continue;
        end

        fprintf(fid_new, 'function %s\n%s\n', func_signature, func_body);
        fclose(fid_new);

        fprintf('Created: %s\n', new_file_name);
    end

    fprintf('All functions successfully split into the Fcn_lib folder.\n');
end
