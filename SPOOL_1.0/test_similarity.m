
% Define the names of the MAT-files to compare
file1 = '/Users/sburton/Library/CloudStorage/OneDrive-Colostate/Work/SPOOL_1.0/Test/Templates/template_345_00000_888_908.mat';
file2 = '/Users/sburton/Work/Extraction_Scripts/Test/template_2008-12-10T12:21:46_data.mat';

% Use visdiff to compare the two MAT-files
visdiff(file1, file2);

data1 = load(file1);
data2 = load(file2);

% Check if the structures are the same
if isequal(data1, data2)
    disp('The two MAT-files are identical.');
else
    disp('The two MAT-files are different.');
end

% Find common variables
commonVars = intersect(fieldnames(data1), fieldnames(data2));

% Display common variables
if isempty(commonVars)
    disp('No common data found between the two .mat files.');
else
    disp('Common data found:');
    disp(commonVars);

    % Check for common values in common variables
    for i = 1:length(commonVars)
        varName = commonVars{i};
        values1 = data1.(varName);
        values2 = data2.(varName);

        % Check if the values are the same
        if isequal(values1, values2)
            fprintf('Values in variable "%s" are identical.\n', varName);
        else
            fprintf('Values in variable "%s" are different.\n', varName);
            % Find common values
            % Check if values are numeric before finding common values
            if isnumeric(values1) && isnumeric(values2)
                commonValues = intersect(values1(:), values2(:)); % Flatten and find common values
                if ~isequal(values1, values2)
                    if numel(values1) == numel(values2)
                        avgDifference = mean(abs(values1(:) - values2(:))); % Calculate average difference
                        fprintf('Average difference in variable "%s": %.4f\n', varName, avgDifference);
                    else
                        fprintf('Numeric values in variable "%s" have different lengths and cannot be compared.\n', varName);
                    end
                end
            elseif isstring(values1) && isstring(values2)
                commonValues = intersect((values1(:)), values2(:)); % Convert to strings and find common values
            elseif iscell(values1) && iscell(values2)
                commonValues = intersect(values1, values2, 'stable'); % Find common values in cell arrays
                if ~isequal(values1, values2)
                    if numel(values1) == numel(values2)
                        fprintf('Cell arrays in variable "%s" have different lengths and cannot be compared.\n', varName);
                    end
                end
            elseif isstruct(values1) && isstruct(values2)
                commonFields = intersect(fieldnames(values1), fieldnames(values2));
                commonValues = [];
                for j = 1:length(commonFields)
                    fieldName = commonFields{j};
                    if isequal(values1.(fieldName), values2.(fieldName))
                        commonValues = [commonValues; values1.(fieldName)]; % Collect common values
                    end
                end
            else
                commonValues = []; % Default to empty if types do not match
            end
            if isempty(commonValues)
                fprintf('No common values found in variable "%s".\n', varName);
            else
                fprintf('%i common values in variable "%s". ',length(commonValues), varName);
                %disp(commonValues);
            end
        end
    end
end

