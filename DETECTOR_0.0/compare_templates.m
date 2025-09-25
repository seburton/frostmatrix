% produces a similarity analysis of the template set
clear
  % Define the maximum range of station allowed from the lava lake center for the analysis
  range_thresh = 3;

  % Load Erika Jaski's list of 17 Eruption templates from her MS Thesis
A = importdata('template_list.txt');

% Determine the number of templates loaded
ntemplates = length(A);
for i = 1:ntemplates
    % Store the template name and load the corresponding data
    template_name_multi{i} = char(A{i});
    eval(['load ', char(template_name_multi{i})]);
    template_multi{:,:,i} = template; % Store the template data
    template_sachdr_multi{i} = template_sachdr; % Store the template header information
    template_range_multi{i} = template_range; % Store the template range information
end
disp([num2str(ntemplates), ' templates loaded'])

kcount = 1; % Initialize a counter for matching templates
for i = 1:ntemplates
    % Convert the template range information from cell to matrix
    ranges = cell2mat(template_range_multi(i));
    nchans1 = length(template_sachdr_multi{i}); % Number of channels in the first template
    sachdr1 = template_sachdr_multi{i}; % Header information for the first template
    
    % Compare the first template with all other templates
    for j = 1:ntemplates
        nchans2 = length(template_sachdr_multi{j}); % Number of channels in the second template
        sachdr2 = template_sachdr_multi{j}; % Header information for the second template
        
        for k = 1:nchans1
            % Find common station and component indices for the current pair of templates
            kout = find(strcmp(sachdr1(k).kstnm, {sachdr2.kstnm}) & strcmp(sachdr1(k).kcmpnm, {sachdr2.kcmpnm}));
            if ~isempty(kout)
                kfnd(kcount) = kout; % Store the index of the matching component
                
                % Save matching information: template indices, channel index, and found index
                kstore(kcount, :) = [i, j, k, kfnd(kcount), ranges(k)];
                
                % Increment the counter for matches found
                kcount = kcount + 1;
            end
        end
        
        % Perform cross-correlation between the two templates
        x1 = cell2mat(template_multi(1, 1, i)); % Data for the first template
        x2 = cell2mat(template_multi(1, 1, j)); % Data for the second template
        
        % Find indices of stored matches that meet the range threshold
        kk = find(kstore(:, 1) == i & kstore(:, 2) == j & kstore(:, 5) <= range_thresh);
        ind1 = kstore(kk, 3); % Indices for the first template
        ind2 = kstore(kk, 4); % Indices for the second template
        
        % Filter out any negative indices
        ind1 = ind1(ind1 > 0);
        ind2 = ind2(ind2 > 0);
        
        % Extract the relevant data segments for cross-correlation
        y1 = x1(:, ind1);
        y2 = x2(:, ind2);
        
        % Reshape the data for cross-correlation
        [a, b] = size(y1);
        y1 = reshape(y1, a * b, 1);
        y2 = reshape(y2, a * b, 1);
        
        % Compute the normalized cross-correlation
        [a, b] = xcorr(y1, y2, 'normalized');
        [amax, IA] = max(a); % Find the maximum correlation value and its index
        
        % Store the maximum correlation value and the corresponding lag
        C(i, j) = max(a); % Store the maximum correlation
        D(i, j) = b(IA) / 200; % Store the lag in seconds
        N(i, j) = length(ind1); % Store the number of channels used for the correlation
    end
end

% Plot the results of the similarity analysis
figure(1000)
colormap(redbluecmap)

% Plot the correlation matrix
subplot(1, 3, 1)
imagesc(C)
clim([-1 1]) % Set color limits for the correlation plot
h = colorbar;
h.Label.String = 'Correlation'; % Label for the colorbar
xlabel('Template Index') % X-axis label
ylabel('Template Index') % Y-axis label
axis square % Set axis to be square
bookfonts_TNR % Apply font settings

% Plot the number of channels used for each correlation
subplot(1, 3, 2)
imagesc(N)
h = colorbar;
h.Label.String = 'Number of Channels'; % Label for the colorbar
xlabel('Template Index') % X-axis label
ylabel('Template Index') % Y-axis label
axis square % Set axis to be square
bookfonts_TNR % Apply font settings

% Plot the correlation lag
subplot(1, 3, 3)
imagesc(D)
h = colorbar;
h.Label.String = 'Correlation Lag (s)'; % Label for the colorbar
xlabel('Template Index') % X-axis label
ylabel('Template Index') % Y-axis label
axis square % Set axis to be square
bookfonts_TNR % Apply font settings
