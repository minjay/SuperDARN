activity_data = readtable('OMNI_Feb29_2012.csv');

% downsample data to every 2 mins
activity_data = activity_data(1:2:size(activity_data, 1), :);

% drop hr = 0, min = 0
activity_data = activity_data(2:size(activity_data, 1), :);

% select the first 4 hours
activity_data = activity_data(1:120, :);

% select column 6 (AE-index, [nT])
ae_index = activity_data.x6;