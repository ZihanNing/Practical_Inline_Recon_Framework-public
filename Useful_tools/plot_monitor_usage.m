% Change this to your actual log filename
filename = '//home/zn23/Gadgetron_Parallel_Framework/Useful_tools/monitor_log_PID30575.txt';

% Open the file
fid = fopen(filename);
lines = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);

% Extract lines (ignore header)
lines = lines{1};
lines(1) = [];  % Remove header line

% Preallocate
num_lines = length(lines);
cpu_used_cores = zeros(num_lines, 1);
ram_used_gb = zeros(num_lines, 1);
gpu_mem_mb = nan(num_lines, 1);  % N/A values become NaN

for i = 1:num_lines
    % Split line using "|" as delimiter
    tokens = strsplit(lines{i}, '|');
    
    if length(tokens) < 7
        continue;  % Skip malformed lines
    end

    try
        % Clean up and parse values
        cpu_used_cores(i) = str2double(strtrim(tokens{3}));
        ram_used_gb(i) = str2double(strtrim(tokens{5}));
        
        gpu_mem_str = strtrim(tokens{7});
        gpu_mem_val = str2double(gpu_mem_str);
        if ~isnan(gpu_mem_val)
            gpu_mem_mb(i) = gpu_mem_val;
        end
    catch
        % Handle possible conversion errors gracefully
        cpu_used_cores(i) = NaN;
        ram_used_gb(i) = NaN;
        gpu_mem_mb(i) = NaN;
    end
end

% Create a time vector: 0, 1, ..., N-1 seconds
time_sec = seconds(0:num_lines-1);

% Plotting
figure;

subplot(3,1,1);
plot(time_sec, cpu_used_cores, 'b', 'LineWidth', 2);
ylabel('CPU Used (cores)');
title('CPU Usage Over Time');
grid on;

subplot(3,1,2);
plot(time_sec, ram_used_gb, 'g', 'LineWidth', 2);
ylabel('RAM Used (GB)');
title('RAM Usage Over Time');
grid on;

subplot(3,1,3);
plot(time_sec, gpu_mem_mb, 'r', 'LineWidth', 2);
ylabel('GPU Memory (MB)');
xlabel('Time (s)');
title('GPU Memory Usage Over Time');
grid on;
