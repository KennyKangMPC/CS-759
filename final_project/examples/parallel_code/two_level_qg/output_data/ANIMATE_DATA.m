function ANIMATE_DATA(dim1_len, dim2_len, file_count, output_freq, level)
% Creates an animation of the output data.
%   ARGUMENTS: - dim1_len, dim2_len: The dimensions of the grid in the
%   Fortran code. Note that MATLAB is row-major and not column major, so
%   dim1_len is the length along the second index and vice-versa in MATLAB.
%   - num_files: The number of output files.
%   - output_freq: The output frequency of the simulation.

% We will store all of the output in a single array, and call them for when
% we plot them.
all_data = zeros(dim2_len, dim1_len , file_count);
% Fill in all_data.
if level == 1
    for k = 0:file_count-1
        csv_file_name=sprintf('level1_%08d.csv', output_freq*k);
        csvdata = csvread(csv_file_name);
        all_data(:,:,k+1) = csvdata(:,1:dim1_len);
    end
elseif level == 2
    for k = 0:file_count-1
        csv_file_name=sprintf('level2_%08d.csv', output_freq*k);
        csvdata = csvread(csv_file_name);
        all_data(:,:,k+1) = csvdata(:,1:dim1_len);
    end
end


% Animate the plot by rewriting it repeatedly.
k = 1;
while k <= file_count
    h = pcolor(all_data(:,:,k));
    set(h, 'EdgeColor', 'none');
    colorbar;
    %caxis([-2 2]);
    title(join(["Level ", int2str(level), " - ", int2str(output_freq*(k-1))]));
    pause(.1);
    k = k + 1;
    if k == file_count
        k = 1;
    end
end
end