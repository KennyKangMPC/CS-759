function ANIMATE_DATA(dim1_len, dim2_len, file_count, output_freq)
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
for k = 0:file_count-1
    csv_file_name=sprintf('out_%08d.csv', output_freq*k);
    csvdata = csvread(csv_file_name);
    all_data(:,:,k+1) = csvdata(:,1:dim1_len);
end

% Animate the plot by rewriting it repeatedly.
k = 1;
while k <= file_count
    h = pcolor(all_data(:,:,k));
    set(h, 'EdgeColor', 'none');
    colorbar;
    colormap(gray);
    caxis([65 67.5]);
    title(int2str(output_freq*(k-1)));
    pause(.1);
    k = k + 1;
    if k == file_count
        k = 1;
    end
end
end

