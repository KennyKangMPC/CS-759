function STAT_ENERG_TROP_DATA = AVG_STAT_ENERGY_PLOT
% Generates a plot of the time-averaged statistical energy versus wavelength, as shown in
% Figure 5 of Di Qi's paper. Given the disturbance potential voriticity
% in the barotropic and baroclinic modes, denoted pot_vort_trop and
% pot_vort_clin, we calculate their corresponding streamfunctions in
% Fourier space, denoted strmfunc_trop and strmfunc_clin, and compute the
% arrays containing the barotropic and baroclinic statistical energies
% (given by the expressions in the caption of Figure 5). We then take the
% radial average of these matrices, and plot them.

% Simulation and output file parameters.
grid_size = 256;
def_wavenum = 10;
num_files = 31;
output_freq = 200;

% Frequency-space operators used to obtain the streamfunctions for the 
% disturbance in potential vorticity in the barotropic and baroclinic
% modes.
wavenumbers = [0:grid_size/2 -grid_size/2+1:-1]';
[x_wavenumbers, y_wavenumbers] = meshgrid(wavenumbers, wavenumbers);
freq_deriv_x = 1i*repmat(wavenumbers',[grid_size 1 2]);
freq_deriv_y = 1i*repmat(wavenumbers,[1 grid_size 2]);
freq_laplacian = freq_deriv_x(:,:,1).^2+freq_deriv_y(:,:,1).^2;
inv_freq_trop = 1./freq_laplacian; inv_freq_trop(1,1) = 0;
inv_freq_clin = 1./(freq_laplacian-def_wavenum^2); inv_freq_clin(1,1) = 0;

% Create arrays for storing the time-averaged statistics.
stat_energ_trop = zeros(grid_size/2+1,1);
stat_energ_clin = zeros(grid_size/2+1,1);

% Loop through all available output files.
for file_num = 15:num_files-1
    lay1_file_name = sprintf('layer1_%08d.csv', output_freq*file_num);
    lay2_file_name = sprintf('layer2_%08d.csv', output_freq*file_num);
    
    
    % Read the files containing the disturbance in potential vorticity in the
    % barotropic and baroclinic modes.
    pot_vort_lay1 = dlmread(lay1_file_name);
    pot_vort_lay2 = dlmread(lay2_file_name);
    pot_vort_lay1 = pot_vort_lay1(:,1:end-1);
    pot_vort_lay2 = pot_vort_lay2(:,1:end-1);
    
    pot_vort_trop = 0.5 * (pot_vort_lay1 + pot_vort_lay2);
    pot_vort_clin = 0.5 * (pot_vort_lay1 - pot_vort_lay2);
    
    % Calculate the corresponding steamfunctions in Fourier Space.
    pot_vort_trop = fft2(pot_vort_trop);
    pot_vort_clin = fft2(pot_vort_clin);
    
    strmfunc_trop = inv_freq_trop.*pot_vort_trop;
    strmfunc_clin = inv_freq_clin.*pot_vort_clin;
    
    % Calculate the radial average of the statisical energy of the barotropic
    % and baroclinic modes.
    for i = 1:grid_size
        for j = 1:grid_size
            wavenum = sqrt(x_wavenumbers(i,j)^2 + y_wavenumbers(i,j)^2);
            if ceil(wavenum) <= grid_size/2
                radius_bin = wavenum - floor(wavenum);
                stat_energ_trop(floor(wavenum)+1) = stat_energ_trop(floor(wavenum)+1) ...
                    + (1-radius_bin)*(wavenum^2)*abs(strmfunc_trop(i,j))^2;
                stat_energ_trop(ceil(wavenum)+1) = stat_energ_trop(ceil(wavenum)+1) ...
                    + radius_bin*(wavenum^2)*abs(strmfunc_trop(i,j))^2;
                stat_energ_clin(floor(wavenum)+1) = stat_energ_clin(floor(wavenum)+1)...
                    + (1-radius_bin)*(wavenum^2 + def_wavenum^2)*abs(strmfunc_clin(i,j))^2;
                stat_energ_clin(ceil(wavenum)+1) = stat_energ_clin(ceil(wavenum)+1)...
                    + radius_bin*(wavenum^2 + def_wavenum^2)*abs(strmfunc_clin(i,j))^2;
            end
        end
    end
    
end

stat_energ_trop = 0.5*stat_energ_trop/(grid_size^4 * num_files);
stat_energ_clin = 0.5*stat_energ_clin/(grid_size^4 * num_files);

STAT_ENERG_TROP_DATA = stat_energ_trop(2:129);
stat_energ_clin_data = stat_energ_clin(2:129);


% Plot the energy matrix for the barotropic part.
%h = pcolor(stat_energ_trop(:,:));
%set(h, 'EdgeColor', 'none');
%disp(stat_energ_trop(1:5,1:5));

% Plot the radial average for the barotropic part.
plot_wavenums = linspace(1,128,128)';
interp_wavenums = linspace(1,128,500)';
  % Interpolate, smooth, then plot.
intperp_trop = interp1(log(plot_wavenums),log(STAT_ENERG_TROP_DATA),log(interp_wavenums));
smth_trop = smooth(intperp_trop);
for i = 1:20
    smth_trop = smooth(smth_trop);
end 
trop_energ = loglog(interp_wavenums, exp(smth_trop),'k-','DisplayName','Barotropic');

hold on
  % Plot linear best fit for wavenumbers 10 - 100.
best_x = plot_wavenums(10:100);
poly_trop = polyfit(log(best_x),log(STAT_ENERG_TROP_DATA(10:100)),1);
best_y_trop = polyval(poly_trop, log(best_x));
loglog(best_x,exp(1).*exp(best_y_trop),'k--');
disp(["Barotropic Best-Fit Slope" poly_trop(1)]);

% Plot the radial average for the baroclinic part.
  % Interpolate, smooth, then plot.
intperp_clin = interp1(log(plot_wavenums),log(stat_energ_clin_data),log(interp_wavenums));
smth_clin = smooth(intperp_clin);
for i = 1:20
    smth_clin = smooth(smth_clin);
end 
clin_energ = loglog(interp_wavenums, exp(smth_clin),'k-.','DisplayName','Baroclinic');
  % Plot linear best fit for wavenumbers 10 - 100.
best_x = plot_wavenums(10:100);
poly_clin = polyfit(log(best_x),log(stat_energ_clin_data(10:100)),1);
best_y_clin = polyval(poly_clin, log(best_x));
loglog(best_x,exp(1).*exp(best_y_clin),'k--');
disp(["Baroclinic Best-Fit Slope" poly_clin(1)]);

% Label the plot.
title("Statistical Energy (Atmosphere, Mid-Latitude)")
xlabel("Wavenumber |k|")
legend([trop_energ clin_energ])

% Change plot axes.
xlim([1,200])
ylim([10^(-6), 4])
hold off
