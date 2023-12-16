%% Open the recording files
% Make sure the current folder contains the files you want to open and
% analyse, since it is written to find all .abf files in current directory.
datafiles = dir ('*abf'); %stores all infos of the .abf in the current directory
Names = {datafiles(:).name}'; %cell containing the names of the .abf
filename = string(Names);
file_num = size(filename, 1);
table_index = NaN(); %matrix to be filled at the end of the loop
data_conc = NaN(); %matrix to be filled with data satisfying condition at line 68
sfreq = 10000; %sampling freq in Hz - HARDCODED
PSDs = NaN (4096, file_num); %HARDCODED

%% Waves definition
labels = ["Slow oscillation" "Delta" "Theta" "Alpha" "Beta" "Gamma"];
slow = [0.5 1];
delta = [2 4];
theta = [4 8];
alpha = [8 14];
beta = [14 30];
gamma = [30 90];
norm_powers = NaN ();
HLF_Ratio = NaN ();

%% Loop for each file found in the directory
for i = 1:file_num
    % Loading files in the directory
    data = abfload(filename(i));
    
    % Definition of the number of chunks in which each recording is divided
    bin = 60; %bin length in seconds
    l_bin = bin*sfreq; % bin length in points
    n_bin = ceil(length (data)/ l_bin); % number of bins in which the recording is divided into
    segm_data = NaN(l_bin, n_bin); %each recording is segmented into a bin-length by number-of-bins matrix
    
    % Filling the segm_data matrix with each chunk
    if length(data) > l_bin
        for k = 1:n_bin
            if k == max(n_bin)
                segm_data (1:length(data)-((k-1)*l_bin),k) = data(((k-1)*l_bin)+1:length(data));
            else
                segm_data (:,k) = data (((k-1)*l_bin)+1:k*l_bin);
            end
        end
    else
        segm_data = data;
    end
    
    % Loop for each chunk
    for j = 1:n_bin
        if j == max(n_bin)
            segm_data = segm_data(1:length(data)-((j-1)*l_bin), n_bin);
            high_f = 100; %lowpass freq
            data_filt = lowpass (segm_data(:,1), high_f, sfreq); %bandpass filter
        else
            % Lowpass filter
            high_f = 100; %lowpass freq
            data_filt = lowpass (segm_data(:,j), high_f, sfreq); %bandpass filter
        end
        
        % Downsample
        num = 1;
        den = 10; %resampling defined as num/den
        resamp_data = resample (data_filt(:,1), num, den);
        n_sf = sfreq*(num/den); %downsampled frequency
        
        % Power spectrum
        [p,f] = pspectrum (resamp_data, n_sf);
        p_area = trapz(f, p); %Calculate the are below the curve
        norm_p = p/p_area; %normalize the powers by the area below the curve
        
        % Save normalized band powers and Synchronization index
        highsynch_p = bandpower(norm_p, f, [0.1 4], 'psd');
        lowsynch_p = bandpower(norm_p, f, [4 100], 'psd');
        table_index (i,j) = highsynch_p/lowsynch_p; %Synchronization index
        
        % If the bin has a synchronization index above a certain threshold, let's
        %concatenate it with the previous (if there is one) satisfying this
        %condition
        if highsynch_p/lowsynch_p < 1.5
            data_conc ((length(resamp_data)*(j-1))+1:j*length(resamp_data),1) = resamp_data;
        end
    end
    
    % Power spectrum of the concatenated chunks
    [WholeP,WholeF] = pspectrum (data_conc, n_sf);
    WholePowArea = trapz(WholeF, WholeP); %Riemann's definite integral
    N_p = WholeP/WholePowArea; %normalized powers
    plot (WholeF, N_p)
    xlim([0 100])
    ylabel("Relative power spectrum (dB)")
    xlabel("Frequency (Hz)")
    hold on
    PSDs (:,i) = N_p;
        
    % Normalized band powers of the concatenated chunks
    norm_powers (i,1) = bandpower(N_p, WholeF, slow, 'psd');
    norm_powers (i,2) = bandpower(N_p, WholeF, delta, 'psd');
    norm_powers (i,3) = bandpower(N_p, WholeF, theta, 'psd');
    norm_powers (i,4) = bandpower(N_p, WholeF, alpha, 'psd');
    norm_powers (i,5) = bandpower(N_p, WholeF, beta, 'psd');
    norm_powers (i,6) = bandpower(N_p, WholeF, gamma, 'psd');
    LowF_Pow = bandpower(N_p, WholeF, [0.5 7.5], 'psd');
    HighF_Pow = bandpower(N_p, WholeF, [7.5 20.5], 'psd');
    HLF_Ratio(i) = HighF_Pow/LowF_Pow;
end

%Plot the mean PSD spectrum
meanPSD = mean (PSDs, 2);
SD_mean = std (PSDs, 0, 2);
SEM = SD_mean/sqrt(file_num);
plot (WholeF, meanPSD);

% Print the outcome in an xls file
xlswrite ('EEG_psd_SI1.5_20221027_Treat.xlsx', norm_powers, 1, 'B2');
xlswrite ('EEG_psd_SI1.5_20221027_Treat.xlsx', labels, 1, 'B1');
xlswrite ('EEG_psd_SI1.5_20221027_Treat.xlsx', HLF_Ratio', 1, 'I2');
xlswrite ('EEG_psd_SI1.5_20221027_Treat.xlsx', table_index, 2, 'B2');
xlswrite ('EEG_psd_SI1.5_20221027_Treat.xlsx', WholeF, 3, 'A2');
xlswrite ('EEG_psd_SI1.5_20221027_Treat.xlsx', PSDs, 3, 'B2');
xlswrite ('EEG_psd_SI1.5_20221027_Treat.xlsx', meanPSD, 3, 'K2');
xlswrite ('EEG_psd_SI1.5_20221027_Treat.xlsx', SD_mean, 3, 'L2');
xlswrite ('EEG_psd_SI1.5_20221027_Treat.xlsx', SEM, 3, 'M2');