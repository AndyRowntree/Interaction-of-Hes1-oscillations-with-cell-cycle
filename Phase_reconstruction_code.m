%% Phase reconstruction

% This code finds the phase readout (between 0 - 2pi) of each cell at each time point throughout its whole trace
% (along with the trace of it's daughter cell), from this we take the phase
% readout at the begninning mitosis. This is later compared with the cell's
% cell cycle length, cluster classification and/or subtype when necessary.

% Assign variable names "Hes1_raw_traces_two_generations" to the imported traces 
% These are 1 x n (n is number of cells, i.e. 145 in our work) arrays where each item is a time series trace vector.

%Hes1_raw_traces_two_generations = Insert data here;

% Preallocation of arrays and vectors

Hes1_smoothened_traces_two_generations = cell(1,lengtH(Hes1_raw_traces_two_generations));
Hes1_twice_smoothened_traces_two_generations = cell(1,lengtH(Hes1_raw_traces_two_generations));
Linear_fit = cell(1,lengtH(Hes1_raw_traces_two_generations));
Hes1_twice_smoothened_and_detrended_traces_two_generations = cell(1,lengtH(Hes1_raw_traces_two_generations));
Hilbert_transform_of__smoothened_detrended_Hes1_trace = cell(1,lengtH(Hes1_raw_traces_two_generations));
Hes1_phase_readout_of_detrended_smoothened_curve = cell(1,lengtH(Hes1_raw_traces_two_generations));
Hes1_phase_readout_at_beginning_of_cell_cycle = zeros(1,lengtH(Hes1_raw_traces_two_generations));

% loop through each cell

for cell_index = 1:length(Hes1_raw_traces_two_generations) 
    
    % Apply a Savitzgy-Golay smoothing filter (slightly different smoothing
    % parameters to dip detection due to differences in trace length). Here
    % we use a moving frame length of half of the trace length and a
    % polynomial order of 4.
    
    % if loop ensures that the moving frame length is odd, as required for the smoothing filter to work 
    if mod(round(length(Hes1_raw_traces_two_generations)/2),2)==1
        smoothing_frame_length = round(length(Hes1_raw_traces_two_generations)/2);
    else
        smoothing_frame_length = round(length(Hes1_raw_traces_two_generations)/2)-1;
    end
 
    Hes1_smoothened_traces_two_generations{cell_index} = sgolayfilt(Hes1_raw_traces_two_generations{cell_index},4,smoothing_frame_length);    
    
    % We apply the smoothing filter a second time with identical parameters
    % to further smooth out any kinks missed by the first smoothing filter
    % and thus provide a better phase readout later on.
    
    Hes1_twice_smoothened_traces_two_generations{cell_index} = sgolayfilt(Hes1_smoothened_traces_two_generations{cell_index},4,smoothing_frame_length);
    
    % We next remove a linear fit to the smoothened trace, this can be done
    % a number of different ways: we choose to again use a Savitzgy-Golay filter
    % with a frame length of the entire trace length and with polynomial
    % order of 1.

    % Again, we ensure that the frame length is odd.
    if mod(length(Hes1_twice_smoothened_traces_two_generations{cell_index}),2)==1
        linear_frame_length = length(Hes1_twice_smoothened_traces_two_generations{cell_index});
    else
        linear_frame_length = length(Hes1_twice_smoothened_traces_two_generations{cell_index})-1;
    end

    Linear_fit{cell_index}=sgolayfilt(Hes1_twice_smoothened_traces_two_generations{cell_index},1,linear_frame_length);
    
    % Obtain smoothened and detrended trace by subtracting our linear trend
    % line from our smoothened trace.

    Hes1_twice_smoothened_and_detrended_traces_two_generations{cell_index}=Hes1_twice_smoothened_traces_two_generations{cell_index}-Linear_fit{cell_index};
    
    % Apply the Hilbert transform on our smoothened, deterended traces.
    % This returns a complex helical sequence, sometimes called the analytic signal, from a real data sequence.
    Hilbert_transform_of__smoothened_detrended_Hes1_trace{cell_index}=hilbert(Hes1_twice_smoothened_and_detrended_traces_two_generations{cell_index});
    
    % Using 'angle', we can obtain the 'phase readout' of our data as a
    % real-valued angle between -pi to pi radians from the complex-valued
    % Hilbert transform output. We chose to have our pahse readout positive
    % to ease understanding for non-mathematicians thus we add pi to each
    % value and obtain a value between 0 to 2pi where 0 and 2pi correspoind
    % to a trough in the wave and pi corresponds to a peak in the wave.
    
    Hes1_phase_readout_of_detrended_smoothened_curve{cell_index} = angle(Hilbert_transform_of__smoothened_detrended_Hes1_trace{cell_index})+pi;
    
    % We take the first value of the phase readout in each trace to give us
    % a final readout of the phase of our Hes1 wave at the beginning of
    % the cell cycle.
    
    Hes1_phase_readout_at_beginning_of_cell_cycle(cell_index) = Hes1_phase_readout_of_detrended_smoothened_curve{cell_index}(1);
    
    
end
