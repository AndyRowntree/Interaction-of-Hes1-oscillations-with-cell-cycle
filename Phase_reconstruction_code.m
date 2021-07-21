%% Phase reconstruction

% This code finds the phase readout (between 0 - 2pi) of each cell at each time point throughout its whole trace
% between mitosis and mitosis, from this we take the phase
% readout at the begninning mitosis. This is compared with the cell's
% cell cycle length, cluster classification and/or subtype throughout Fig 8.

% Assign variable names "Hes1_raw_traces_two_generations" to the imported traces 
% These are 1 x n (n is number of cells, i.e. 152 in our work for our first generation) arrays where each item is a time series trace vector. 
% These can be any length but for our work we used pseudo-time traces to ensure our smoothing kept the same moving window length.

Hes1_psuedo_time_traces = Insert data here;

% Preallocation of arrays and vectors

Gaussian_smooth_cells = cell(1,length(Hes1_psuedo_time_traces));
Phase_of_cell{cell_index} = cell(1,length(Hes1_psuedo_time_traces));
Hes1_phase_readout_at_beginning_of_cell_cycle = zeros(1,length(Hes1_psuedo_time_traces));

% loop through each cell

window = 20;

for cell_index = 1:length(Hes1_raw_traces_two_generations) 
 
    % Apply a Gaussian smoothing filter (slightly different smoothing
    % method to dip detection as we found that this particular one yields 
    % a more accurate readout at the beginning and end of traces). Here
    % we use a moving window of 20 time points. We recall that the pseudo-time traces are 144 time points in length.
   
    pre_smooth =  smoothdata(Hes1_psuedo_time_traces,'gaussian',window);
    Gaussian_smooth_cells{cell_index} =  smoothdata(pre_smooth,'gaussian',window);

    % We apply the smoothing filter a second time with identical parameters
    % to further smooth out any kinks missed by the first smoothing filter
    % and thus provide a better phase readout later on.

    % Apply the Hilbert transform on our smoothened traces.
    % This returns a complex helical sequence, sometimes called the analytic signal, from a real data sequence.
    
    Phase_of_cell{cell_index} = mod(angle(hilbert(Gaussian_smooth_cells{cell_index})),2*pi);
         
    % Using 'angle', we can obtain the 'phase readout' of our data as a
    % real-valued angle between -pi to pi radians from the complex-valued
    % Hilbert transform output. We chose to have our phase readout positive
    % to ease understanding for non-mathematicians thus we apply modulo 2pi to the output 
    % value and obtain a value between 0 to 2pi where 0 and 2pi correspoind
    % to a peak in the wave and pi corresponds to a trough in the wave.
    
    % We take the first value of the phase readout in each trace to give us
    % a final readout of the phase of our Hes1 wave at the beginning of
    % the cell cycle.
    
    Hes1_phase_readout_at_beginning_of_cell_cycle(cell_index) = Phase_of_cell{cell_index}(1);
end
