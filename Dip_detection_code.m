% Dip detection method for Phase register of Hes1 oscillations 
% with mitoses underlies cell-cycle heterogeneity in ER+ breast cancer
% cells. (2020) by Nitin Sabherwal, Andrew Rowntree, Jochen Kursawe and Nancy Papalopulu*

% This code detects a dip in expression of Hes1 via smoothing with a Savitzgy
% Golay filter and looking for the minimum turning point (after a certain time) of such
% smoothening. 

% Assign variable names "Hes1_raw_traces" to the imported traces 
% These are 1 x n (n is number of cells, i.e. 152 in our work) arrays where each item is a time series trace vector.

load array_of_cells
Hes1_raw_traces = array_of_cells{2};

% Preallocate vectors and arrays

Hes1_smoothened_traces = cell(1,length(Hes1_raw_traces));
Turning_points = cell(1,length(Hes1_raw_traces));
Detected_dip = zeros(1,length(Hes1_raw_traces));
Time_from_start_of_cell_cycle_until_dip = zeros(1,length(Hes1_raw_traces));
Time_from_dip_until_end_of_cell_cyle = zeros(1,length(Hes1_raw_traces));

% Loop through all cells

for cell_index = 1:length(Hes1_raw_traces)
    
    % Set a moving frame length for the Savitzgy-Golay filter, we choose
    % this to be at a quarter of the length of each trace. The fram length
    % must be odd in order to work, hence the if loop. We also set the
    % order of polynmomial as 3, this means a local cubic fit is performed
    % for each time point as the frame moves across the trace.
    
    if mod(round((1/4)*(length(Hes1_raw_traces{cell_index}))),2) == 1
        % if a quarter of the trace length is odd, then use this value
        frame_length_for_filter = round((1/4)*(length(Hes1_raw_traces{cell_index})));
    else
        % if a quarter of the trace length is even, then use one less than this value
        frame_length_for_filter = round((1/4)*(length(Hes1_raw_traces{cell_index})))-1; %MUST BE ODD!
    end
    
    % Apply the Savitzgy-Golay filter
    
    Hes1_smoothened_traces{cell_index} = sgolayfilt(Hes1_raw_traces{cell_index},3,frame_length_for_filter);
    
    % Chop off the ends of the smoothened trace, to abandon smoothing artefacts 
    
    Truncated_smoothened_trace = Hes1_smoothened_traces{cell_index}(3:end-2);
    
    Turning_points{cell_index} = [];
    
    % Find the turning points of the smooth trace. This is done by
    % multiplying the difference of three consecutive time points.
    % For instance, if we have an upward trend and a downward trend either
    % side of a time point then this product will be negative, i.e. a
    % turning point. We record every time this is the case.
    
    for time_index = 2:length(Truncated_smoothened_trace)-2
        if (Truncated_smoothened_trace(time_index) - Truncated_smoothened_trace(time_index+1))*...
                (Truncated_smoothened_trace(time_index+1) - Truncated_smoothened_trace(time_index+2)) < 0
            Turning_points{cell_index} = [Turning_points{cell_index},(time_index+3)];
        end
    end
    
    turning_points_not_at_start_of_cell_cycle = [];

    % Some turning points are captured early in the cell cycle, we know that
    % cells sometimes have an early dip in expression, yet this not the
    % only dip in the whole trace and not the one we are looking for. To
    % account for this, we disgard turning points in the first 30% of the
    % cell cycle.
    
    for t_point = 1:length(Turning_points{cell_index})
        if  Turning_points{cell_index}(t_point) > (length(Hes1_smoothened_traces{cell_index}))*(3/10)
            turning_points_not_at_start_of_cell_cycle = [turning_points_not_at_start_of_cell_cycle, Turning_points{cell_index}(t_point)];
        end
    end
    
    % Our detected dip is now considered to be the lowest (minimum) of such remaining
    % turning points.
    
    Detected_dip(cell_index) = find(Hes1_smoothened_traces{cell_index} == min(Hes1_smoothened_traces{cell_index}(turning_points_not_at_start_of_cell_cycle))); %    
    Time_from_start_of_cell_cycle_until_dip(cell_index) = Detected_dip(cell_index)/4;
    Time_from_dip_until_end_of_cell_cyle(cell_index) = (length(Hes1_smoothened_traces{cell_index})-Detected_dip(cell_index))/4;
    
end
