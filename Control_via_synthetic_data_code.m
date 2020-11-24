% Clustering and interpolation control via synthetic data for "Phase register of Hes1 oscillations
% with mitoses underlies cell-cycle heterogeneity in ER+ breast cancer
% cells." (2020) by Nitin Sabherwal, Andrew Rowntree, Jochen Kursawe and Nancy Papalopulu*

% This code is for a synthetic control which ensures that the clustering
% performed within the paper does not produce clusters based solely on stretched data.
% I.e. it is the underlying dynamic shapes which are being picked up by the
% clustering algortihm and not merely data which looks similart just
% because it's been stretched or squashed.

% Preallocate an array for our 5 dynamic shapes

array_of_control_dynamic_shapes={};

% Set three different length scales (akin to the the lengths of our data,
% this represents heterogeneity in cell cycle lengths).

short_control_length = 80;
med_control_length = 140;
long_control_length = 200;
three_lengths = [short_control_length,med_control_length,long_control_length];

% Three of our dynamic shapes will be based on actual traces from our data,
% we choose these by randomly selecting three from our 152 generation 1
% cells.

three_example_data_traces_indices = randi(152,1,3);

% Loop through each of our three control trace lengths our three data base
% traces so we develop 9 control traces of varying length and shape. This
% is done by applying a polynomial fit to the trace and interpolating to
% each length.

% We also add anti-phase sin waves, each at the three different lengths (6
% total) to our control data. The idea here is that these two shapes are
% completely oppposite in dynamics and should never cluster together.

% This gives us 15 control traces

% We can import our data here
% Hes1_raw_traces = Insert data here;

for len = 1:length(three_lengths)
    for data_trace = 1 : length(three_example_data_traces_indices)        
        f=polyfit((1:numel(Hes1_raw_traces{three_example_data_traces_indices(data_trace)}))',Hes1_raw_traces{three_example_data_traces_indices(data_trace)} ,11);
        m=polyval(f,(1:numel(Hes1_raw_traces{three_example_data_traces_indices(data_trace)}))'); 
        vector_of_steps = linspace(1,numel(m),three_lengths(len));
        array_of_control_dynamic_shapes{end+1} =  interp1(1:numel(m),m,vector_of_steps);
    end    
    t=linspace(0,2*pi,three_lengths(len));
    array_of_control_dynamic_shapes{end+1} = sin(t);
    array_of_control_dynamic_shapes{end+1} = -sin(t);
end

% We now duplicate each of our 15 control traces 10 times each and add
% noise. This noise is taken from a standard normal distribution but with
% standard deviation of a half instead of 1. This is similar to the noise
% found in our actual data on average. This will then yield 150 control
% traces which nearly matches our experimental data set (n=152).

dup_array_of_control_normalised_noisy_dynamic_shapes={};

for trace_index = 1:length(array_of_control_dynamic_shapes)
    for duplicate = 1:10
        noise_vector = (1/2).*randn(1,length(array_of_control_dynamic_shapes{trace_index}));
        dup_array_of_control_normalised_noisy_dynamic_shapes{end+1} = zscore(array_of_control_dynamic_shapes{trace_index})+noise_vector;
    end
end

% We now stretch the data so that they are all the same length and can be
% clustered.

time_values_for_stretching = length(dup_array_of_control_normalised_noisy_dynamic_shapes)-1;
% Preallocate matrix for clustering
Matrix_of_synth_data = zeros(time_values_for_stretching,length(dup_array_of_control_normalised_noisy_dynamic_shapes));

% Stretching is done using interp1
for synth_cell_index = 1:length(dup_array_of_control_normalised_noisy_dynamic_shapes)
    vector_of_steps = linspace(1,numel(dup_array_of_control_normalised_noisy_dynamic_shapes{synth_cell_index}),time_values_for_stretching);
    Matrix_of_synth_data(:,synth_cell_index) = interp1(1:numel(dup_array_of_control_normalised_noisy_dynamic_shapes{synth_cell_index}),...
        dup_array_of_control_normalised_noisy_dynamic_shapes{synth_cell_index},vector_of_steps);
end

% We now randomly mix te data so that lcustering does not start of with
% ordered data

Matrix_of_synth_data=Matrix_of_synth_data(:,randperm(length(dup_array_of_control_normalised_noisy_dynamic_shapes)));

% This matrix of randomly ordered synthetic traces can be inputted
% into the following analysis.

% We next performed a Silhouette analysis for many different cluster numbers
% and all possible Gaussian covairiance matrix types to determine the most
% suitable GMM clustering for the synthetic data we have created here. As
% the results in the paper show, this analysis provided the answer as k=5
% with a shared and diagonal covariance matrix as the most suitable mixture
% model for this synthetic data. Clustering with this showed that the 5
% original dynamical shapes do indeed cluster separately, justifying our
% claim that clustering is not an artifact of stretching

% The code used for this analysis is identical to that used for the Silhouette and
% clustering anlaysis of our experimental data and thus is described in
% more detail in the "Silhouette_analysis_and_clustering_code.m" file.




