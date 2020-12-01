% Silhouette analysis and clustering code for "Phase register of Hes1 oscillations 
% with mitoses underlies cell-cycle heterogeneity in ER+ breast cancer
% cells." (2020) by Nitin Sabherwal, Andrew Rowntree, Jochen Kursawe and Nancy Papalopulu*

% This code performs a Silhouette analysis which can provide an answer to
% which number of clusters is most suitable to use for a given set of data.
% We also use this method to choose the most suitable restrictions to the
% covariance matrices of the Gaussian components which will make up our
% mixture model in the clustering.

% Assign variable names "Hes1_interpolated_traces" to the imported traces 
% These are 1 x n (n is number of cells, i.e. 152 in our work) arrays where
% each item is a peudo-time series trace vector.

%       Hes1_interpolated_traces = Insert data here;

% Covariance matrix options to be looped through 
% Gaussian components identical (shared) or not (unshared).
Shared_Covariance = {true,false};

% Gaussian components only orthoganally stretchable (diagonal) or not
% (full).
Sigma = {'full','diagonal'}; 

% Setting number of EM iterations
options = statset('MaxIter', 75536);   

% Max number of mixtures/clusters to check
kmax = 8; 

% Preallocate a three dimensional array (one for each component number we check, 
% ones for shared/unshared and one for full/diagonal covariance matrices), 
% this is where each GMM shall go, 
This_GMM = cell(kmax-1,2,2);

% Create a vector of cell index values
Vector_of_indices = 1:length(Hes1_interpolated_traces);

% Set a number of repetitions to run (we run many repetitions because
% the GMM process initially chooses ranomdised means and thus local 
% opitma may be foun with the mixture models meaning slightly different
% clustering results may come out each time. We plan to average out
%all repetitions to get a more reliable result).
no_of_reps = 10;

% Preallocate a vector with rows being each repetition and column being
% component numbers

Silhouette_score = zeros(no_of_reps,kmax);

% The above vector will be averages for each covarianve matrix type,
% this is preallocated into an array here.
Silhouette_score_per_cov_type = cell(2,2);

% We now loop through our choice of four covariance matrix types
for SCov = 1:2
    for Sigloop = 1:2
        
        % We loop through our set number of repetitions
        for repeat = 1:no_of_reps
            
            
            % For each repetition we randomise our order of beginging
            % indices. This makes for fairer results
            random_order = Vector_of_indices(randperm(length(Vector_of_indices)));
            
            % We loop through each component number (K=1 is omitted since
            % one cluster is trivial to us).
            for k_choice = 2:kmax               
                
                % Fits a k_choice-component GMM to our randomly ordered
                % traces (NB, the transpose is required otherwise time
                % will be clustered and not traces).
                This_GMM = fitgmdist(Hes1_interpolated_traces(random_order)',...
                    k_choice,'RegularizationValue',0.01,'CovarianceType',Sigma{Sigloop},...
                    'SharedCovariance',Shared_Covariance{SCov},'Options',options);
                
                % Assigns each cell a number based on which Gaussian
                % component the above fit decides
                GMM_clusterX = cluster(This_GMM,Hes1_interpolated_traces');
                
                % Preallocates an array where the index of each cluster 
                % of cells shall be placed.
                Cluster_index_arrays = cell(1,k_choice);
                
                % Place all cell indices into an element of the above array
                % depending on it's newly assigned cluster
                for cell_index = 1:size(Gen_1_pseudo_time,2)
                    for this_k= 1:k_choice
                        if GMM_clusterX(cell_index) == this_k
                            Cluster_index_arrays{this_k}(end+1) = cell_index;
                        end
                    end
                end
                
                % Next we perform the Silhouette method to develop a 
                % "Silhouette score" for each covarince matirx type
                % and component number 
                
                % Preallocate and fill in an array which will contain a
                % mean traces of each cluster
                
                mean_of_clusters = cell(1,k_choice);
                for k_now = 1:k_choice                    
                    mean_of_clusters{k_now} = mean(Hes1_interpolated_traces(:,Cluster_index_arrays{k_now})');
                end
                
                % Set up a vector for Silhoette coefficients, found later
                Silhouette_coefficient=[];
                
                % Preallocate a vector which will be filled with
                % information of which cluster is closest in similarity to
                % our cluster of interest as we loop through
                closest_cluster_to_this=zeros(1,k_choice);
                
                % Loop through all clusters we are checking
                for k_now = 1:k_choice
                    % Set up a vector which will determine the distance
                    % (euclidean) between this cluster and each the others
                    Dist_of_this_clust_from_other=[];
                    for k_others = 1:k_choice
                        Dist_of_this_clust_from_other(end+1) = sqrt(sum((mean_of_clusters{k_now}-...
                            mean_of_clusters{k_others}).^ 2));
                    end
                    
                    % Here we find which cluster is closest to our cluster
                    % of interest
                    closest_cluster_to_this(k_now) = find(Dist_of_this_clust_from_other==...
                        min(Dist_of_this_clust_from_other(Dist_of_this_clust_from_other>0)));
                    
                    % Set up a vector which will be filled with the
                    % Silhouette coefficients for each cluster
                    Silhouette_coeffficient_per_cluster=zeros(1,length(Cluster_index_arrays{k_now}));
                    
                    % For each clucter, we loop through each cell in that
                    % cluster to obtain the inputs for the Silhouette
                    % coefficient formula
                    for cell_index = 1:length(Cluster_index_arrays{k_now})
                        % Assign a name to the cluster we are working with
                        cluster_without_this_one = Cluster_index_arrays{k_now};
                        
                        % Remove our cell of interest from the cluster 
                        cluster_without_this_one(cell_index) = [];
                        
                        % Set up a vector which will be filled with the
                        % distance (Euclidean) between our cell of interest
                        % and all other cells in that cluster as we loop
                        % through
                        Distance_from_same_cluster_cell = [];
                        for other_cell = 1:length(cluster_without_this_one)                            
                            Distance_from_same_cluster_cell(end+1) = sqrt(sum((Hes1_interpolated_traces(:,Cluster_index_arrays{k_now}(cell_index)) - ...
                                Hes1_interpolated_traces(:,cluster_without_this_one(other_cell))).^ 2));                                                        
                        end
                        
                        % Calculate the mean of distances between our cell of
                        % interest and each other cell in the same cluster
                        avg_for_same_cluster(cell_index) = mean(Distance_from_same_cluster_cell);
                        
                        % We now, similarly, calculate the mean of
                        % distances  between our cell of interesst and all
                        % cells in the nearest cluster (nearest cluster was
                        % determined earlier)
                        Distance_from_other_cluster_cell=[];                       
                        for other_cell = 1:length(Cluster_index_arrays{closest_cluster_to_this(k_now)})                            
                            Distance_from_other_cluster_cell(end+1) = sqrt(sum((Hes1_interpolated_traces(:,Cluster_index_arrays{k_now}(cell_index)) - ...
                                Hes1_interpolated_traces(:,Cluster_index_arrays{closest_cluster_to_this(k_now)}(other_cell))).^ 2));
                        end
                        avg_for_other_cluster(cell_index) = mean(Distance_from_other_cluster_cell);
                        
                        
                        % We now have values for the average distance of
                        % each cell from all other cells in its cluster and
                        % all cells its nearest cluster. This is enough to
                        % find the Silhouette coefficient for each cluster
                        % number choice.

                        Silhouette_coeffficient_per_cluster(cell_index) = (avg_for_other_cluster(cell_index)-avg_for_same_cluster(cell_index))/...
                            max(avg_for_other_cluster(cell_index),avg_for_same_cluster(cell_index));
                        
                    end

                    % Add each silhouette coefficient to the following
                    % vector. 
                    Silhouette_coefficient=[Silhouette_coefficient,Silhouette_coeffficient_per_cluster];    
                end
                
                % We add the mean silhouette score for each cell across
                % clusters and a row vector of Silhouettte scores for each
                % cluster choice from column 1 being 2 clusters to the end
                % column being the number of clusters we check. Each row is
                % now a different repeat
                Silhouette_score(repeat,k_choice) = nanmean(Silhouette_coefficient);
            end  
        end          
        Silhouette_score_per_cov_type{Sigloop,SCov} = Silhouette_score;   
    end
end

% We now have many repeats fot Sihouette scores for each given cluster
% choice at each of four different covariance matirix types. These can be
% averaged out and the component number choice and covariance type for which the
% Silhouette score is greatest will indicate our most suitable clustering
% options for GMM clustering of our data.



