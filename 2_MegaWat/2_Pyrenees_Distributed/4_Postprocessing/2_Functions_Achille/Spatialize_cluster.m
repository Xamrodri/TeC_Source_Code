function [VAR_map, Date] = Spatialize_cluster(Cluster_map,table_cluster,newTimeStep,var_to_spatialize)
%SPATIALIZE_CLUSTER Summary of this function goes here
%   Detailed explanation goes here

no_clusters = size(table_cluster,1); % number of clusters

col_var = strcmp(table_cluster{1,1}.Properties.VariableNames,var_to_spatialize); % Columns of the variable to spatialize

% Compute number of elements after retime

TT_newtime = retime(table_cluster{1,1}(:,col_var),newTimeStep,'mean');
Date = TT_newtime.Date;
num_el = numel(Date);


VAR_map = repmat(Cluster_map.*NaN,[1 1 num_el]);

table_cluster_c = cell(no_clusters,1); % allocate memory

% First, retime all the clusters timetable to the map frequency resolution
for c = 1:no_clusters
    if strcmp(var_to_spatialize,'Pr_poi') || strcmp(var_to_spatialize,'Pr_sno_poi')
       table_cluster_c{c} = retime(table_cluster{c,1}(:,col_var),newTimeStep,@nansum);
    else
       table_cluster_c{c} = retime(table_cluster{c,1}(:,col_var),newTimeStep,@nanmean);
    end 
end 

parfor dt = 1:numel(Date)
  VAR_i = Cluster_map.*NaN;
  for c = 1:no_clusters
    %Retime according to the vector given
     cur = Cluster_map == table_cluster{c,4};
     table_cluster_i = table_cluster_c{c};
     VAR_i(cur) = table_cluster_i.(var_to_spatialize)(dt);
   end
   VAR_map(:,:,dt) = flipud(VAR_i);
end

end

