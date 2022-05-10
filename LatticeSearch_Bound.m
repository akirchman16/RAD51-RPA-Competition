%Function to search the ssDNA lattice to find the number of bound proteins
%of each type and where they're bound. This will (MAYBE - Check!) remove
%the need for BoundAtSpot arrays.

% Location Arrays: Lists of all locations where each type of protein is
% bound. The protein type is in the variable name.
% Counts: Number of bound (or hinged open) proteins of each type.
%     1 - RPA-A Bound
%     2 - RPA-A Open
%     3 - RPA-D Bound
%     4 - RPA-D Open
%     5 - RAD51 Monomer
%     6 - RAD51 Dimer

function [x] = LatticeSearch_Bound(DNA,n_RAD51,n_A,n_D)
    n_RPA = n_A+n_D;
    
    %Initialize vectors for position arrays.
    RPA_A_Bound_Pos = [];       RPA_A_Open_Pos = [];
    RPA_D_Bound_Pos = [];       RPA_D_Open_Pos = [];
    RAD51_Mon_Pos = [];         RAD51_Dim_Pos = [];
    
    % 1: RPA-A (Bound)
    if ~isempty(find(DNA(2,:) == 1, 1))
        % RPA-A clusters that are bound to ssDNA and their size
        RPA_A_Bound_Cluster = logical(DNA(2,:) == 1);
        bb_RPA_A_Bound = regionprops(RPA_A_Bound_Cluster,'BoundingBox'); bb_RPA_A_Bound = cat(1,bb_RPA_A_Bound.BoundingBox);
        RPA_A_Bound_Cluster = [bb_RPA_A_Bound(:,1)+0.5, bb_RPA_A_Bound(:,3)];  %[Left edge of cluster, Cluster Size]
        % Determine the number of RPA-A proteins bound together in each
        % cluster
        RPA_A_Bound_NumProteins = (RPA_A_Bound_Cluster(:,2)/n_A).';
        % Record positions of each bound RPA-A protein
        for i = 1:numel(RPA_A_Bound_NumProteins)
            RPA_A_Bound_Pos = sort([RPA_A_Bound_Pos,RPA_A_Bound_Cluster(i,1)+(n_A*(0:1:RPA_A_Bound_NumProteins(i)-1))]);
        end
    end
    % 2: RPA-A (Open)
    if ~isempty(find(DNA(1,:) == 1,1))
        % RPA-A cluster which are hinged open (in first row of ssDNA
        % lattice)
        RPA_A_Open_Cluster = logical(DNA(1,:) == 1);
        bb_RPA_A_Open = regionprops(RPA_A_Open_Cluster,'BoundingBox'); bb_RPA_A_Open = cat(1,bb_RPA_A_Open.BoundingBox);
        RPA_A_Open_Cluster = [bb_RPA_A_Open(:,1)+0.5, bb_RPA_A_Open(:,3)];  %[Left edge of cluster, Cluster Size]
        % Determine the number of RPA-A proteins that are hinged open
        % together in each cluster
        RPA_A_Open_NumProteins = (RPA_A_Open_Cluster(:,2)/n_A).';
        % Record the positions of each hinged open RPA-A protein
        for i = 1:numel(RPA_A_Open_NumProteins)
            RPA_A_Open_Pos = sort([RPA_A_Open_Pos,RPA_A_Open_Cluster(i,1)+(n_A*(0:1:RPA_A_Open_NumProteins(i)-1))]);
        end
    end
    % 3: RPA-D (Bound)
    if ~isempty(find(DNA(2,:) == 3, 1))
        % RPA-D clusters which are bound to the ssDNA lattice
        RPA_D_Bound_Cluster = logical(DNA(2,:) == 3);
        bb_RPA_D_Bound = regionprops(RPA_D_Bound_Cluster,'BoundingBox'); bb_RPA_D_Bound = cat(1,bb_RPA_D_Bound.BoundingBox);
        RPA_D_Bound_Cluster = [bb_RPA_D_Bound(:,1)+0.5, bb_RPA_D_Bound(:,3)];   %[Left edge of cluster, Cluster Size]
        % Determine the number of RPA-D proteins that are bound together in
        % each cluster along the ssDNA lattice
        RPA_D_Bound_NumProteins = (RPA_D_Bound_Cluster(:,2)/n_D).';
        % Record positions of each bound RPA-D protein
        for i = 1:numel(RPA_D_Bound_NumProteins)
            RPA_D_Bound_Pos = sort([RPA_D_Bound_Pos,RPA_D_Bound_Cluster(i,1)+(n_D*(0:1:RPA_D_Bound_NumProteins(i)-1))]);
        end
    end
    % 4: RPA-D (Open)
    if ~isempty(find(DNA(1,:) == 3,1))
        % RPA-D clusters which are hinged open together
        RPA_D_Open_Cluster = logical(DNA(1,:) == 3);
        bb_RPA_D_Open = regionprops(RPA_D_Open_Cluster,'BoundingBox'); bb_RPA_D_Open = cat(1,bb_RPA_D_Open.BoundingBox);
        RPA_D_Open_Cluster = [bb_RPA_D_Open(:,1)+0.5, bb_RPA_D_Open(:,3)];  %[Left edge of cluster, Cluster Size]
        % Determine the number of RPA-D proteins bound together in each
        % cluster along the ssDNA lattice
        RPA_D_Open_NumProteins = (RPA_D_Open_Cluster(:,2)/n_D).';
        % Record the positions of each hinged open RPA-D protein
        for i = 1:numel(RPA_D_Open_NumProteins)
            RPA_D_Open_Pos = sort([RPA_D_Open_Pos,RPA_D_Open_Cluster(i,1)+(n_D*(0:1:RPA_D_Open_NumProteins(i)-1))]);
        end
    end
    % 5: RAD51 (this will be used to split up into monomers vs. dimers)
    if ~isempty(find(DNA(2,:) == 51,1))
        % RAD51 clusters whcih are bound together on ssDNA lattice
        RAD51_Cluster = logical(DNA(2,:) == 51);
        bb_RAD51 = regionprops(RAD51_Cluster,'BoundingBox'); bb_RAD51 = cat(1,bb_RAD51.BoundingBox);
        RAD51_Cluster = [bb_RAD51(:,1)+0.5, bb_RAD51(:,3)]; %[Left edge of cluster, Cluster Size]
        % Determine the number of RAD51 monomers and dimers which are bound
        % together in each cluster along the ssDNA lattice
        RAD51_Mon_NumProteins = (RAD51_Cluster(:,2)/n_RAD51).';
        RAD51_Dim_NumProteins = (RAD51_Cluster(:,2)/(2*n_RAD51)).';
        RAD51_Dim_NumProteins(RAD51_Dim_NumProteins > 1) = ceil(RAD51_Dim_NumProteins(RAD51_Dim_NumProteins > 1));  %takes care of 3 monomers bound together being 2 possible dimers
        RAD51_Dim_NumProteins(RAD51_Dim_NumProteins < 1) = floor(RAD51_Dim_NumProteins(RAD51_Dim_NumProteins < 1)); %rounds number of proteins to whole numbers when less than 1 (0.5 -> Monomer)
        % Record the positions of each RAD51 monomer
        for i = 1:numel(RAD51_Mon_NumProteins)
            RAD51_Mon_Pos = sort([RAD51_Mon_Pos,RAD51_Cluster(i,1)+(n_RAD51*(0:1:RAD51_Mon_NumProteins(i)-1))]);
        end
        %Record the positions of each RAD51 dimer
        RAD51_Dim_NumProteins = RAD51_Dim_NumProteins(RAD51_Dim_NumProteins ~= 0);
        RAD51_Dim_Cluster = RAD51_Cluster(RAD51_Cluster(:,2) > n_RAD51);
        for j = 1:numel(RAD51_Dim_NumProteins)
            RAD51_Dim_Pos = sort([RAD51_Dim_Pos,RAD51_Dim_Cluster(j,1)+(n_RAD51*(0:1:RAD51_Dim_NumProteins(j)-1))]);
        end
    end
    
    %Record the number of bound proteins for each type
    RPA_A_Bound_Count = numel(RPA_A_Bound_Pos);     RPA_A_Open_Count = numel(RPA_A_Open_Pos);
    RPA_D_Bound_Count = numel(RPA_D_Bound_Pos);     RPA_D_Open_Count = numel(RPA_D_Open_Pos);
    RAD51_Mon_Count = numel(RAD51_Mon_Pos);         RAD51_Dim_Count = numel(RAD51_Dim_Pos);
    
    %Check "Open" proteins to see if they can hinge close (micro-binding).
    %Record where they are and count them.
    if ~isempty(RPA_A_Open_Pos)
        RPA_A_Open_Bindable = RPA_A_Open_Pos(all(DNA(2,RPA_A_Open_Pos:RPA_A_Open_Pos+(n_A-1)) == 0,'all'));
    else
        RPA_A_Open_Bindable = [];
    end
    if ~isempty(RPA_D_Open_Pos)
        RPA_D_Open_Bindable = RPA_D_Open_Pos(all(DNA(2,RPA_D_Open_Pos:RPA_D_Open_Pos+(n_D-1)) == 0,'all'));
    else
        RPA_D_Open_Bindable = [];
    end
    RPA_A_Bindable_Count = numel(RPA_A_Open_Bindable);  RPA_D_Bindable_Count = numel(RPA_D_Open_Bindable);
    
    %Check "Bound" proteins to see if they can hinge open
    %(micro-unbinding). Record where these are at and count them.
    if ~isempty(RPA_A_Bound_Pos)
        RPA_A_Bound_Hingeable = RPA_A_Bound_Pos(all(DNA(1,RPA_A_Bound_Pos:RPA_A_Bound_Pos+(n_A-1)) == 0,'all'));
    else
        RPA_A_Bound_Hingeable = [];
    end
    if ~isempty(RPA_D_Bound_Pos)
        RPA_D_Bound_Hingeable = RPA_D_Bound_Pos(all(DNA(1,RPA_D_Bound_Pos:RPA_D_Bound_Pos+(n_D-1)) == 0,'all'));
    else
        RPA_D_Bound_Hingeable = [];
    end
    RPA_A_Hingeable_Count = numel(RPA_A_Bound_Hingeable);   RPA_D_Hingeable_Count = numel(RPA_D_Bound_Hingeable);
    
    %Organize everything into neat arrays to record final cada
end