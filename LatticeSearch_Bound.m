%Function to search the ssDNA lattice to find the number of bound proteins
%of each type and where they're bound. This will (MAYBE - Check!) remove
%the need for BoundAtSpot arrays.

% I THINK I HAVE TO DO EACH PROTEIN IN IT'S OWN INDIVIDUAL IF STATEMENT IN
% CASE THE LATTICE IS EMPTY OR FREE OF THAT PROTEIN!!!!!!

function [x] = LatticeSearch_Bound(DNA,n_RAD51,n_A,n_D)
    n_RPA = n_A+n_D;
    %First, search again for clusters of each protein type. Record how
    %large the strands of bound proteins are and where their edges are.
    % 1: RPA-A (Bound)
    if ~isempty(find(DNA(2,:) == 1, 1))
        RPA_A_Bound_Cluster = logical(DNA(2,:) == 1);
        bb_RPA_A_Bound = regionprops(RPA_A_Bound_Cluster,'BoundingBox'); bb_RPA_A_Bound = cat(1,bb_RPA_A_Bound.BoundingBox);
    end
    % 2: RPA-A (Open)
    RPA_A_Open_Cluster = logical(DNA(1,:) == 1);
    bb_RPA_A_Open = regionprops(RPA_A_Open_Cluster,'BoundingBox'); bb_RPA_A_Open = cat(1,bb_RPA_A_Open.BoundingBox);
    % 3: RPA-D (Bound)
    RPA_D_Bound_Cluster = logical(DNA(2,:) == 3);
    bb_RPA_D_Bound = regionprops(RPA_D_Bound_Cluster,'BoundingBox'); bb_RPA_D_Bound = cat(1,bb_RPA_D_Bound.BoundingBox);
    % 4: RPA-D (Open)
    RPA_D_Open_Cluster = logical(DNA(1,:) == 3);
    bb_RPA_D_Open = regionprops(RPA_D_Open_Cluster,'BoundingBox'); bb_RPA_D_Open = cat(1,bb_RPA_D_Open.BoundingBox);
    % 5: RAD51 (this will be used to split up into monomers vs. dimers)
    RAD51_Cluster = logical(DNA(2,:) == 51);
    bb_RAD51 = regionprops(RAD51_Cluster,'BoundingBox'); bb_RAD51 = cat(1,bb_RAD51.BoundingBox);
    
    %Determine the number of proteins bound together in each cluster
    RPA_A_Bound_NumProteins = bb_RPA_A_Bound(:,3)./n_A; RPA_A_Open_NumProteins = bb_RPA_A_Open(:,3)./n_A;
    RPA_D_Bound_NumProteins = bb_RPA_D_Bound(:,3)./n_D; RPA_D_Open_NumProteins = bb_RPA_D_Open(:,3)./n_D;
    RAD51_NumProteins = bb_RAD51(:,3)./n_RAD51;
    %RAD51 Clusters which are large enough to contain dimers (same
    %structure as the BoundingBox data). Also renaming the dataset for just
    %monomers.
    bb_RAD51_Dim = bb_RAD51(bb_RAD51(3) >= 2*n_RAD51,:);
    bb_RAD51_Mon = bb_RAD51;
    
    %Now determine where each type of protein is bound (or hinged open)
    %along the lattice.
end