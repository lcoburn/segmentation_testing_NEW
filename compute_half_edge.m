%% Add the repositories
repPath = [filesep 'mDrives' filesep 'storage4' filesep 'Guillermo' filesep 'guillermo_functions'];
repPath = {repPath, [filesep 'mDrives' filesep 'storage4' filesep 'Guillermo' filesep 'software_repository']};

for rep=repPath
    if ispc
        addpath(genpath(['X:\' rep{1}]));
    elseif isunix
        addpath(genpath(['mnt/LocalStorage' rep{1}]));
    end
end

%%
tic
% take 6192 s in total
clearvars -except 'param' 'param1' 'Cells' 'count_t' 'cells' 'cx' 'cy' 'bx' 'by' 'c1' 'half_edges' 'C' 'NHE' 'HE_C' 'verts' 'half_edges1' 'verts1' 'num_he'
already_loaded = 0;
close all;
clc;

%% load the experimental data
% this takes 2 mins
if(already_loaded == 0)
    %     tracklet_file_path='X:\mDrives\storage4\Guillermo\guillermo_functions\half_edge_correction\test_data_set_segmentation_parameters_vertices_1_200_';   %filename without '0_.mat'
    %     param = load_tracking_param_from_pieces(tracklet_file_path);
    load 'Y:\mDrives\storage4\Guillermo\segmentation_testing\output\testOutputTracking\param.mat'
%     load param.mat
end
already_loaded = 1;

if(1 == 1)
    param1 = param;
end
% param1 = param2;
% save param1
%% load param1
if(1 == 0)
    pieces = 10;
    param1.tracks = {};
    for i = 1:275
        i
        %         load(['mat_files/test/param_piece_' num2str(i) '.mat']) % , 'A', '-v7.3'
        for j = 1:length(A)
            param1.tracks((i - 1)*pieces + j).t = A(j).t;
            param1.tracks((i - 1)*pieces + j).A = A(j).A;
            param1.tracks((i - 1)*pieces + j).cent = A(j).cent;
            param1.tracks((i - 1)*pieces + j).ellipse = A(j).ellipse;
            param1.tracks((i - 1)*pieces + j).perim = A(j).perim;
            param1.tracks((i - 1)*pieces + j).birth = A(j).birth;
            param1.tracks((i - 1)*pieces + j).death = A(j).death;
            param1.tracks((i - 1)*pieces + j).daughters = A(j).daughters;
            param1.tracks((i - 1)*pieces + j).neighs = A(j).neighs;
            param1.tracks((i - 1)*pieces + j).vertices = A(j).vertices;
        end
    end
end
%%
Ns = max([param1.tracks.t]);
Ntot = length(param1.tracks);
plotting = 1;

%% REMOVE CELLS: find 1 cell, 2 cell and cells with weird number of neighbours
% takes 226 s to go through all cells
% takes 1297 s to go through all cells and correct errors (20 mins)
if(1 == 1)
    % there 69 1 cells
    % there 19197 2 cells
    % no other cells
    count1 = 0;
    count2 = 0;
    count3 = 0;
    count4 = 0;
    
    for t = 1:Ns
        disp(['timestep = ' num2str(t)])
        %         cc = Cells(t, :);
        %         cc(cc == 0) = [];
        count = 0;
        
        %% remove 1 cell and 2 cell
        for id = 1:Ntot
            
            % get cell life times
            time1 = param1.tracks(id).t;
            
            % get frame for cell (this does not correpond to frame in experiment)
            t1 = find(time1 == t);
            
            % get cell position, neighbours and vertices on this frame
            % [xc, yc, neighs, vertices, num_n, num_v, index] = get_cent_neigh_verts(param1, id, t1);
            
            xc = double(param1.tracks(id).cent(t1, 1));
            yc = double(param1.tracks(id).cent(t1, 2));
            neighs = cell2mat(param1.tracks(id).neighs(t1));
            vertices = cell2mat(param1.tracks(id).vertices(t1));
            
            [h, num_v] = size(vertices);
            num_n = length(neighs);
            
            % this variable will change if we have 1 cell or 2 cell or bad cell
            case1 = 0;
            
            %% 3 cases: 1 cell, 2 cell or different number of verts and neighs (the third case seems to not happen)
            %% case 1
            if(length(neighs) == 1)
                disp([num2str(t), ' : 1 cell: ', num2str(id), ' ' , num2str(num_n), ' ' , num2str(num_v)])
                % remove the one cell from param and neighbours
                % param1 = null_cell(param1, id, t1);
                
                param1.tracks(id).cent(t1, 1) = 0;
                param1.tracks(id).cent(t1, 2) = 0;
                param1.tracks(id).neighs(t1) = mat2cell(0, 1);
                param1.tracks(id).vertices(t1) = mat2cell([0; 0], 2);
                param1.tracks(id).A(t1) = 0;
                param1.tracks(id).perim(t1) = 0;
                
                % remove this cell from the neighbours
                cell1 = neighs;
                if(cell1 > 0)
                    % param1 = remove_bad_cell(id, param1, cell1, t);
                    
                    % look at cell1 and neighs of cell1 to remove errors related to the one cell only
                    time2 = param1.tracks(cell1).t;
                    t2 = find(time2 == t);
                    
                    % [xc1, yc1, neighs1, vertices1, num_n1, num_v1, index1] = get_cent_neigh_verts(param1, cell1, t2);
                    
                    xc1 = double(param1.tracks(cell1).cent(t2, 1));
                    yc1 = double(param1.tracks(cell1).cent(t2, 2));
                    neighs1 = cell2mat(param1.tracks(cell1).neighs(t2));
                    vertices1 = cell2mat(param1.tracks(cell1).vertices(t2));
                    
                    % id, param1, cell1, t
                    % param1 = remove_bad_cell(id, param1, cell1, t);
                    % now find the one cell in neighs1
                    
                    if(length(neighs1) > 0) && (sum(neighs1 ~= 0) ~= 0)
                        ind1 = find(neighs1 == id);
                        
                        % remove one cell from cell1 if it is there
                        if(length(ind1) > 0)
                            [neighs1, vertices1] = remove_neigh(neighs1, vertices1, ind1);
                        end
                        
                        % put neighs1 and vertices1 back in
                        param1.tracks(cell1).neighs(t2) = mat2cell(neighs1, 1);
                        param1.tracks(cell1).vertices(t2) = mat2cell(vertices1, 2);
                        
                        % go through the neighbours of the neighbour
                        for kk = 1:length(neighs1)
                            
                            cell2 = neighs1(kk);
                            if(cell2 > 0)
                                time3 = param1.tracks(cell2).t;
                                t3 = find(time3 == t);
                                
                                xc2 = double(param1.tracks(cell2).cent(t3, 1));
                                yc2 = double(param1.tracks(cell2).cent(t3, 2));
                                neighs2 = cell2mat(param1.tracks(cell2).neighs(t3));
                                vertices2 = cell2mat(param1.tracks(cell2).vertices(t3));
                                
                                if(length(neighs2) > 0) && (sum(neighs2 ~= 0) ~= 0)
                                    ind2 = find(neighs2 == id);
                                    
                                    % remove one cell from cell2 if it is there
                                    if(length(ind2) > 0)
                                        [neighs2, vertices2] = remove_neigh(neighs2, vertices2, ind2);
                                    end
                                    
                                    % put neighs1 and vertices1 back in
                                    param1.tracks(cell2).neighs(t3) = mat2cell(neighs2, 1);
                                    param1.tracks(cell2).vertices(t3) = mat2cell(vertices2, 2);
                                end
                            end
                        end
                    end
                end
                
                case1 = 1;
                count = count + 1;
                count1 = count1 + 1;
            end
            
            %% case 2
            if(length(neighs) == 2)
                
                [h, k] = size(vertices);
                disp([num2str(t), ' : 2 cell: ', num2str(id), ' ' , num2str(num_n), ' ' , num2str(num_v), ' ' , num2str(neighs(2) - neighs(1))])
                
                % remove the one cell from param and neighbours
                % param1 = null_cell(param1, id, t1);
                
                param1.tracks(id).cent(t1, 1) = 0;
                param1.tracks(id).cent(t1, 2) = 0;
                param1.tracks(id).neighs(t1) = mat2cell(0, 1);
                param1.tracks(id).vertices(t1) = mat2cell([0; 0], 2);
                param1.tracks(id).A(t1) = 0;
                param1.tracks(id).perim(t1) = 0;
                
                % remove this cell from the neighbours
                cell1 = neighs(1);
                if(cell1 > 0)
                    
                    % param1 = remove_bad_cell(id, param1, cell1, t);
                    % look at cell1 and neighs of cell1 to remove errors related to the two cell only
                    time2 = param1.tracks(cell1).t;
                    t2 = find(time2 == t);
                    % [xc1, yc1, neighs1, vertices1, num_n1, num_v1, index1] = get_cent_neigh_verts(param1, cell1, t2);
                    
                    xc1 = double(param1.tracks(cell1).cent(t2, 1));
                    yc1 = double(param1.tracks(cell1).cent(t2, 2));
                    neighs1 = cell2mat(param1.tracks(cell1).neighs(t2));
                    vertices1 = cell2mat(param1.tracks(cell1).vertices(t2));
                    
                    if(length(neighs1) > 0) && (sum(neighs1 ~= 0) ~= 0)
                        
                        % now find the one cell in neighs1
                        ind1 = find(neighs1 == id);
                        
                        % remove one cell from cell1 if it is there
                        if(length(ind1) > 0)
                            [neighs1, vertices1] = remove_neigh(neighs1, vertices1, ind1);
                        end
                        
                        % put neighs1 and vertices1 back in
                        param1.tracks(cell1).neighs(t2) = mat2cell(neighs1, 1);
                        param1.tracks(cell1).vertices(t2) = mat2cell(vertices1, 2);
                        
                        % go through the neighbours of the neighbour
                        for kk = 1:length(neighs1)
                            cell2 = neighs1(kk);
                            if(cell2 > 0)
                                time3 = param1.tracks(cell2).t;
                                t3 = find(time3 == t);
                                % [xc2, yc2, neighs2, vertices2, num_n2, num_v2, index2] = get_cent_neigh_verts(param1, cell2, t3);
                                xc2 = double(param1.tracks(cell2).cent(t3, 1));
                                yc2 = double(param1.tracks(cell2).cent(t3, 2));
                                neighs2 = cell2mat(param1.tracks(cell2).neighs(t3));
                                vertices2 = cell2mat(param1.tracks(cell2).vertices(t3));
                                
                                if(length(neighs2) > 0) && (sum(neighs2 ~= 0) ~= 0)
                                    ind2 = find(neighs2 == id);
                                    
                                    % remove one cell from cell2 if it is there
                                    if(length(ind2) > 0)
                                        [neighs2, vertices2] = remove_neigh(neighs2, vertices2, ind2);
                                    end
                                end
                                % put neighs1 and vertices1 back in
                                param1.tracks(cell2).neighs(t3) = mat2cell(neighs2, 1);
                                param1.tracks(cell2).vertices(t3) = mat2cell(vertices2, 2);
                            end
                        end
                        case1 = 2;
                        count = count + 1;
                        count2 = count2 + 1;
                    end
                end
            end
            
            %% case 3
            if(num_n ~= num_v) && (case1 == 0)
                disp([num2str(t), ' : bad cell: ', num2str(id), ' ' , num2str(num_n), ' ' , num2str(num_v)])
                
                % remove the info saved in id
                param1.tracks(id).cent(t1, 1) = 0;
                param1.tracks(id).cent(t1, 2) = 0;
                param1.tracks(id).neighs(t1) = mat2cell(0, 1);
                param1.tracks(id).vertices(t1) = mat2cell([0; 0], 2);
                param1.tracks(id).A(t1) = 0;
                param1.tracks(id).perim(t1) = 0;
                
                % remove this cell from the neighbours
                for cell1 = neighs
                    if(cell1 > 0)
                        time2 = param1.tracks(cell1).t;
                        t2 = find(time2 == t);
                        
                        xc1 = double(param1.tracks(cell1).cent(t2, 1));
                        yc1 = double(param1.tracks(cell1).cent(t2, 2));
                        neighs1 = cell2mat(param1.tracks(cell1).neighs(t2));
                        vertices1 = cell2mat(param1.tracks(cell1).vertices(t2));
                        
                        if(length(neighs1) > 0) && (sum(neighs1 ~= 0) ~= 0)
                            % now find the one cell in neighs1
                            ind1 = find(neighs1 == id);
                            
                            % remove one cell from cell1 if it is there
                            if(length(ind1) > 0)
                                [neighs1, vertices1] = remove_neigh(neighs1, vertices1, ind1);
                            end
                            
                            % put neighs1 and vertices1 back in
                            param1.tracks(cell1).neighs(t2) = mat2cell(neighs1, 1);
                            param1.tracks(cell1).vertices(t2) = mat2cell(vertices1, 2);
                            
                            % go through the neighbours of the neighbour
                            for kk = 1:length(neighs1)
                                
                                cell2 = neighs1(kk);
                                if(cell2 > 0)
                                    time3 = param1.tracks(cell2).t;
                                    t3 = find(time3 == t);
                                    
                                    xc2 = double(param1.tracks(cell2).cent(t3, 1));
                                    yc2 = double(param1.tracks(cell2).cent(t3, 2));
                                    neighs2 = cell2mat(param1.tracks(cell2).neighs(t3));
                                    vertices2 = cell2mat(param1.tracks(cell2).vertices(t3));
                                    
                                    if(length(neighs2) > 0) && (sum(neighs2 ~= 0) ~= 0)
                                        ind2 = find(neighs2 == id);
                                        
                                        % remove one cell from cell2 if it is there
                                        if(length(ind2) > 0)
                                            [neighs2, vertices2] = remove_neigh(neighs2, vertices2, ind2);
                                        end
                                        
                                        % put neighs1 and vertices1 back in
                                        param1.tracks(cell2).neighs(t3) = mat2cell(neighs2, 1);
                                        param1.tracks(cell2).vertices(t3) = mat2cell(vertices2, 2);
                                    end
                                end
                            end
                            
                        end
                        
                        case1 = 3;
                        count = count + 1;
                        count3 = count3 + 1;
                    end
                end
            end
            
        end
        disp('---------------------')
        
        %% remove repeated neighs from cells
        for id = 1:Ntot
            % get cell life times
            time1 = param1.tracks(id).t;
            
            % get frame for cell (this does not correpond to frame in experiment)
            t1 = find(time1 == t);
            
            % get cell position, neighbours and vertices on this frame
            xc = double(param1.tracks(id).cent(t1, 1));
            yc = double(param1.tracks(id).cent(t1, 2));
            neighs = cell2mat(param1.tracks(id).neighs(t1));
            vertices = cell2mat(param1.tracks(id).vertices(t1));
            %             length(neighs) - length(vertices)
            
            [h, num_v] = size(vertices);
            num_n = length(neighs);
            
            u = unique(neighs);
            nn = histc(neighs, u);
            
            if(find(nn > 1))
                index = find(nn > 1);
            else
                index = 0;
            end
            
            % repeated cells
            if(index ~= 0)
                
                disp([num2str(t), ' : repeat neigh: ', num2str(id), ' ' , num2str(num_n), ' ' , num2str(num_v), ' ' , num2str(length(index))])
                
                % index
                count4 = count4 + 1;
                u = unique(neighs);
                
                for j = 1:length(index)
                    cell1 = u(index(j));
                    % [param1, neighs, vertices] = remove_second_to_last_one(neighs, vertices, cell1, param1, t1, id);
                    ind = find(neighs == cell1);
                    ind1 = ind;
                    
                    if(length(ind) == 2)
                        if(ind(2) - ind(1) == 1)
                            ind1(1) = [];
                        else
                            ind1(2) = [];
                        end
                    end
                    
                    if(length(ind) == 3)
                        if( (ind(2) - ind(1) == 1) && (ind(3) - ind(2) == 1) )
                            ind1(1) = [];
                        end
                        if( (ind(2) - ind(1) > 1) && (ind(3) - ind(2) == 1) )
                            ind1(2) = [];
                        end
                        if( (ind(2) - ind(1) == 1) && (ind(3) - ind(2) > 1) )
                            ind1(3) = [];
                        end
                    end
                    % ind1
                    neighs(ind1) = [];
                    vertices(:, ind1) = [];
                    
                    param1.tracks(id).neighs(t1) = mat2cell(neighs, 1);
                    param1.tracks(id).vertices(t1) = mat2cell(vertices, 2);
                    
                    % plot
                    if(1 == 0)
                        x = vertices(1,:);
                        y = vertices(2,:);
                        plot([x,x(1)], [y, y(1)])
                    end
                    
                end
            end
        end
        disp('******************')
        
    end
end

%% GET CELLS: exclude cells outside bounds
% this takes 1003 s (~20 mins)
if(1 == 1)
    %     load('bounds1.mat');
    bounds = [0 0; 0 600; 600 600; 600 0]
    % these bounds have been created for F38_237
    % this should be redone for other experiments
    bx = bounds(:, 1);
    by = bounds(:, 2);
    
    count_t = zeros(1,Ns);
    Cells = zeros(Ns, Ntot);
    for id = 1: Ntot
        if(mod(id,100) == 0);
            disp(id/Ntot);
        end
        time = param.tracks(id).t;
        for t = 1:length(time)
            xc = double(param.tracks(id).cent(t, 1));
            yc = double(param.tracks(id).cent(t, 2));
            in = inpolygon(xc, yc, bx, by);
            
            if(in == 0)
                % remove the cell from its neighbour and the corresponding vertex
                %                 neighs = cell2mat(param1.tracks(id).neighs(t));
                %
                %                 for i = 1:length(neighs)
                %                     if(neighs(i) ~= 0)
                %                         time1 = param.tracks(neighs(i)).t;
                %                         t1 = find(time1 == t);
                %                         neighs1 = cell2mat(param1.tracks(neighs(i)).neighs(t1));
                %                         vertices1 = cell2mat(param1.tracks(id).vertices(t1));
                %                         ind1 = id;
                %                         [neighs1, vertices1] = remove_neigh(neighs1, vertices1, ind1);
                %
                %                         % put neighs and vertices back in
                %                         param1.tracks(id).neighs(t1) = mat2cell(neighs1, 1);
                %                         param1.tracks(id).vertices(t1) = mat2cell(vertices1, 2);
                %                     end
                %                 end
                
                % remove the info saved in id
                %                 param1.tracks(id).cent(t, 1) = 0;
                %                 param1.tracks(id).cent(t, 2) = 0;
                %                 param1.tracks(id).neighs(t) = mat2cell(0, 1);
                %                 param1.tracks(id).vertices(t) = mat2cell([0; 0], 2);
                %                 param1.tracks(id).A(t) = 0;
                %                 param1.tracks(id).perim(t) = 0;
            else
                count_t(time(t)) = count_t(time(t)) + 1;
                Cells(time(t), count_t(time(t))) = id;
            end
            
        end
    end
    save('cells_in_box_during_each_frame_1.mat', 'Cells', 'count_t')
    % this is for experiment F38_237
    %     param2 = param1;
    % reset if any issues
end

%% save param1 in 270 pieces
if(1 == 0)
    pieces = 1000;
    % D = {};
    for i = 1:ceil(length(param1.tracks)/pieces)
        i
        if(i*pieces < length(param1.tracks))
            A = param1.tracks((i - 1)*pieces + 1: i*pieces);
            %         D((i - 1)*pieces + 1: i*pieces) = A;
            length(A)
        else
            A = param1.tracks((i - 1)*pieces + 1: end);
            %         D((i - 1)*pieces + 1: (i - 1)*pieces + length(A)) = A;
            length(A)
        end
        save(['mat_files/test/param_piece_' num2str(i) '.mat'], 'A', '-v7.3')
    end
end

%% CREATE HALF EDGE DATA STRUCTURE
if(1 == 1)
    load('cells_in_box_during_each_frame_1.mat') % , 'Cells', 'count_t'
    % this takes 4303 s
    if(1 == 1)
        half_edges = {};
        verts = {};
        NHE = zeros(1, Ns);
        HE_C = zeros(1, Ns);
        C = zeros(Ntot + 1, Ns);
        
        for t = 1:Ns
            t
            cc = Cells(t, :);
            cc(cc == 0) = [];
            
            he_count = 0;
            c_count = 0;
            HE = zeros(1000000, 8);
            V  = zeros(1000000, 4);
            
            % go thru cells and define half edges
            for id = cc
                if(mod(c_count, 1000) == 0)
                    disp(['a: ', num2str(t), ' ', num2str(c_count/length(cc))])
                end
                c_count = c_count + 1;
                
                % get everything about the cell on the given time step
                % get cell life times
                time1 = param1.tracks(id).t;
                
                % get frame for cell (this does not correpond to frame in experiment)
                t1 = find(time1 == t);
                xc = double(param1.tracks(id).cent(t1, 1));
                yc = double(param1.tracks(id).cent(t1, 2));
                neighs   = cell2mat(param1.tracks(id).neighs(t1));
                vertices = cell2mat(param1.tracks(id).vertices(t1));
                
                if(length(neighs) > 2) % just in case
                    
                    % sort verts in ACW direction
                    x = double(vertices(1, :));
                    y = double(vertices(2, :));
                    angle = atan2(y - mean(y), x - mean(x));
                    
                    cw_acw = 0;
                    for j = 1:length(angle)
                        j1 = mod(j, length(angle)) + 1;
                        if(angle(j1) - angle(j) > 0)
                            cw_acw = cw_acw + 1;
                        else
                            cw_acw = cw_acw - 1;
                        end
                    end
                    
                    if(cw_acw < 0)
                        neighs = fliplr(neighs);
                        vertices = fliplr(vertices);
                        x = double(vertices(1, :));
                        y = double(vertices(2, :));
                        % perhaps put back into param1 here ??
                    end
                    
                    % replace cell with ghost cell (Ntot + 1) if not in cc (not in box) (this takes a while)
                    for j = 1:length(neighs)
                        if(ismember(neighs(j), cc) == 0)
                            neighs(j) = Ntot + 1;
                        end
                    end
                    
                    % put data in half edge and cell arrays
                    % fill in half edge data
                    he   = he_count + 1: he_count + length(neighs);
                    ind  = 1:length(he);
                    ind1 = mod(ind, length(he)) + 1;
                    he1  = he(ind1);
                    ind2 = mod(ind - 2, length(he)) + 1;
                    he2  = he(ind2);
                    
                    % half edge index
                    HE(he, 1) = he';
                    % cell index
                    HE(he, 2) = id*ones(length(neighs), 1);
                    % cell neighbour index
                    HE(he, 3) = neighs';
                    % next half edge
                    HE(he, 4) = he1';
                    % prev half edge
                    HE(he, 5) = he2';
                    
                    % vertices
                    V(he, 1:4) = [vertices' , vertices(:, ind2)'];
                    
                    % fill in data into cells
                    C(id, t) = he(1);                       % put one half edge in
                    he_count = he_count + length(neighs);
                end
            end
            
            % define ghost cell
            C(Ntot + 1, 2) = he_count + 1; % put one half edge in
            num_he = he_count;
            ghost_he = 0;
            
            % find the opposite (this is giving trouble)
            % go thru half edges
            for i = 1:num_he
                if(mod(i, 1000) == 0)
                    disp(['b: ', num2str(t), ' ', num2str(i/num_he)])
                end
                c1 = HE(i, 2);
                c2 = HE(i, 3);
                
                if(c2 == Ntot + 1)
                    % put he into ghost cell
                    ghost_he = ghost_he + 1;
                    he_count = he_count + 1;
                    HE(he_count, 1) = he_count;
                    HE(he_count, 2) = Ntot + 1;
                    HE(he_count, 3) = c1;
                    HE(he_count, 4) = he_count + 1;
                    HE(he_count, 5) = he_count - 1;
                    HE(he_count, 6) = HE(i, 1);
                    HE(i, 6) = he_count;
                    V(he_count, 1:4) = V(i, [3, 4, 1, 2]);
                    
                else
                    
                    he1 = C(c2, t); % a half edge belonging to c2
                    he1_start = he1;
                    if(he1 ~= 0)
                        
                        c1a = HE(he1, 2);
                        c2a = HE(he1, 3);
                        count = 0;
                        if(c2a == c1)
                            
                            HE(i, 6) = he1;
                            HE(he1, 6) = i;
                            a = V(i, :);
                            b = V(he1, :);
                            b = [b(3:4), b(1:2)];
                            
                            c = abs(a - b);
                            % this means the vertices have been improperly assigned
                            % to the half_edges and have to be changed
                            if(sum(c) ~= 0)
                                
                                % get the indices of HE that correspond to the half edges of c1
                                i_start = i;
                                ind1 = i_start;
                                i1 = HE(i_start, 4);
                                while(i1 ~= i_start)
                                    ind1 = [ind1, i1];
                                    i1 = HE(i1, 4);
                                end
                                
                                he1_start1 = he1;
                                ind2 = he1_start1;
                                he2 = HE(he1_start1, 4);
                                while(he2 ~= he1_start1)
                                    ind2 = [ind2, he2];
                                    he2 = HE(he2, 4);
                                end
                                
                                ind1a = ind1;
                                ind2a = ind2;
                                a = V(ind1a(1), :);
                                b = V(ind2a(1), :);
                                b = [b(3:4), b(1:2)];
                                c = abs(a - b);
                                cont1 = 0;
                                % disp('-')
                                for ii = 0:length(ind1)
                                    ind1a = ind1(mod(ii:ii + length(ind1) - 1, length(ind1))  + 1);
                                    for jj = 1:length(ind2);
                                        ind2a = ind2(mod(jj:jj + length(ind2) - 1, length(ind2))  + 1);
                                        a = V(ind1a(1), :);
                                        b = V(ind2a(1), :);
                                        b = [b(3:4), b(1:2)];
                                        c = abs(a - b);
                                        if(sum(c) == 0)
                                            cont1 = 1;
                                            break;
                                        end
                                    end
                                    if(cont1 == 1)
                                        break;
                                    end
                                end
                                % disp('*')
                                ind1a;
                                ind2a;
                                V(ind1', :) = V(ind1a',:);
                                V(ind2', :) = V(ind2a',:);
                                
                            end
                        else
                            he1 = HE(he1, 4);
                            while( (c2a ~= c1)  && (he1 ~= he1_start) ) % && (he1 ~= he1_start)
                                
                                count = count + 1;
                                c1a = HE(he1, 2);
                                c2a = HE(he1, 3);
                                
                                if(c2a == c1)
                                    HE(i, 6) = he1;
                                    HE(he1, 6) = i;
                                    a = V(i, :);
                                    b = V(he1, :);
                                    b = [b(3:4), b(1:2)];
                                    
                                    c = abs(a - b);
                                    % this means the vertices have been improperly assigned
                                    % to the half_edges and have to be changed
                                    if(sum(c) ~= 0)
                                        
                                        % get the indices of HE that correspond to the half edges of c1
                                        i_start = i;
                                        ind1 = i_start;
                                        i1 = HE(i_start, 4);
                                        while(i1 ~= i_start)
                                            ind1 = [ind1, i1];
                                            i1 = HE(i1, 4);
                                        end
                                        
                                        he1_start1 = he1;
                                        ind2 = he1_start1;
                                        he2 = HE(he1_start1, 4);
                                        while(he2 ~= he1_start1)
                                            ind2 = [ind2, he2];
                                            he2 = HE(he2, 4);
                                        end
                                        
                                        ind1a = ind1;
                                        ind2a = ind2;
                                        a = V(ind1a(1), :);
                                        b = V(ind2a(1), :);
                                        b = [b(3:4), b(1:2)];
                                        c = abs(a - b);
                                        cont1 = 0;
                                        
                                        for ii = 0:length(ind1)
                                            ind1a = ind1(mod(ii:ii + length(ind1) - 1, length(ind1))  + 1);
                                            for jj = 1:length(ind2);
                                                ind2a = ind2(mod(jj:jj + length(ind2) - 1, length(ind2))  + 1);
                                                a = V(ind1a(1), :);
                                                b = V(ind2a(1), :);
                                                b = [b(3:4), b(1:2)];
                                                c = abs(a - b);
                                                if(sum(c) == 0)
                                                    cont1 = 1;
                                                    break;
                                                end
                                            end
                                            if(cont1 == 1)
                                                break;
                                            end
                                        end
                                        % disp('*')
                                        ind1a;
                                        ind2a;
                                        V(ind1', :) = V(ind1a',:);
                                        V(ind2', :) = V(ind2a',:);
                                        
                                    end
                                end
                                
                                he1 = HE(he1, 4);
                            end
                        end
                    end
                end
            end
            
            HE(num_he + 1, 5) = he_count;
            HE(he_count, 4) = num_he + 1;
            
            % crop HE and V to last he_count
            ind = min(find(HE(:, 1) == 0));
            HE(ind:end, :) = [];
            V(ind:end, :) = [];
            
            % put HE into struct and save number of number of half edges on frame
            half_edges(t).HE = HE;
            verts(t).V = V;
            NHE(t) = num_he;
            HE_C(t) = he_count;
        end
        %         save('test_he_4.mat', 'HE', 'V', 'C')
    end
    
    % remove any half edge without opposite
    if(1 == 1)
        for t = 1:Ns-1
            t
            % get half edge array on current and next time step
            HE = half_edges(t).HE;
            V = verts(t).V;
            I = find(HE(:, 6) == 0);
            
            if(length(I) > 0)
                for i = I'
                    next_he = HE(i, 4);
                    prev_he = HE(i, 5);
                    HE(next_he, 5) = prev_he;
                    HE(prev_he, 4) = next_he;
                    HE(i, :) = zeros(1, 8);
                end
            end
            half_edges(t).HE = HE;
        end
    end
    
    % remove any half whose opposite is 0
    if(1 == 1)
        for t = 1:Ns-1
            t
            % get half edge array on current and next time step
            HE = half_edges(t).HE;
            V = verts(t).V;
            for i = 1:length(HE)
                he1 = HE(i, 6);
                if(he1 > 0)
                    if(HE(he1, 6) == 0)
                        disp('error')
                    end
                end
            end
        end
    end
    
    %% find vertex number and create vertex data
    if(1 == 1)
        VN = {};
        
        for t = 1:Ns
            t
            
            % load HE and V
            HE = double(half_edges(t).HE);
            V = double(verts(t).V);
            num_he = NHE(t);
            
            % arrays for vertices
            HE_vertex = zeros(length(HE), 1);
            Vertices = zeros(2*num_he, 8);
            v_pos = zeros(2*num_he, 2);
            v_count = 0;
            
            % create visited half edge array
            visited = zeros(length(HE), 1);
            
            % got thru edges
            for he = 1:num_he
                c1 = HE(he, 2);
                c2 = HE(he, 3);
                edge = 0;
                
                % if not visited
                if(visited(he) == 0) && (c1 ~= Ntot + 1) && (c2 ~= Ntot + 1)
                    
                    %% plot vertex and vertex bonds
                    counter = 1;
                    % vertex co-ords
                    vx = V(he, 1);
                    vy = V(he, 2);
                    % cells
                    c1 = HE(he, 2);
                    c2 = HE(he, 3);
                    % start vertex edges array
                    e = he;
                    
                    % get opposite, previous, visited and v(1) v(2).
                    % Increment counter
                    he1 = HE(he, 6);
                    if(he1 ~= 0)
                        he1 = HE(he1, 5);
                        c1 = HE(he1, 2);
                        c2 = HE(he1, 3);
                        if(c1 == Ntot + 1) || (c2 == Ntot + 1)
                            edge = 1;
                        end
                    end
                    
                    while(he1 ~= he) && (he1 ~= 0) && (c1 ~= Ntot + 1) && (c2 ~= Ntot + 1)
                        counter = counter + 1;
                        vx = [vx, V(he1, 1)];
                        vy = [vy, V(he1, 2)];
                        v = V(he1, :);
                        e = [e, he1];
                        
                        % get opposite, previous, visited and v(3) v(4). Increment counter
                        he1 = HE(he1, 6);
                        if(he1 ~= 0)
                            he1 = HE(he1, 5);
                            c1 = HE(he1, 2);
                            c2 = HE(he1, 3);
                            if(c1 == Ntot + 1) || (c2 == Ntot + 1)
                                edge = 1;
                            end
                        end
                    end
                    
                    % put vertex into vertices array
                    if(c1 ~= Ntot + 1) && (edge == 0) && (length(e) <= 6)
                        % save number of vertices
                        v_count = v_count + 1;
                        Vertices(v_count, 1:length(e)) = e;
                        v_pos(v_count, 1) = vx(1);
                        v_pos(v_count, 2) = vy(1);
                        % assign vertex to edge and mark
                        for j = 1:length(e)
                            HE_vertex(e(j)) = v_count;
                            visited(e(j)) = 1;
                        end
                    end
                    
                end
            end
            
            % make opposite he of unvivtsed he unvisted
            for he = 1:num_he;
                if(visited(he) == 0)
                    visited(HE(he, 6)) = 0;
                end
            end
            
            Vertices(v_count + 1:end,:) = [];
            v_pos(v_count + 1:end, :) = [];
            
            % change C1 back to uint
            half_edges(t).HE = uint32(HE);
            half_edges(t).HE_vertex = uint32(HE_vertex);
            half_edges(t).visited = visited;
            verts(t).V = uint16(V);
            verts(t).Vertices = uint32(Vertices);
            verts(t).v_pos = uint16(v_pos);
        end
        
               
    end
    
    % remove short
    if(1 == 0)
        threshold = 5; % pxls
        for t = 1:Ns;
            t
            
            % load vertices and half edges
            HE = double(half_edges(t).HE);
            V = double(verts(t).V);
            HE_vertex = double(half_edges(t).HE_vertex);
            visited = half_edges(t).visited;
            Vertices = double(verts(t).Vertices);
            v_pos = double(verts(t).v_pos);
            
            num_he = NHE(t);
            good_bond = visited;
            
            for he = 1:num_he;
                
                v = V(he, :);
                c1 = HE(he, 2);
                c2 = HE(he, 3);
                len = sqrt( (v(1) - v(3))^2 + (v(2) - v(4))^2 );
                
                % get opposite
                heo = HE(he, 6);
                
                % if bond is shorted than threshold get rid of it
                if(len <= threshold) && (heo ~= 0)
                    
                    % get vertices
                    v1 = HE_vertex(he, 1);
                    v2 = HE_vertex(heo, 1);
                    
                    % dont go in unless both edges hhave assigned vertex
                    % somehow one is not being assigned ? (ignore)
                    if(v1 > 0) && (v2 > 0)
                        % remove bond and c1 and c2 if necessary
                        [HE, V, C, v_pos, Vertices, HE_vertex, good_bond] = remove_half_edge_1(he, heo, v1, v2, HE, V, C, t, Ntot, v_pos, Vertices, HE_vertex, good_bond, c1, c2);
                    end
                end
            end
            
            for he = 1:num_he;
                v = V(he, :);
                len = sqrt((v(1) - v(3))^2 + (v(2) - v(4))^2);
                heo = HE(he, 6);
                c1 = HE(he, 2);
                c2 = HE(he, 3);
                
                if(len <= threshold)
                    good_bond(he) = 0;
                end
                if(len > 40)
                    good_bond(he) = 0;
                end
                if(heo == 0)
                    good_bond(he) = 0;
                end
                if(c1 == 0) || (c1 == Ntot + 1) || (c2 == 0) || (c2 == Ntot + 1)
                    good_bond(he) = 0;
                end
                
            end
            
            if(t == 1)
                ind2 = (HE(:, 8) == 0);
                ind = find(ind2);
                good_bond(ind) = 0;
            end
            
            if(t > 1)
                ind1 = (HE(:, 7) == 0);
                ind2 = (HE(:, 8) == 0);
                ind = find(ind1.*ind2);
                good_bond(ind) = 0;
            end
            
            % change C1 back to uint
            half_edges(t).HE = uint32(HE);
            half_edges(t).HE_vertex = uint32(HE_vertex);
            half_edges(t).visited = visited;
            half_edges(t).good_bond = good_bond;
            verts(t).V = uint16(V);
            verts(t).Vertices = uint32(Vertices);
            verts(t).v_pos = uint16(v_pos);
        end
    end    
    
    %% connect frames: find next half edge index and previous index
    if(1 == 1)
        for t = 1:Ns - 1
            t
            %% get half edge array on current and next time step
            t1 = t + 1;                 % next time step
            HE = double(half_edges(t).HE);      % half edge on current time step
            HE1 = double(half_edges(t1).HE);    % half edge on next time step
            num_he = NHE(t);
            he_count = HE_C(t);
            
            %% go thru half edges and get next and previous
            for i = 1:num_he
                c1 = HE(i, 2);
                c2 = HE(i, 3);
                
                if(c1 > 0) && (c2 ~= Ntot + 1)
                    % get a half edge of c1 on next time step
                    he1 = C(c1, t1);
                    he1_start = he1;        % first value for he1
                    
                    if(he1 > 0) % this will mean that the cell is still present
                        he2 = HE1(he1, 6);        % pair he
                        
                        if(he2 > 0) && (he2 <= num_he)
                            
                            if(HE1(he2, 2) == c2)
                                % he1
                                HE(i, 8) = he1;
                                HE1(he1, 7) = i;
                                HE(HE(i, 6),8) = HE1(he1, 6);
                                HE1(HE1(he1, 6),7) = HE(i, 6);
                            end
                            
                            if(HE1(he2, 2) ~= c2)
                                he1 = HE1(he1, 4);    % next he
                                he2 = HE1(he1, 6);    % and its pair he
                                
                                while( (he2 > 0) && (HE1(he2, 2) ~= c2) && (he1 ~= he1_start) )   % while the pair cell is not c2
                                    he1 = HE1(he1, 4);    % next he
                                    he2 = HE1(he1, 6);    % and its pair he
                                end
                                
                                if(he1 ~= he1_start) && (he2 > 0)
                                    % the bond is still there
                                    % he1
                                    HE(i, 8) = he1;
                                    HE1(he1, 7) = i;
                                    HE(HE(i, 6),8) = HE1(he1, 6);
                                    HE1(HE1(he1, 6),7) = HE(i, 6);
                                    HE1(he1, [1 2 3 4 5 6 7 8]);
                                end
                                
                            end
                            
                        end
                    end
                end
            end
            
            % put half edges back in
            half_edges(t).HE = HE;
            half_edges(t1).HE = HE1;
        end
    end
    
    % convert C, HE and V to uint16 and uint32
    if(1 == 1)
        half_edges1 = {};
        verts1 = {};
        for t = 1:Ns
            half_edges1(t).HE = uint32(half_edges(t).HE);
            verts1(t).V = uint16(verts(t).V);
        end
        C1 = uint32(C);
    end
    
    % save
    if(1 == 1)
        % this takes 60 s
        save(['half_edge_data_cells_9.mat'], 'C1', 'Cells', 'Ntot', 'Ns', '-v7.3');
        save(['half_edge_data_he_9.mat'],'half_edges1','NHE','HE_C','num_he','he_count','ghost_he','-v7.3');
        save(['half_edge_data_vert_9.mat'],'verts1', '-v7.3');
    end
end

toc