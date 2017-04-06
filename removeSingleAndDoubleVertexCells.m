function removeSingleAndDoubleVertexCells(param, time)
% remove cells with 1 or 2 vertices, and fix the neighbourhood.

    % there 69 1 cells
    % there 19197 2 cells
    % no other cells
    count1 = 0;
    count2 = 0;
    count3 = 0;
    count4 = 0;
    
    disp(['timestep = ' num2str(time)])
    %         cc = Cells(t, :);
    %         cc(cc == 0) = [];
    count = 0;

    %% remove 1 cell and 2 cell
    for id = 1:Ntot

        % get cell life times
        time1 = param.tracks(id).t;

        % get frame for cell (this does not correpond to frame in experiment)
        t1 = find(time1 == time);

        % get cell position, neighbours and vertices on this frame
        % [xc, yc, neighs, vertices, num_n, num_v, index] = get_cent_neigh_verts(param1, id, t1);

        xc = double(param.tracks(id).cent(t1, 1));
        yc = double(param.tracks(id).cent(t1, 2));
        neighs = cell2mat(param.tracks(id).neighs(t1));
        vertices = cell2mat(param.tracks(id).vertices(t1));

        [h, num_v] = size(vertices);
        num_n = length(neighs);

        % this variable will change if we have 1 cell or 2 cell or bad cell
        case1 = 0;

        %% 3 cases: 1 cell, 2 cell or different number of verts and neighs (the third case seems to not happen)
        %% case 1
        if(length(neighs) == 1)
            disp([num2str(time), ' : 1 cell: ', num2str(id), ' ' , num2str(num_n), ' ' , num2str(num_v)])
            % remove the one cell from param and neighbours
            % param1 = null_cell(param1, id, t1);

            param.tracks(id).cent(t1, 1) = 0;
            param.tracks(id).cent(t1, 2) = 0;
            param.tracks(id).neighs(t1) = mat2cell(0, 1);
            param.tracks(id).vertices(t1) = mat2cell([0; 0], 2);
            param.tracks(id).A(t1) = 0;
            param.tracks(id).perim(t1) = 0;

            % remove this cell from the neighbours
            cell1 = neighs;
            if(cell1 > 0)
                % param1 = remove_bad_cell(id, param1, cell1, t);

                % look at cell1 and neighs of cell1 to remove errors related to the one cell only
                time2 = param.tracks(cell1).t;
                t2 = find(time2 == time);

                % [xc1, yc1, neighs1, vertices1, num_n1, num_v1, index1] = get_cent_neigh_verts(param1, cell1, t2);

                xc1 = double(param.tracks(cell1).cent(t2, 1));
                yc1 = double(param.tracks(cell1).cent(t2, 2));
                neighs1 = cell2mat(param.tracks(cell1).neighs(t2));
                vertices1 = cell2mat(param.tracks(cell1).vertices(t2));

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
                    param.tracks(cell1).neighs(t2) = mat2cell(neighs1, 1);
                    param.tracks(cell1).vertices(t2) = mat2cell(vertices1, 2);

                    % go through the neighbours of the neighbour
                    for kk = 1:length(neighs1)

                        cell2 = neighs1(kk);
                        if(cell2 > 0)
                            time3 = param.tracks(cell2).t;
                            t3 = find(time3 == time);

                            xc2 = double(param.tracks(cell2).cent(t3, 1));
                            yc2 = double(param.tracks(cell2).cent(t3, 2));
                            neighs2 = cell2mat(param.tracks(cell2).neighs(t3));
                            vertices2 = cell2mat(param.tracks(cell2).vertices(t3));

                            if(length(neighs2) > 0) && (sum(neighs2 ~= 0) ~= 0)
                                ind2 = find(neighs2 == id);

                                % remove one cell from cell2 if it is there
                                if(length(ind2) > 0)
                                    [neighs2, vertices2] = remove_neigh(neighs2, vertices2, ind2);
                                end

                                % put neighs1 and vertices1 back in
                                param.tracks(cell2).neighs(t3) = mat2cell(neighs2, 1);
                                param.tracks(cell2).vertices(t3) = mat2cell(vertices2, 2);
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
            disp([num2str(time), ' : 2 cell: ', num2str(id), ' ' , num2str(num_n), ' ' , num2str(num_v), ' ' , num2str(neighs(2) - neighs(1))])

            % remove the one cell from param and neighbours
            % param1 = null_cell(param1, id, t1);

            param.tracks(id).cent(t1, 1) = 0;
            param.tracks(id).cent(t1, 2) = 0;
            param.tracks(id).neighs(t1) = mat2cell(0, 1);
            param.tracks(id).vertices(t1) = mat2cell([0; 0], 2);
            param.tracks(id).A(t1) = 0;
            param.tracks(id).perim(t1) = 0;

            % remove this cell from the neighbours
            cell1 = neighs(1);
            if(cell1 > 0)

                % param1 = remove_bad_cell(id, param1, cell1, t);
                % look at cell1 and neighs of cell1 to remove errors related to the two cell only
                time2 = param.tracks(cell1).t;
                t2 = find(time2 == time);
                % [xc1, yc1, neighs1, vertices1, num_n1, num_v1, index1] = get_cent_neigh_verts(param1, cell1, t2);

                xc1 = double(param.tracks(cell1).cent(t2, 1));
                yc1 = double(param.tracks(cell1).cent(t2, 2));
                neighs1 = cell2mat(param.tracks(cell1).neighs(t2));
                vertices1 = cell2mat(param.tracks(cell1).vertices(t2));

                if(length(neighs1) > 0) && (sum(neighs1 ~= 0) ~= 0)

                    % now find the one cell in neighs1
                    ind1 = find(neighs1 == id);

                    % remove one cell from cell1 if it is there
                    if(length(ind1) > 0)
                        [neighs1, vertices1] = remove_neigh(neighs1, vertices1, ind1);
                    end

                    % put neighs1 and vertices1 back in
                    param.tracks(cell1).neighs(t2) = mat2cell(neighs1, 1);
                    param.tracks(cell1).vertices(t2) = mat2cell(vertices1, 2);

                    % go through the neighbours of the neighbour
                    for kk = 1:length(neighs1)
                        cell2 = neighs1(kk);
                        if(cell2 > 0)
                            time3 = param.tracks(cell2).t;
                            t3 = find(time3 == time);
                            % [xc2, yc2, neighs2, vertices2, num_n2, num_v2, index2] = get_cent_neigh_verts(param1, cell2, t3);
                            xc2 = double(param.tracks(cell2).cent(t3, 1));
                            yc2 = double(param.tracks(cell2).cent(t3, 2));
                            neighs2 = cell2mat(param.tracks(cell2).neighs(t3));
                            vertices2 = cell2mat(param.tracks(cell2).vertices(t3));

                            if(length(neighs2) > 0) && (sum(neighs2 ~= 0) ~= 0)
                                ind2 = find(neighs2 == id);

                                % remove one cell from cell2 if it is there
                                if(length(ind2) > 0)
                                    [neighs2, vertices2] = remove_neigh(neighs2, vertices2, ind2);
                                end
                            end
                            % put neighs1 and vertices1 back in
                            param.tracks(cell2).neighs(t3) = mat2cell(neighs2, 1);
                            param.tracks(cell2).vertices(t3) = mat2cell(vertices2, 2);
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
            disp([num2str(time), ' : bad cell: ', num2str(id), ' ' , num2str(num_n), ' ' , num2str(num_v)])

            % remove the info saved in id
            param.tracks(id).cent(t1, 1) = 0;
            param.tracks(id).cent(t1, 2) = 0;
            param.tracks(id).neighs(t1) = mat2cell(0, 1);
            param.tracks(id).vertices(t1) = mat2cell([0; 0], 2);
            param.tracks(id).A(t1) = 0;
            param.tracks(id).perim(t1) = 0;

            % remove this cell from the neighbours
            for cell1 = neighs
                if(cell1 > 0)
                    time2 = param.tracks(cell1).t;
                    t2 = find(time2 == time);

                    xc1 = double(param.tracks(cell1).cent(t2, 1));
                    yc1 = double(param.tracks(cell1).cent(t2, 2));
                    neighs1 = cell2mat(param.tracks(cell1).neighs(t2));
                    vertices1 = cell2mat(param.tracks(cell1).vertices(t2));

                    if(length(neighs1) > 0) && (sum(neighs1 ~= 0) ~= 0)
                        % now find the one cell in neighs1
                        ind1 = find(neighs1 == id);

                        % remove one cell from cell1 if it is there
                        if(length(ind1) > 0)
                            [neighs1, vertices1] = remove_neigh(neighs1, vertices1, ind1);
                        end

                        % put neighs1 and vertices1 back in
                        param.tracks(cell1).neighs(t2) = mat2cell(neighs1, 1);
                        param.tracks(cell1).vertices(t2) = mat2cell(vertices1, 2);

                        % go through the neighbours of the neighbour
                        for kk = 1:length(neighs1)

                            cell2 = neighs1(kk);
                            if(cell2 > 0)
                                time3 = param.tracks(cell2).t;
                                t3 = find(time3 == time);

                                xc2 = double(param.tracks(cell2).cent(t3, 1));
                                yc2 = double(param.tracks(cell2).cent(t3, 2));
                                neighs2 = cell2mat(param.tracks(cell2).neighs(t3));
                                vertices2 = cell2mat(param.tracks(cell2).vertices(t3));

                                if(length(neighs2) > 0) && (sum(neighs2 ~= 0) ~= 0)
                                    ind2 = find(neighs2 == id);

                                    % remove one cell from cell2 if it is there
                                    if(length(ind2) > 0)
                                        [neighs2, vertices2] = remove_neigh(neighs2, vertices2, ind2);
                                    end

                                    % put neighs1 and vertices1 back in
                                    param.tracks(cell2).neighs(t3) = mat2cell(neighs2, 1);
                                    param.tracks(cell2).vertices(t3) = mat2cell(vertices2, 2);
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
        time1 = param.tracks(id).t;

        % get frame for cell (this does not correpond to frame in experiment)
        t1 = find(time1 == time);

        % get cell position, neighbours and vertices on this frame
        xc = double(param.tracks(id).cent(t1, 1));
        yc = double(param.tracks(id).cent(t1, 2));
        neighs = cell2mat(param.tracks(id).neighs(t1));
        vertices = cell2mat(param.tracks(id).vertices(t1));
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

            disp([num2str(time), ' : repeat neigh: ', num2str(id), ' ' , num2str(num_n), ' ' , num2str(num_v), ' ' , num2str(length(index))])

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

                param.tracks(id).neighs(t1) = mat2cell(neighs, 1);
                param.tracks(id).vertices(t1) = mat2cell(vertices, 2);

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