function [half_edges, verts, NHE, HE_C, C] = createHalfEdge(t, param, Cells, half_edges, verts, NHE, HE_C, C)
    
    %
    cc = Cells(t, :);
    cc(cc == 0) = [];

    he_count = 0;
    c_count = 0;
    HE = zeros(10000, 8);
    V  = zeros(10000, 4);

    % go thru cells and define half edges
    for id = cc
        if(mod(c_count, 1000) == 0)
            disp(['a: ', num2str(t), ' ', num2str(c_count/length(cc))])
        end
        c_count = c_count + 1;

        % get everything about the cell on the given time step
        % get cell life times
        time1 = param.tracks(id).t;

        % get frame for cell (this does not correpond to frame in experiment)
        t1 = find(time1 == t);
        xc = double(param.tracks(id).cent(t1, 1));
        yc = double(param.tracks(id).cent(t1, 2));
        neighs   = cell2mat(param.tracks(id).neighs(t1));
        vertices = cell2mat(param.tracks(id).vertices(t1));

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

            % replace cell with ghost cell (1) if not in cc (not in box) (this takes a while)
            for j = 1:length(neighs)
                if(ismember(neighs(j), cc) == 0)
                    neighs(j) = -1;
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

    % define ghost cell (index of the ghost cell MUST BE 1)
    C(1, 2) = he_count + 1; % put one half edge in
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

        if(c2 == 1)
            % put he into ghost cell
            ghost_he = ghost_he + 1;
            he_count = he_count + 1;
            HE(he_count, 1) = he_count;
            HE(he_count, 2) = 1;
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
    half_edges(t).HE = uint32(HE);
    verts(t).V = uint16(V);
    NHE(t) = num_he;
    HE_C(t) = he_count;
    
end