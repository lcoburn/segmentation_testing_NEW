close all;
clc;
for t = 1:Ns;
        clearvars -except t threshold plotting half_edges3 verts3 NHE C1 Cells Ntot Ns HE_C num_he he_count ghost_he half_edges2 verts2
        t
        
        % load vertices and half edges
        HE = double(half_edges2(t).HE);
        V = double(verts2(t).V);
        HE_vertex = double(half_edges2(t).HE_vertex);
        visited = half_edges2(t).visited;
        Vertices = double(verts2(t).Vertices);
        v_pos = double(verts2(t).v_pos);

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
                    [HE, V, C1, v_pos, Vertices, HE_vertex, good_bond] = remove_half_edge_1(he, heo, v1, v2, HE, V, C1, t, Ntot, v_pos, Vertices, HE_vertex, good_bond, c1, c2);
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
        half_edges3(t).HE = uint32(HE);
        half_edges3(t).HE_vertex = uint32(HE_vertex);
        half_edges3(t).visited = visited;
        half_edges3(t).good_bond = good_bond;
        verts3(t).V = uint16(V);
        verts3(t).Vertices = uint32(Vertices);
        verts3(t).v_pos = uint16(v_pos);
    end
    
tic
%% load
if(1 == 1)
    clear all;
    load('half_edge_data_cells_10.mat')
    load('half_edge_data_he_10.mat')
    load('half_edge_data_vert_10.mat')
end

if(1 == 1)
    %figure('units','normalized','outerposition',[0 0 0.66 1]);
    threshold = 5;
    plotting = 0;
    half_edges3 = {};
    verts3 = {};
    
    for t = 1:Ns;
        clearvars -except t threshold plotting half_edges3 verts3 NHE C1 Cells Ntot Ns HE_C num_he he_count ghost_he half_edges2 verts2
        t
        
        % load vertices and half edges
        HE = double(half_edges2(t).HE);
        V = double(verts2(t).V);
        HE_vertex = double(half_edges2(t).HE_vertex);
        visited = half_edges2(t).visited;
        Vertices = double(verts2(t).Vertices);
        v_pos = double(verts2(t).v_pos);

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
                    [HE, V, C1, v_pos, Vertices, HE_vertex, good_bond] = remove_half_edge_1(he, heo, v1, v2, HE, V, C1, t, Ntot, v_pos, Vertices, HE_vertex, good_bond, c1, c2);
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
        half_edges3(t).HE = uint32(HE);
        half_edges3(t).HE_vertex = uint32(HE_vertex);
        half_edges3(t).visited = visited;
        half_edges3(t).good_bond = good_bond;
        verts3(t).V = uint16(V);
        verts3(t).Vertices = uint32(Vertices);
        verts3(t).v_pos = uint16(v_pos);
    end
    
    % change C1 back to uint
    C1 = uint32(C1);
    
    % save
    if(1 == 1)
        % this takes 60 s
        save(['half_edge_data_cells_13.mat'], 'C1', 'Cells', 'Ntot', 'Ns', '-v7.3');
        save(['half_edge_data_he_13.mat'], 'half_edges3', 'NHE', 'HE_C', 'num_he', 'he_count', 'ghost_he', '-v7.3');
        save(['half_edge_data_vert_13.mat'],'verts3',  '-v7.3');
    end
    
end
