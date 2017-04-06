close all;
clc;

tic
%% load
if(1 == 1)
    clear all;
    % figure('units','normalized','outerposition',[0 0 0.66 1]);
    load('half_edge_data_cells_9.mat')
    load('half_edge_data_he_9.mat')
    load('half_edge_data_vert_9.mat')
    load('param1')
    load bond_index_life_time_1.mat; %','BONDS','-v7.3');
    load bond_life_time_data_1.mat; % ','lifetime','-v7.3'
end

%% find vertex number
if(1 == 1)
    VN = {};
    plotting = 1;
    half_edges2 = {};
    verts2 = {};
    C1 = double(C1);
    
    for t = 1:Ns
        t
        
        % load HE and V
        HE = double(half_edges1(t).HE);
        V = double(verts1(t).V);
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
        half_edges2(t).HE = uint32(HE);
        half_edges2(t).HE_vertex = uint32(HE_vertex);
        half_edges2(t).visited = visited;
        verts2(t).V = uint16(V);
        verts2(t).Vertices = uint32(Vertices);
        verts2(t).v_pos = uint16(v_pos);
    end
    
    % change C1 back to uint
    C1 = uint32(C1);
    half_edges1 = half_edges2;
    verts1 = verts2;
    
    % save 
    if(1 == 1)
        % this takes 60 s
        save(['half_edge_data_cells_10.mat'], 'C1', 'Cells', 'Ntot', 'Ns', '-v7.3');
        save(['half_edge_data_he_10.mat'], 'half_edges2', 'NHE', 'HE_C', 'num_he', 'he_count', 'ghost_he', '-v7.3');        
        save(['half_edge_data_vert_10.mat'],'verts2',  '-v7.3');        
    end
    
end

%% plot vertices to check if all is good
if(1 == 0)
    figure('units','normalized','outerposition',[0 0 0.66 1]);
    for t = 1:Ns
        cla
        hold on;
        
        % load vertices and half edges
        HE = double(half_edges2(t).HE);
        V = double(verts2(t).V);        
        HE_vertex = double(half_edges2(t).HE_vertex);
        visited = double(half_edges2(t).visited);        
        Vertices = double(verts2(t).Vertices);   
        v_pos = double(verts2(t).v_pos);   

        num_he = NHE(t); 
        
        % plot 1
        if(1 == 1)
            for he = 1:num_he
                if(visited(he) == 1)
                    v = V(he,:);
                    plot([v(1) v(3)], [v(2) v(4)], '-', 'color', 'r')
                end
            end
            title(t)
            axis([-10 611 -10 611]);
            
        end
        
        % plot 2
        if(1 == 1)
            for ve = 1:length(Vertices)
                plot(v_pos(ve, 1), v_pos(ve, 2), '.', 'color', 'b');
%                 e = Vertices
%                 if(visited(he) == 1)
%                     v = V(he,:);
%                     plot([v(1) v(3)], [v(2) v(4)], '.-', 'color', 'r')
%                 end
            end
            title(t)
            axis([-10 611 -10 611]);
            drawnow
        end
    end
end




toc