function [HE, V, C1, v_pos, Vertices, HE_vertex, good_bond] = remove_half_edge_1(he, heo, v1, v2, HE, V, C1, t, Ntot, v_pos, Vertices, HE_vertex, good_bond, c1, c2)

% get position
p1 = v_pos(v1, :);
p2 = v_pos(v2, :);
p = ([mean([p1(1), p2(1)]) , mean([p1(2), p2(2)])]);

% get bonds
e1 = Vertices(v1, :);
e1(e1 == 0) = [];
e1(e1 == he) = [];
e2 = Vertices(v2, :);
e2(e2 == 0) = [];
e2(e2 == heo) = [];

% bond list for v_new
e = [e1, e2];
e(e == 0) = [];
if(length(e) > 8)
    he;
end

% if the new vertex has less than or equal 8 edges then change the data
if(length(e) <= 8)
    
    % new vertex given index of second vertex of pair and input position
    v_new = v2;
    v_pos(v1, :) = [0, 0];
    v_pos(v2, :) = [0, 0];
    v_pos(v_new, :) = p;
    Vertices(v1, :) = zeros(1, 8);
    Vertices(v2, :) = zeros(1, 8);
    Vertices(v_new, :) = zeros(1, 8);
    Vertices(v_new, 1:length(e)) = e;
    
    % give appropriate bonds the new vertex
    HE_vertex(he) = 0;
    HE_vertex(heo) = 0;
    for i = e1
        HE_vertex(i) = v_new;
    end
    for i = e2
        HE_vertex(i) = v_new;
    end
    
    % change next and previous for next and previous of he and heo
    HE(HE(he, 4), 5) = HE(he, 5);
    HE(HE(he, 5), 4) = HE(he, 4);
    HE(HE(heo, 4), 5) = HE(heo, 5);
    HE(HE(heo, 5), 4) = HE(heo, 4);
    
    % change endpoints of edge
    for i = e1
        V(i, 1:2) = p;
        if(HE(i, 6) > 0)
            V(HE(i, 6), 3:4)= p;
        end
    end
    for i = e2
        V(i, 1:2) = p;
        if(HE(i, 6) > 0)
            V(HE(i, 6), 3:4)= p;
        end
    end
    
    % remove he and heo
    HE(he, :) = zeros(1, 8);
    HE(heo, :) = zeros(1, 8);
    V(he, :) = zeros(1, 4);
    V(heo, :) = zeros(1, 4);
    
    % remove cell1 if now two sided
    if(c1 ~= Ntot + 1) && (C1(c1, t) > 0)
        he1 = C1(c1, t);
        e_count = 1;
        e = he1;
        
        he2 = HE(he1, 4);
        while(he2 ~= he1) && (he2 > 0)
            e_count = e_count + 1;
            e = [e, he2];
            he2 = HE(he2, 4);
        end
        
        if(e_count == 2)
            % get rid this mofo
            for he1 = e
                heo1 = HE(he1, 6);
                
                % get vertices
                v1 = HE_vertex(he1, 1);
                v2 = HE_vertex(heo1, 1);
                
                % remove edges from cell
                [HE, V, C1, v_pos, Vertices, HE_vertex] = remove_half_edge_1(he1, heo1, v1, v2, HE, V, C1, t, Ntot, v_pos, Vertices, HE_vertex, good_bond, c1, c2);
            end
            C1(c1, t) = 0;
        end
    end
    
    % remove cell2 if now two sided
    if(c2 ~= Ntot + 1) && (C1(c2, t) > 0)
        he1 = C1(c2, t);
        e_count = 1;
        e = he1;
        
        he2 = HE(he1, 4);
        while(he2 ~= he1) && (he2 > 0)
            e_count = e_count + 1;
            e = [e, he2];
            he2 = HE(he2, 4);
        end
        
        if(e_count == 2)
            % get rid this mofo
            for he1 = e
                heo1 = HE(he1, 6);
                
                % get vertices
                v1 = HE_vertex(he1, 1);
                v2 = HE_vertex(heo1, 1);
                
                % remove edges from cell
                [HE, V, C1, v_pos, Vertices, HE_vertex] = remove_half_edge_1(he1, heo1, v1, v2, HE, V, C1, t, Ntot, v_pos, Vertices, HE_vertex, good_bond, c1, c2);
            end
            C1(c2, t) = 0;
        end
    end
    
else
    good_bond(he) = 0;
end