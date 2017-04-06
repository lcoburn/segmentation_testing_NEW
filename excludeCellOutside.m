function [Cells, count_t] = excludeCellOutside(param, bounds, t, Cells, count_t)

    % Exclude cells from the boundaries. Return Cells and count_t.
    % Cells: array [timePoints x Segments] containing the IDs of alive cells per
    % time point
    % count_t: number of cells alive per timePoint inside the boundaries
    
    bx = bounds(:, 1);
    by = bounds(:, 2);

    aliveIDs=cellfun(@(x) ismember(t,x),{param.tracks.t});

    for id = aliveIDs
        tIdx = param.tracks(1).t==t;
        xc = double(param.tracks(id).cent(tIdx, 1));
        yc = double(param.tracks(id).cent(tIdx, 2));
        in = inpolygon(xc, yc, bx, by);
        
        if(in == 1)
            count_t(t) = count_t(t) + 1;
            Cells(t, count_t(t)) = id;
        end
    end

end