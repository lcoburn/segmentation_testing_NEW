function[neighs, vertices] = remove_neigh(neighs, vertices, ind)


neighs(ind) = [];
vertices(:, ind) = [];