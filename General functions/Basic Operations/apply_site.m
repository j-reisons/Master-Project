function [mps_matrix] = apply_site(mps_matrix,mpo_matrix)
%apply_site Does the same thing as apply on single site, but only returns a site tensor, not entire mps
result = contract(mps_matrix,3,mpo_matrix,4);
R_size = size(result);
result = reshape(permute(result,[3 1 4 2 5]),[R_size(1)*R_size(3),R_size(2)*R_size(4),R_size(5)]);
mps_matrix = result;
end

