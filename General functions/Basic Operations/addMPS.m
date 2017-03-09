function mps_out = addMPS(mps_1,mps_2)

N = length(mps_1);
D = size(mps_1{1},3);
mps_out = cell(1,N);

for k = 1:D
	mps_out{1}(1,:,k) = [mps_1{1}(1,:,k),mps_2{1}(1,:,k)];
end

for site = 2:(N-1)
	for k = 1:D
		mps_out{site}(:,:,k) = blkdiag(mps_1{site}(:,:,k),mps_2{site}(:,:,k));
	end
end

for k = 1:D
	mps_out{N}(:,1,k) = [mps_1{N}(:,1,k);mps_2{N}(:,1,k)];
end

end
