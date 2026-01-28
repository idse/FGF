function [znew, bw] = local_gaussian_smoothing(x,y,z,bw)

xyz = [x(:),y(:),z(:)];

if ~exist('bw','var')
    [~,~,bw] = ksdensity(xyz(:,[1,2]));
    bw = sqrt(bw);
end

znew = NaN(size(z));

for zi = 1:size(xyz,1)
    mu = xyz(zi,[1,2]);
    weights = exp(-0.5*(((xyz(:,1) - mu(1))/bw(1)).^2 + ((xyz(:,2) - mu(2))/bw(2)).^2))/(2*pi*bw(1)*bw(2));
    weights = weights/sum(weights);
    znew(zi) = sum(weights.*xyz(:,3),'omitnan');
end

end