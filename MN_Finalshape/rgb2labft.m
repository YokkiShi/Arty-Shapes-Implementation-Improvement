function ft = rgb2labft(t)
    [r, c] = size(t);
    ft = zeros(r, c);
    index = find(t > 0.008856);
    ft(index) = t(index) .^ (1 / 3);
    index = setdiff((1 : r * c)', index);
    ft(index) = 7.787 .* t(index) + 16 / 116;
end
