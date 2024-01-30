function centraldiff(x)
    # calculate central difference of a vector
    if ndims(x) != 1
        error("wrong number of dimensions, only 1 allowed")
    else
        a = similar(x)
        dx = diff(x) ./ 2.0
        a .= [dx; dx[end]]
        a .+= [dx[1]; dx]
        return a
    end

end

function centralderiv(y, x, dims)
    # calculates the appriximate derivative of array y on 
    # mesh x along dimension dims. The derivative is based on the central difference
    dy = mapslices(centraldiff, y, dims = dims)
    dx = mapslices(centraldiff, x, dims = dims)
    return dy ./ dx
end
