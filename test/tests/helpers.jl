function is_named_struct(x)
    t = typeof(x)
    return isstructtype(t) && !isa(x, AbstractArray) && !isempty(fieldnames(t))
end

function intersection(test, current)
    test_names = propertynames(test)
    current_names = propertynames(current)
    intersection_names = intersect(test_names, current_names)
    missing_in_current = setdiff(test_names, current_names)
    missing_in_test = setdiff(current_names, test_names)

    # Warn about missing fields
    if !isempty(missing_in_current)
        @warn "Fields missing in current:" missing_in_current
    end
    if !isempty(missing_in_test)
        @warn "Fields missing in test:" missing_in_test
    end

    return intersection_names
end

function wrapper_test(test, current; tol = 1e-6)
    if all(typeof.(test) .== Symbol) || all(typeof.(test) .== String)
        return @test test == current
    else
        return @test all(isapprox.(test, current; atol = tol))
    end
end
