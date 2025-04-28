"""
    find_field_with_value(struct_instance, field_value)

Find the names of fields in a struct instance that contain a specific value.

Parameters:

  - `struct_instance`: An instance of a Julia struct.
  - `field_value`: The value to search for within the fields of the struct.

Returns:

  - A list of field names that contain the specified value, or an empty list if no such
    field is found.
"""
function find_field_with_value(struct_instance, field_value)
    matching_fields = Symbol[]

    for field_name in fieldnames(typeof(struct_instance))
        field_val = getfield(struct_instance, field_name)
        if field_value == field_val || field_value in field_val
            push!(matching_fields, field_name)
        end
    end

    if isempty(matching_fields)
        @warn "No field with value $field_value found!"
    end

    return matching_fields
end

"""
    find_field_with_value(struct_instance, field_value, ss)

Wrapper function for `find_field_with_value` that allows for selection based on the
presence of "SS" in the field name.
"""
function find_field_with_value(struct_instance, field_value, ss)
    matching_fields = find_field_with_value(struct_instance, field_value)
    if length(matching_fields) == 0
        @warn "No field with value $field_value found!"
    elseif length(matching_fields) == 1
        return matching_fields[1]
    elseif length(matching_fields) == 2
        if ss == true
            for field in matching_fields
                if occursin("SS", string(field))
                    return field
                end
            end
        else
            for field in matching_fields
                if !occursin("SS", string(field))
                    return field
                end
            end
        end
    elseif length(matching_fields) > 2
        @error "More than two fields with value $field_value found!"
    end
end

"""
    mapround(b; digits = 4)

The `mapround` function rounds every element of a matrix `b` to `digits` number of decimal
places.

# Arguments

  - `b`: a matrix of numbers to be rounded
  - `digits`: number of decimal places to round to (default is 4)

# Returns

A matrix with the same shape as `b` where every element has been rounded to `digits` number
of decimal places.
"""
function mapround(b; digits = 4)
    map(x -> round.(x; digits = digits), b)
end
