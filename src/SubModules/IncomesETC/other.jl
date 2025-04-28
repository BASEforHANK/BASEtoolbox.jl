function mutil(c::Union{Real,AbstractArray}, m_par)
    if m_par.ξ == 1.0
        mutil = 1.0 ./ c
    elseif m_par.ξ == 2.0
        mutil = 1.0 ./ (c .^ 2)
    elseif m_par.ξ == 4.0
        mutil = 1.0 ./ ((c .^ 2) .^ 2)
    else
        mutil = c .^ (-m_par.ξ)
    end
    return mutil
end

function mutil!(mu::AbstractArray, c::AbstractArray, m_par)
    if m_par.ξ == 1.0
        mu .= 1.0 ./ c
    elseif m_par.ξ == 2.0
        mu .= 1.0 ./ (c .^ 2)
    elseif m_par.ξ == 4.0
        mu .= 1.0 ./ ((c .^ 2) .^ 2)
    else
        mu .= c .^ (-m_par.ξ)
    end
    return mu
end

function invmutil(mu, m_par)
    if m_par.ξ == 1.0
        c = 1.0 ./ mu
    elseif m_par.ξ == 2.0
        c = 1.0 ./ (sqrt.(mu))
    elseif m_par.ξ == 4.0
        c = 1.0 ./ (sqrt.(sqrt.(mu)))
    else
        c = 1.0 ./ mu .^ (1.0 ./ m_par.ξ)
    end
    return c
end

function invmutil!(c, mu, m_par)
    if m_par.ξ == 1.0
        c .= 1.0 ./ mu
    elseif m_par.ξ == 2.0
        c .= 1.0 ./ (sqrt.(mu))
    elseif m_par.ξ == 4.0
        c .= 1.0 ./ (sqrt.(sqrt.(mu)))
    else
        c .= 1.0 ./ mu .^ (1.0 ./ m_par.ξ)
    end
    return c
end
