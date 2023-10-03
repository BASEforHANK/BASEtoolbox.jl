# contains
# - beta_pars
# - ig_pars



@doc raw"""
    beta_pars(betamean,betavariance)

Compute the location and shape parameter of the Beta distribution from the mean and variance.

# Arguments
- `betamean`: prior mean of beta distribution [scalar]
- `betavariance`: prior variance of beta distribution [scalar]

# Returns
- `a`: location parameter of beta distribution [scalar]
- `b`: shape parameter of beta distribution [scalar]
"""
function beta_pars(betamean, betavariance)
    b =
        (betamean - 2 * betamean^2 + betamean^3 - betavariance + betamean * betavariance) /
        betavariance
    a = -betamean * b / (betamean - 1)

    return a, b
end

@doc raw"""
    gamma_pars(gammamean,gammavariance)

Compute the location and shape parameter of the Gamma distribution from the mean and variance.

# Arguments
- `gammamean`: prior mean of gamma distribution [scalar]
- `gammavariance`: prior variance of gamma distribution [scalar]

# Returns
- `a`: location parameter of gamma distribution [scalar]
- `b`: shape parameter of gamma distribution [scalar]
"""
function gamma_pars(gammamean, gammavariance)
    a = gammamean^2 / gammavariance
    b = gammavariance / gammamean

    return a, b
end

@doc raw"""
    ig_pars(igmean,igvariance)

Compute the location and shape parameter of the Inverse Gamma distribution from the mean and variance.

# Arguments
- `igmean`: prior mean of inverse gamma distribution [scalar]
- `igvariance`: prior variance of inverse gamma distribution [scalar]

Ouputs:
- `a`: location parameter of inverse gamma distribution [scalar]
- `b`: shape parameter of inverse gamma distribution [scalar]
"""
function ig_pars(igmean, igvariance)
    a = igmean^2 / igvariance + 2
    b = igmean * (igmean^2 / igvariance + 1)

    return a, b
end
