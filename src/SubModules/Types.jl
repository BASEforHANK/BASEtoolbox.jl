module Types

abstract type AbstractMacroModel end
struct OneAsset <: AbstractMacroModel end
struct TwoAsset <: AbstractMacroModel end
struct CompleteMarkets <: AbstractMacroModel end
export AbstractMacroModel, TwoAsset, OneAsset, CompleteMarkets

end
