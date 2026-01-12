"""
Calibration utilities, labeled data containers, and SAM helpers for JCGE.
"""
module JCGECalibrate

using CSV
using DataFrames

export rho_from_sigma, sigma_from_rho
export calibrate_ces_share_scale
export LabeledVector, LabeledMatrix
export SAMTable, StartingValues, ModelParameters
export load_sam_table, compute_starting_values, compute_calibration_params
export load_canonical_sets, load_canonical_labels, load_canonical_subsets
export load_canonical_params, load_canonical_sam
export load_labeled_matrix, load_labeled_vector

"""
Convert CES elasticity sigma to rho (rho = (sigma - 1)/sigma).
"""
rho_from_sigma(sigma::Real) = (sigma - 1) / sigma

"""
Convert CES rho back to sigma (sigma = 1/(1-rho)).
"""
sigma_from_rho(rho::Real) = 1 / (1 - rho)

"""
Calibrate CES share/scale parameters (placeholder).
Return a NamedTuple so callers can evolve without breaking.
"""
function calibrate_ces_share_scale(; shares::AbstractVector, scale::Real=1.0)
    return (α = collect(shares), A = scale)
end

"""
Vector with labels and lookup indices.
"""
struct LabeledVector{T}
    data::Vector{T}
    labels::Vector{Symbol}
    index::Dict{Symbol,Int}
end

"""
Matrix with row/column labels and lookup indices.
"""
struct LabeledMatrix{T}
    data::Matrix{T}
    row_labels::Vector{Symbol}
    col_labels::Vector{Symbol}
    row_index::Dict{Symbol,Int}
    col_index::Dict{Symbol,Int}
end

"""
Construct a labeled vector from data and labels.
"""
LabeledVector(data::Vector{T}, labels::Vector{Symbol}) where {T} =
    LabeledVector{T}(data, labels, Dict(l => i for (i, l) in pairs(labels)))

"""
Construct a labeled matrix from data and row/column labels.
"""
LabeledMatrix(data::Matrix{T}, row_labels::Vector{Symbol}, col_labels::Vector{Symbol}) where {T} =
    LabeledMatrix{T}(
        data,
        row_labels,
        col_labels,
        Dict(l => i for (i, l) in pairs(row_labels)),
        Dict(l => i for (i, l) in pairs(col_labels)),
    )

"""
    load_labeled_matrix(path::AbstractString; label_col::String="label") -> LabeledMatrix

Load a CSV with a row label column and return a `LabeledMatrix`.
All data columns are coerced to `Float64`.
"""
function load_labeled_matrix(path::AbstractString; label_col::String="label")
    df = DataFrame(CSV.File(path))
    row_labels = Symbol.(df[:, label_col])
    col_labels = Symbol.(names(df)[names(df) .!= label_col])
    data = Matrix{Float64}(df[:, names(df) .!= label_col])
    return LabeledMatrix(data, row_labels, col_labels)
end

"""
    load_labeled_vector(path::AbstractString; label_col::String="label", value_col::String="value") -> LabeledVector

Load a CSV with a label and value column and return a `LabeledVector`.
Values are coerced to `Float64`.
"""
function load_labeled_vector(path::AbstractString; label_col::String="label", value_col::String="value")
    df = DataFrame(CSV.File(path))
    labels = Symbol.(df[:, label_col])
    values = Vector{Float64}(df[:, value_col])
    return LabeledVector(values, labels)
end

"""
Sum of labeled vector values.
"""
Base.sum(v::LabeledVector) = sum(v.data)

"""
Index into a labeled vector by symbol.
"""
Base.getindex(v::LabeledVector, label::Symbol) = v.data[v.index[label]]
"""
Index into a labeled vector by a list of symbols.
"""
Base.getindex(v::LabeledVector, labels::Vector{Symbol}) = v.data[[v.index[l] for l in labels]]

"""
Index into a labeled matrix by row and column symbols.
"""
Base.getindex(m::LabeledMatrix, r::Symbol, c::Symbol) = m.data[m.row_index[r], m.col_index[c]]
"""
Index into a labeled matrix by row/column symbol vectors.
"""
Base.getindex(m::LabeledMatrix, rows::Vector{Symbol}, cols::Vector{Symbol}) =
    m.data[[m.row_index[r] for r in rows], [m.col_index[c] for c in cols]]
"""
Index into a labeled matrix by a single row symbol and column vector.
"""
Base.getindex(m::LabeledMatrix, r::Symbol, cols::Vector{Symbol}) =
    vec(m.data[m.row_index[r], [m.col_index[c] for c in cols]])
"""
Index into a labeled matrix by row vector and single column symbol.
"""
Base.getindex(m::LabeledMatrix, rows::Vector{Symbol}, c::Symbol) =
    vec(m.data[[m.row_index[r] for r in rows], m.col_index[c]])

"""
Social Accounting Matrix with labeled accounts and key labels.
"""
struct SAMTable
    goods::Vector{Symbol}
    factors::Vector{Symbol}
    numeraire_factor_label::Symbol
    indirectTax_label::Symbol
    tariff_label::Symbol
    households_label::Symbol
    government_label::Symbol
    investment_label::Symbol
    restOfTheWorld_label::Symbol
    sam::LabeledMatrix{Float64}
end

"""
Starting values computed from a SAM table.
"""
struct StartingValues
    Td0::Float64
    Tz0::LabeledVector{Float64}
    Tm0::LabeledVector{Float64}
    F0::LabeledMatrix{Float64}
    Y0::LabeledVector{Float64}
    X0::LabeledMatrix{Float64}
    Z0::LabeledVector{Float64}
    M0::LabeledVector{Float64}
    tau_z::LabeledVector{Float64}
    tau_m::LabeledVector{Float64}
    Xp0::LabeledVector{Float64}
    FF::LabeledVector{Float64}
    Xg0::LabeledVector{Float64}
    Xv0::LabeledVector{Float64}
    E0::LabeledVector{Float64}
    Q0::LabeledVector{Float64}
    D0::LabeledVector{Float64}
    Sp0::Float64
    Sg0::Float64
    Sf::Float64
    pWe::LabeledVector{Float64}
    pWm::LabeledVector{Float64}
    pf0::LabeledVector{Float64}
    py0::LabeledVector{Float64}
    pz0::LabeledVector{Float64}
    pq0::LabeledVector{Float64}
    pe0::LabeledVector{Float64}
    pm0::LabeledVector{Float64}
    pd0::LabeledVector{Float64}
    epsilon0::Float64
end

"""
Calibrated model parameters derived from starting values.
"""
struct ModelParameters
    sigma::LabeledVector{Float64}
    psi::LabeledVector{Float64}
    eta::LabeledVector{Float64}
    phi::LabeledVector{Float64}
    alpha::LabeledVector{Float64}
    beta::LabeledMatrix{Float64}
    b::LabeledVector{Float64}
    ax::LabeledMatrix{Float64}
    ay::LabeledVector{Float64}
    mu::LabeledVector{Float64}
    lambda::LabeledVector{Float64}
    delta_m::LabeledVector{Float64}
    delta_d::LabeledVector{Float64}
    gamma::LabeledVector{Float64}
    xid::LabeledVector{Float64}
    xie::LabeledVector{Float64}
    theta::LabeledVector{Float64}
    ssp::Float64
    ssg::Float64
    tau_d::Float64
end

"""
Convert a vector of values to symbols.
"""
to_symbols(values::Vector) = Symbol.(values)

"""
    load_sam_table(file_path::AbstractString; kwargs...) -> SAMTable

Load a SAM CSV from `file_path` and return a `SAMTable`.
Missing values are replaced with zeros.
"""
function load_sam_table(file_path::AbstractString; goods::Vector{String} = ["BRD", "MLK"],
    factors::Vector{String} = ["CAP", "LAB"],
    numeraire_factor_label::String = "LAB",
    indirectTax_label::String = "IDT",
    tariff_label::String = "TRF",
    households_label::String = "HOH",
    government_label::String = "GOV",
    investment_label::String = "INV",
    restOfTheWorld_label::String = "EXT",
    label_col::String = "label")
    df = DataFrame(CSV.File(file_path))
    for col in eachcol(df)
        replace!(col, missing => 0)
    end
    label_name = label_col in names(df) ? label_col : ("Column1" in names(df) ? "Column1" : names(df)[1])
    row_labels = Symbol.(df[:, label_name])
    data_cols = names(df)[names(df) .!= label_name]
    col_labels = Symbol.(data_cols)
    sam = LabeledMatrix(Matrix{Float64}(df[:, data_cols]), row_labels, col_labels)
    return SAMTable(
        to_symbols(goods),
        to_symbols(factors),
        Symbol(numeraire_factor_label),
        Symbol(indirectTax_label),
        Symbol(tariff_label),
        Symbol(households_label),
        Symbol(government_label),
        Symbol(investment_label),
        Symbol(restOfTheWorld_label),
        sam,
    )
end

"""
    load_sam_table(io::IO; kwargs...) -> SAMTable

Load a SAM CSV from an `IO` stream and return a `SAMTable`.
Missing values are replaced with zeros.
"""
function load_sam_table(io::IO; goods::Vector{String} = ["BRD", "MLK"],
    factors::Vector{String} = ["CAP", "LAB"],
    numeraire_factor_label::String = "LAB",
    indirectTax_label::String = "IDT",
    tariff_label::String = "TRF",
    households_label::String = "HOH",
    government_label::String = "GOV",
    investment_label::String = "INV",
    restOfTheWorld_label::String = "EXT",
    label_col::String = "label")
    df = DataFrame(CSV.File(io))
    for col in eachcol(df)
        replace!(col, missing => 0)
    end
    label_name = label_col in names(df) ? label_col : ("Column1" in names(df) ? "Column1" : names(df)[1])
    row_labels = Symbol.(df[:, label_name])
    data_cols = names(df)[names(df) .!= label_name]
    col_labels = Symbol.(data_cols)
    sam = LabeledMatrix(Matrix{Float64}(df[:, data_cols]), row_labels, col_labels)
    return SAMTable(
        to_symbols(goods),
        to_symbols(factors),
        Symbol(numeraire_factor_label),
        Symbol(indirectTax_label),
        Symbol(tariff_label),
        Symbol(households_label),
        Symbol(government_label),
        Symbol(investment_label),
        Symbol(restOfTheWorld_label),
        sam,
    )
end

"""
    load_canonical_sets(dir::AbstractString) -> Dict{Symbol,Vector{Symbol}}

Load canonical sets from `sets.csv` (columns: set,item).
"""
function load_canonical_sets(dir::AbstractString)
    path = joinpath(dir, "sets.csv")
    df = DataFrame(CSV.File(path))
    sets = Dict{Symbol,Vector{Symbol}}()
    for row in eachrow(df)
        name = Symbol(row.set)
        push!(get!(sets, name, Symbol[]), Symbol(row.item))
    end
    return sets
end

"""
    load_canonical_labels(dir::AbstractString) -> Dict{Tuple{Symbol,Symbol},String}

Load canonical labels from `labels.csv` (columns: set,item,label).
"""
function load_canonical_labels(dir::AbstractString)
    path = joinpath(dir, "labels.csv")
    if !isfile(path)
        return Dict{Tuple{Symbol,Symbol},String}()
    end
    df = DataFrame(CSV.File(path))
    labels = Dict{Tuple{Symbol,Symbol},String}()
    for row in eachrow(df)
        labels[(Symbol(row.set), Symbol(row.item))] = String(row.label)
    end
    return labels
end

"""
    load_canonical_subsets(dir::AbstractString) -> Dict{Symbol,Vector{Symbol}}

Load canonical subsets from `subsets.csv` (columns: subset,parent_set,item).
"""
function load_canonical_subsets(dir::AbstractString)
    path = joinpath(dir, "subsets.csv")
    if !isfile(path)
        return Dict{Symbol,Vector{Symbol}}()
    end
    df = DataFrame(CSV.File(path))
    subsets = Dict{Symbol,Vector{Symbol}}()
    for row in eachrow(df)
        name = Symbol(row.subset)
        push!(get!(subsets, name, Symbol[]), Symbol(row.item))
    end
    return subsets
end

"""
    load_canonical_params(dir::AbstractString) -> DataFrame

Load canonical parameters from `params.csv` (columns: name,set1,index1,set2,index2,set3,index3,value,section).
"""
function load_canonical_params(dir::AbstractString)
    path = joinpath(dir, "params.csv")
    return DataFrame(CSV.File(path))
end

"""
    load_canonical_sam(dir::AbstractString; kwargs...) -> SAMTable

Load canonical SAM from `sam.csv` (label column: `label`).
"""
function load_canonical_sam(dir::AbstractString; kwargs...)
    return load_sam_table(joinpath(dir, "sam.csv"); kwargs..., label_col="label")
end

"""
    compute_starting_values(sam_table::SAMTable) -> StartingValues

Compute calibrated starting values from a `SAMTable`.

The output includes base flows, tax rates, and benchmark prices.
"""
function compute_starting_values(sam_table::SAMTable)
    sam = sam_table.sam
    goods = sam_table.goods
    factors = sam_table.factors
    Td0 = sam[sam_table.government_label, sam_table.households_label]
    Tz0 = sam[sam_table.indirectTax_label, goods]
    Tm0 = sam[sam_table.tariff_label, goods]
    F0 = sam[factors, goods]
    Y0 = vec(sum(F0, dims=1))
    X0 = sam[goods, goods]
    Z0 = vec(sum(X0, dims=1)) .+ Y0
    M0 = sam[sam_table.restOfTheWorld_label, goods]
    tau_z = Tz0 ./ Z0
    tau_m = Tm0 ./ M0
    Xp0 = sam[goods, sam_table.households_label]
    FF = sam[sam_table.households_label, factors]
    Xg0 = sam[goods, sam_table.government_label]
    Xv0 = sam[goods, sam_table.investment_label]
    E0 = sam[goods, sam_table.restOfTheWorld_label]
    Q0 = Xp0 .+ Xg0 .+ Xv0 .+ vec(sum(X0, dims=2))
    D0 = (1 .+ tau_z) .* Z0 .- E0
    Sp0 = sam[sam_table.investment_label, sam_table.households_label]
    Sg0 = sam[sam_table.investment_label, sam_table.government_label]
    Sf = sam[sam_table.investment_label, sam_table.restOfTheWorld_label]
    pWe = ones(length(goods))
    pWm = ones(length(goods))
    pf0 = ones(length(factors))
    py0 = ones(length(goods))
    pz0 = ones(length(goods))
    pq0 = ones(length(goods))
    pe0 = ones(length(goods))
    pm0 = ones(length(goods))
    pd0 = ones(length(goods))
    epsilon0 = 1.0
    return StartingValues(
        Td0,
        LabeledVector(Tz0, goods),
        LabeledVector(Tm0, goods),
        LabeledMatrix(F0, factors, goods),
        LabeledVector(Y0, goods),
        LabeledMatrix(X0, goods, goods),
        LabeledVector(Z0, goods),
        LabeledVector(M0, goods),
        LabeledVector(tau_z, goods),
        LabeledVector(tau_m, goods),
        LabeledVector(Xp0, goods),
        LabeledVector(FF, factors),
        LabeledVector(Xg0, goods),
        LabeledVector(Xv0, goods),
        LabeledVector(E0, goods),
        LabeledVector(Q0, goods),
        LabeledVector(D0, goods),
        Sp0,
        Sg0,
        Sf,
        LabeledVector(pWe, goods),
        LabeledVector(pWm, goods),
        LabeledVector(pf0, factors),
        LabeledVector(py0, goods),
        LabeledVector(pz0, goods),
        LabeledVector(pq0, goods),
        LabeledVector(pe0, goods),
        LabeledVector(pm0, goods),
        LabeledVector(pd0, goods),
        epsilon0,
    )
end

"""
    compute_calibration_params(sam_table::SAMTable, start::StartingValues) -> ModelParameters

Compute calibrated model parameters from a `SAMTable` and starting values.

The parameters follow standard CGE calibration formulas for CES/CET nests.
"""
function compute_calibration_params(sam_table::SAMTable, start::StartingValues)
    goods = sam_table.goods
    factors = sam_table.factors
    sigma = fill(2.0, length(goods))
    psi = fill(2.0, length(goods))
    eta = (sigma .- 1.0) ./ sigma
    phi = (psi .+ 1.0) ./ psi
    alpha = start.Xp0.data ./ sum(start.Xp0.data)
    beta = start.F0.data ./ sum(start.F0.data, dims=1)
    b = start.Y0.data ./ vec(prod(start.F0.data .^ beta, dims=1))
    ax = start.X0.data ./ transpose(start.Z0.data)
    ay = start.Y0.data ./ start.Z0.data
    mu = start.Xg0.data ./ sum(start.Xg0.data)
    lambda = start.Xv0.data ./ (start.Sp0 + start.Sg0 + start.Sf)
    delta_m = (1 .+ start.tau_m.data) .* start.M0.data .^ (1 .- eta) ./
              ((1 .+ start.tau_m.data) .* start.M0.data .^ (1 .- eta) .+ start.D0.data .^ (1 .- eta))
    delta_d = start.D0.data .^ (1 .- eta) ./
              ((1 .+ start.tau_m.data) .* start.M0.data .^ (1 .- eta) .+ start.D0.data .^ (1 .- eta))
    gamma = start.Q0.data ./ (delta_m .* start.M0.data .^ eta .+ delta_d .* start.D0.data .^ eta) .^ (1 ./ eta)
    xie = start.E0.data .^ (1 .- phi) ./ (start.E0.data .^ (1 .- phi) .+ start.D0.data .^ (1 .- phi))
    xid = start.D0.data .^ (1 .- phi) ./ (start.E0.data .^ (1 .- phi) .+ start.D0.data .^ (1 .- phi))
    theta = start.Z0.data ./ (xie .* start.E0.data .^ phi .+ xid .* start.D0.data .^ phi) .^ (1 ./ phi)
    ssp = start.Sp0 / sum(start.FF.data)
    ssg = start.Sg0 / (start.Td0 + sum(start.Tz0.data) + sum(start.Tm0.data))
    tau_d = start.Td0 / sum(start.FF.data)
    return ModelParameters(
        LabeledVector(sigma, goods),
        LabeledVector(psi, goods),
        LabeledVector(eta, goods),
        LabeledVector(phi, goods),
        LabeledVector(alpha, goods),
        LabeledMatrix(beta, factors, goods),
        LabeledVector(b, goods),
        LabeledMatrix(ax, goods, goods),
        LabeledVector(ay, goods),
        LabeledVector(mu, goods),
        LabeledVector(lambda, goods),
        LabeledVector(delta_m, goods),
        LabeledVector(delta_d, goods),
        LabeledVector(gamma, goods),
        LabeledVector(xid, goods),
        LabeledVector(xie, goods),
        LabeledVector(theta, goods),
        ssp,
        ssg,
        tau_d,
    )
end

end # module
