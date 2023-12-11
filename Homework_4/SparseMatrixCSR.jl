import Base: size, getindex
import LinearAlgebra: mul!

struct SparseMatrixCSR{T} <: AbstractMatrix{T}
    rows::Vector{Int64}
    cols::Vector{Int64}
    values::Vector{T}
    m::Int
    n::Int
    function SparseMatrixCSR{T}(rows::Vector{Int64}, cols::Vector{Int64}, values::Vector{T}, m::Int, n::Int) where T
        new{T}(rows, cols, values, m, n)
    end
end

function size(A::SparseMatrixCSR)
    (A.m, A.n)
end


function getindex(A::SparseMatrixCSR{T}, i::Int, j::Int) where T
end
function mul!(out::AbstractVector, A::SparseMatrixCSR{T}, x::AbstractVector{T}) where T
    fill!(out, zero(T))
    for k in 1:length(A.values)
        out[A.rows[k]] += A.values[k] * x[A.cols[k]]
    end
    return out
end

using LinearAlgebra

function euclidean_norm(v::AbstractVector)
    sum_sq = sum(x^2 for x in v)
    return sqrt(sum_sq)
end

function plot_solutions(initial_u::AbstractVector, final_u::AbstractVector)
    indices = 1:length(initial_u)  # Assuming initial_u and final_u have the same length
    plot(indices, initial_u, label="Initial u", xlabel="Index", ylabel="Value", title="Initial vs Final Solution")
    plot!(indices, final_u, label="Final u")  # 'plot!' adds to the existing plot
end

function plot_solutions(initial_u::AbstractVector, final_u::AbstractVector)
    indices = 1:length(initial_u)  # Assuming initial_u and final_u have the same length
    plot(indices, initial_u, label="Initial u", xlabel="Index", ylabel="Value", title="Initial vs Final Solution")
    plot!(indices, final_u, label="Final u")  # 'plot!' adds to the existing plot
    display(plot!())
end

function iterate(A::SparseMatrixCSR{T}, b::Vector{T}, u::Vector{T}, a::T; tol::T = 1e-3) where T
    Au = A * u
    if isnothing(Au)
        error("Matrix-vector multiplication returned Nothing")
    end
    r = b - Au
    iterations = 0
    while norm(r) >= tol
        u += a * r
        Au = A * u
        if isnothing(Au)
            error("Matrix-vector multiplication returned Nothing")
        end
        r = b - Au
        norm = [euclidean_norm(final_u[1:i]) for i in 1:length(final_u)]
        iterations += 1
    end
    return u, iterations, norm(r)
end

n = 100
h = 1.0 / (n + 1)
rows_tridiagonal = Int64[]
cols_tridiagonal = Int64[]
values_tridiagonal = Float64[]
for i = 1:n+1
    push!(rows_tridiagonal, i)
    push!(cols_tridiagonal, i)
    push!(values_tridiagonal, 2 + h^2)
    if i > 1
        push!(rows_tridiagonal, i)
        push!(cols_tridiagonal, i-1)
        push!(values_tridiagonal, -1)
        push!(rows_tridiagonal, i-1)
        push!(cols_tridiagonal, i)
        push!(values_tridiagonal, -1)
    end
end

A_tridiagonal = SparseMatrixCSR{Float64}(rows_tridiagonal, cols_tridiagonal, values_tridiagonal, n+1, n+1)
b_tridiagonal = [h^2 * i for i in 1:n+1]
u_tridiagonal = zeros(n+1)
final_u = zeros(n+1)
final_u, iterations, final_norm = iterate(A_tridiagonal, b_tridiagonal, u_tridiagonal, 1.0)
using Plots
function plot_norms(norms)
    plot(norms, title="Norm of r over Iterations", xlabel="Iteration", ylabel="Norm of r", legend=false)
end
plot_solutions(u_tridiagonal, final_u)
plot_norms(norms)


# Verification
n = 100
h = 1.0 / (n + 1)
A_example = SparseMatrixCSR{Float64}(rows_tridiagonal, cols_tridiagonal, values, n+1, n+1)  
b_example = [h^2 * i for i in 1:n+1]
u_example = zeros(n+1)
a_example = 1.0


@code_warntype iterate(A_example, b_example, u_example, a_example)
out_example = zeros(n+1)
x_example = rand(n+1)  
@code_warntype mul!(out_example, A_example, x_example)
@time iterate(A_tridiagonal, b_tridiagonal, u_tridiagonal, 1.0) 




