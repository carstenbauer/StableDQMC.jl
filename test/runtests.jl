using StableDQMC
using Test, Random, LinearAlgebra




@testset "StableDQMC.jl" begin

    @testset "UDT Type" begin
        Random.seed!(1234)
        x = rand(10,10)

        @test typeof(udt(x)) == UDT{Float64,Float64,Array{Float64,2}}
        @test typeof(udt!(copy(x))) == UDT{Float64,Float64,Array{Float64,2}}

        F = udt(x)
        u,d,t = F

        @test isapprox(u * u', I) # unitarity of u
        @test isreal(d)
        @test sort!(count.(iszero,eachcol(t))) == 0:9 # triangularity of t (istriu(t) == false due to pivoting)
        @test isapprox(Matrix(F), x)
        @test isapprox(u * Diagonal(d) * t, x)

        @test size(F) == size(x)
        @test typeof(similar(F)) == UDT{Float64,Float64,Array{Float64,2}}

        N = udt(3)
        @test N.U == reshape([1.0], 1,1)
        @test N.D == [3.0]
        @test N.T == reshape([1.0], 1,1)


        # operations
        @test isapprox(inv(F), inv(x))
        G = fact_mult(F, F)
        @test isapprox(Matrix(G), x*x)
    end


    @testset "SVD decomposition variants" begin
        for dtype in (Float64, ComplexF64)
            for decomp in (gesdd, gesvd, gesvj, genericsvd)
                x = rand(dtype, 5,5)

                F = decomp(x)
                @test typeof(F) <: SVD
                @test isapprox(Matrix(F), x)
                @test isapprox(F.U * F.U', I)
                @test isreal(F.S)
                @test isapprox(F.Vt * F.Vt', I)
            end
        end

        a, b = rand(5,5), rand(5,5)
        A, B = svd(a), svd(b)
        @test isapprox(Matrix(fact_mult(A,B)), a*b)
    end









    @testset "QR / UDT operations" begin
        for dtype in (Float64, ComplexF64)
            x = rand(dtype, 5,5)
            X = udt(x)
            y = rand(dtype, 5,5)
            Y = udt(y)
            res = similar(x)

            inv_I_plus_x = inv(I + x)
            inv_I_plus_xydagger = inv(I + x * y')
            inv_sum_xy = inv(x + y)

            @test isapprox(inv_one_plus(X), inv_I_plus_x)
            @test isapprox(inv_one_plus!(res, X), inv_I_plus_x)
            @test isapprox(res, inv_I_plus_x)

            @test typeof(udt_inv_one_plus(X, Y)) <: UDT
            @test isapprox(Matrix(udt_inv_one_plus(X, Y)), inv_I_plus_xydagger)
            @test isapprox(inv_one_plus(X,Y), inv_I_plus_xydagger)
            @test isapprox(inv_one_plus!(res, X,Y), inv_I_plus_xydagger)
            @test isapprox(res, inv(I + x * y'))

            @test typeof(udt_inv_sum(X,Y)) <: UDT
            @test isapprox(Matrix(udt_inv_sum(X,Y)), inv_sum_xy)
            @test isapprox(inv_sum(X,Y), inv_sum_xy)
            @test isapprox(inv_sum!(res, X,Y), inv_sum_xy)
            @test isapprox(res, inv_sum_xy)

            # Loh et al variants
            @test typeof(udt_inv_one_plus_loh(X)) <: UDT
            @test isapprox(Matrix(udt_inv_one_plus_loh(X)), inv_I_plus_x)
            @test isapprox(inv_one_plus_loh(X), inv_I_plus_x)
            @test isapprox(inv_one_plus_loh!(res,X), inv_I_plus_x)
            @test isapprox(res, inv_I_plus_x)

            @test typeof(udt_inv_sum_loh(X,Y)) <: UDT
            @test isapprox(Matrix(udt_inv_sum_loh(X,Y)), inv_sum_xy)
            @test isapprox(inv_sum_loh(X,Y), inv_sum_xy)
            @test isapprox(inv_sum_loh!(res,X,Y), inv_sum_xy)
            @test isapprox(res, inv_sum_xy)
        end
    end



    @testset "SVD / USVt operations" begin
        for dtype in (Float64, ComplexF64)
            x = rand(dtype, 5,5)
            X = svd(x)
            y = rand(dtype, 5,5)
            Y = svd(y)
            res = similar(x)

            inv_I_plus_x = inv(I + x)
            inv_sum_xy = inv(x + y)

            @test isapprox(inv_one_plus(X), inv_I_plus_x)
            @test isapprox(Matrix(svd_inv_one_plus_loh(X)), inv_I_plus_x)
            @test isapprox(inv_one_plus_loh(X), inv_I_plus_x)
            @test isapprox(inv_one_plus_loh!(res, X), inv_I_plus_x)
            @test isapprox(res, inv_I_plus_x)

            @test isapprox(Matrix(svd_inv_sum_loh(X, Y)), inv_sum_xy)
            @test isapprox(inv_sum_loh(X, Y), inv_sum_xy)
            @test isapprox(inv_sum_loh!(res, X, Y), inv_sum_xy)
            @test isapprox(res, inv_sum_xy)
        end
    end

end
