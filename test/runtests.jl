using StableDQMC
using Test, Random, LinearAlgebra



@testset "StableDQMC.jl" begin

    @testset "Opt-in JacobiSVD/GenericSVD" begin
        @test_throws UndefVarError gesvj
        @test_throws UndefVarError gesvj!
        @test_throws UndefVarError genericsvd
        @test_throws UndefVarError genericsvd!
        StableDQMC.addJacobiSVD()
        StableDQMC.addGenericSVD()
    end

    @testset "UDT Type" begin
        Random.seed!(1234)
        x = rand(10,10)

        @test typeof(udt(x)) == UDT{Float64,Float64,Array{Float64,2}}
        @test typeof(udt!(copy(x))) == UDT{Float64,Float64,Array{Float64,2}}

        F = udt(x)
        u,d,t = F

        @test u * u' ≈ I # unitarity of u
        @test isreal(d)
        @test sort!(count.(iszero,eachcol(t))) == 0:9 # triangularity of t (istriu(t) == false due to pivoting)
        @test Matrix(F) ≈ x
        @test u * Diagonal(d) * t ≈ x

        @test size(F) == size(x)
        @test typeof(similar(F)) == UDT{Float64,Float64,Array{Float64,2}}

        N = udt(3)
        @test N.U == reshape([1.0], 1,1)
        @test N.D == [3.0]
        @test N.T == reshape([1.0], 1,1)


        # operations
        @test inv(F) ≈ inv(x)
        G = fact_mult(F, F)
        @test Matrix(G) ≈ x*x
    end


    @testset "SVD decomposition variants" begin
        for dtype in (Float64, ComplexF64)
            for decomp in (gesdd, gesvd, gesvj, genericsvd)
                x = rand(dtype, 5,5)

                F = decomp(x)
                @test typeof(F) <: SVD
                @test Matrix(F) ≈ x
                @test F.U * F.U' ≈ I
                @test isreal(F.S)
                @test F.Vt * F.Vt' ≈ I
            end
        end

        a, b = rand(5,5), rand(5,5)
        A, B = svd(a), svd(b)
        @test Matrix(fact_mult(A,B)) ≈ a*b
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

            @test inv_one_plus(X) ≈ inv_I_plus_x
            @test inv_one_plus!(res, X) ≈ inv_I_plus_x
            @test res ≈ inv_I_plus_x

            @test typeof(udt_inv_one_plus(X, Y)) <: UDT
            @test Matrix(udt_inv_one_plus(X, Y)) ≈ inv_I_plus_xydagger
            @test inv_one_plus(X,Y) ≈ inv_I_plus_xydagger
            @test inv_one_plus!(res, X,Y) ≈ inv_I_plus_xydagger
            @test res ≈ inv(I + x * y')

            @test typeof(udt_inv_sum(X,Y)) <: UDT
            @test Matrix(udt_inv_sum(X,Y)) ≈ inv_sum_xy
            @test inv_sum(X,Y) ≈ inv_sum_xy
            @test inv_sum!(res, X,Y) ≈ inv_sum_xy
            @test res ≈ inv_sum_xy

            # Loh et al variants
            @test typeof(udt_inv_one_plus_loh(X)) <: UDT
            @test Matrix(udt_inv_one_plus_loh(X)) ≈ inv_I_plus_x
            @test inv_one_plus_loh(X) ≈ inv_I_plus_x
            @test inv_one_plus_loh!(res,X) ≈ inv_I_plus_x
            @test res ≈ inv_I_plus_x

            @test typeof(udt_inv_sum_loh(X,Y)) <: UDT
            @test Matrix(udt_inv_sum_loh(X,Y)) ≈ inv_sum_xy
            @test inv_sum_loh(X,Y) ≈ inv_sum_xy
            @test inv_sum_loh!(res,X,Y) ≈ inv_sum_xy
            @test res ≈ inv_sum_xy
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

            @test inv_one_plus(X) ≈ inv_I_plus_x
            @test Matrix(svd_inv_one_plus_loh(X)) ≈ inv_I_plus_x
            @test inv_one_plus_loh(X) ≈ inv_I_plus_x
            @test inv_one_plus_loh!(res, X) ≈ inv_I_plus_x
            @test res ≈ inv_I_plus_x

            @test Matrix(svd_inv_sum(X, Y)) ≈ inv_sum_xy
            @test inv_sum(X, Y) ≈ inv_sum_xy
            @test inv_sum!(res, X, Y) ≈ inv_sum_xy
            @test res ≈ inv_sum_xy

            @test Matrix(svd_inv_sum_loh(X, Y)) ≈ inv_sum_xy
            @test inv_sum_loh(X, Y) ≈ inv_sum_xy
            @test inv_sum_loh!(res, X, Y) ≈ inv_sum_xy
            @test res ≈ inv_sum_xy
        end
    end



    @testset "Helpers" begin
        x = [1 2; 3 4]
        @test cond(svdvals(x)) ≈ cond(x)

        R, svs = calc_Bchain(x, 5)
        @test R ≈ x^5
        for i in 1:5
            @test svs[i,:] ≈ log.(svdvals(x^i))
        end


        # SVD / UDV
        R, svs = calc_Bchain_svd(x, 5)
        xsvd = svd(x^5)
        @test R.U ≈ xsvd.U
        @test R.S ≈ xsvd.S
        @test R.Vt ≈ xsvd.Vt
        for i in 1:5
            @test svs[i,:] ≈ log.(svdvals(x^i))
        end


        # QR / UDT
        R, svs = calc_Bchain_qr(x, 5)
        xudt = udt(x^5)
        @test R.U ≈ xudt.U
        @test R.D ≈ xudt.D
        @test R.T ≈ xudt.T
        for i in 1:5
            @test svs[i,:] ≈ log.(udt(x^i).D)
        end


        g = calc_tdgf_qr(x,2,3)
        @test g ≈ inv(x^-2 + x^3) # [B^-N1 + B^N2]^-1
        g = calc_tdgf_svd(x,2,3)
        @test g ≈ inv(x^-2 + x^3) # [B^-N1 + B^N2]^-1
    end

end
