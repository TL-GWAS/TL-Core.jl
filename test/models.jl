using Test
using GenesInteraction
using MLJ

@testset "Test InteractionTransformer" begin
    X = (rs1234=[1, 2, 3], rs455=[4, 5, 6], rs4489=[7, 8, 9], rstoto=[1, 2, 3])
    t = GenesInteraction.InteractionTransformer(r"^rs[0-9]+")
    mach = machine(t, X)
    fit!(mach)
    Xt = transform(mach, X)

    @test Xt == (
        rs1234 = [1, 2, 3],
        rs455 = [4, 5, 6],
        rs4489 = [7, 8, 9],
        rstoto = [1, 2, 3],
        rs1234_rs455 = [4.0, 10.0, 18.0],
        rs1234_rs4489 = [7.0, 16.0, 27.0],
        rs455_rs4489 = [28.0, 40.0, 54.0]
    )
    @test mach.fitresult.ninter == 3
    @test mach.fitresult.interaction_pairs == [:rs1234 => :rs455, :rs1234 => :rs4489, :rs455 => :rs4489]

end

@testset "InteractionLM pipeline" begin
    model = GenesInteraction.InteractionLM()
end