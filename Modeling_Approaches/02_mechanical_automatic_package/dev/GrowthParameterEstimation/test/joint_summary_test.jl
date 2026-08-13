using Test
using DataFrames
using GrowthParameterEstimation

@testset "Joint comparison summaries" begin
    joint = DataFrame(
        model = repeat(["A", "B"], 2),
        cell_line = repeat(["naive", "resistant"], inner = 2),
        density = fill("20k", 4),
        bic = [10.0, 13.0, 20.0, 19.0],
    )
    summary = summarize_joint_bic(joint)
    @test nrow(summary.aggregate) == 2
    @test summary.aggregate.model[1] == "A"
    @test summary.aggregate.total_bic[1] == 30.0
    @test summary.aggregate.delta_total_bic[1] == 0.0
    @test nrow(summary.winners) == 2

    parameters = DataFrame(
        model = ["A", "A", "A", "A"],
        parameter = ["kill", "kill", "delay", "delay"],
        value = [0.20, 0.22, 1.0, 4.0],
        lower_bound = [0.0, 0.0, 0.0, 0.0],
        upper_bound = [1.0, 1.0, 5.0, 5.0],
    )
    stability = summarize_joint_parameter_stability(parameters)
    kill = only(filter(row -> row.parameter == "kill", stability))
    delay = only(filter(row -> row.parameter == "delay", stability))
    @test kill.stability_class == "tight"
    @test delay.stability_class == "broad"
    @test isapprox(kill.bound_range_fraction, 0.02)

    grouped = DataFrame(
        model = repeat(["A", "B", "C"], 4),
        cell_line = repeat(["naive", "naive", "cis", "cis"], inner = 3),
        density = repeat(["20k", "30k", "20k", "30k"], inner = 3),
        bic = [10.0, 12.0, 18.0, 11.0, 15.0, 19.0, 14.0, 13.0, 20.0, 16.0, 15.0, 21.0],
    )
    top_by_cell = summarize_joint_bic_by_group(grouped; top_n = 2)
    @test nrow(top_by_cell) == 4
    @test all(combine(groupby(top_by_cell, :cell_line), nrow => :count).count .== 2)
    @test top_by_cell.model[top_by_cell.cell_line .== "naive"] == ["A", "B"]
    @test top_by_cell.model[top_by_cell.cell_line .== "cis"] == ["B", "A"]
    @test all(first(group).delta_total_bic == 0.0 for group in groupby(top_by_cell, :cell_line))
end
