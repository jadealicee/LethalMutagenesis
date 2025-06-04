using Markdown
using InteractiveUtils
using Random, Distributions, StatsBase, Parameters, DataFrames, Plots
using CSV, Tables
using Pkg

#@with_kw struct Params

# physical constants
const kt_boltzmann = 0.001987204118 * (273.15+37)
const ΔΔG = Normal(1.0, 1.7)

# fixed parameters
const G = 10 # number_of_genes
const sim_length = 500 # number of generations 500
const F = -5.0 # initial_free_energy (of all proteins)

# variable parameters
const U = Poisson(1.0) # mutation_rate - 4.5
const L = 0.0 # lethal_fraction - 0.1
const N = 500 # start_popsize
const K = 500 # carrying_capacity
const R = 1 # fecundity - 4
const Rep = 2 #Replicate No. #We don't really need this if we are incorporating seeds into file names. What do you think?
#description::String # manual entry
#replicate::Char # manual entry
const start_fitness = prod([1 / (1 + ℯ^(ΔG/kt_boltzmann)) for ΔG in fill(F, G)])

mutable struct Virus
    μ_counts::Vector{Int64}
    ΔG_list::Vector{Float64}
    fitness::Float64
    zombie::Bool
end

function update_fitness!(virus::Virus) #, p::Params
    #@unpack kt_boltzmann = p
    # I could check if any ΔG value > 1
    virus.fitness = prod([1 / (1 + ℯ^(ΔG/kt_boltzmann)) for ΔG in virus.ΔG_list])
end
   #begin
#p = Params(description = "preparatory", replicate = '1')
#end
function mutate!(virus::Virus) #p::Params
    #@unpack G, L, U, ΔΔG = p
    number_of_mutations = only(rand(U, 1))
    ΔΔG_values = rand(ΔΔG, number_of_mutations)
    mutgene_coord = rand((1:G), number_of_mutations)
    for (index, gene_id) in enumerate(mutgene_coord)
        virus.μ_counts[gene_id] += 1
        virus.ΔG_list[gene_id] = virus.ΔG_list[gene_id] + ΔΔG_values[index]
        if only(rand(Float64, 1)) < L
            virus.fitness = 0
        end
    end
    ( virus.fitness > 0  ) && ( update_fitness!(virus) )
    if virus.fitness <= 0 # make me efficient!
        virus.zombie = true
    end
    return virus.fitness
end

function reproduce(parent::Virus) #, p::Params
    sprog = deepcopy(parent)
    mutate!(sprog)
    return sprog
end

function initialize_population()
    #@unpack F, G, N = p
    initial_population = [Virus(zeros(Int, G), fill(F, G), start_fitness, false) for _ in 1:N]
    return initial_population
end
    
function get_weights(populace)
    weights = [v.fitness for v in populace]
    return Weights(weights/sum(weights))
end
    
function probabilistic_round(number)
    frac = abs(number - floor(number))
    #frac = number % floor(number)
    if only(rand(Float64, 1)) < frac
        return ceil(Int, number)
    end
    return floor(Int, number)
end
    
function initialize_report()
    report = DataFrame(psiz = Int[], q1fit = Float64[], meanfit = Float64[],
        q2fit = Float64[], maxfit = Float64[], minfree = Float64[], 
        meanfree = Float64[], maxfree = Float64[], minmut = Float64[], 
        meanmut = Float64[], maxmut = Float64[])
    return report
end
    
function report_update!(populace, report)
    push!(report, 
        [length(populace), #psiz
        quantile([v.fitness for v in populace], 0.25), #q1fit
        mean([v.fitness for v in populace]), #meanfit
        median([v.fitness for v in populace]), #q2fit
        maximum([v.fitness for v in populace]), #maxfit
        mean([minimum(v.ΔG_list) for v in populace]), #minfree
        mean([mean(v.ΔG_list) for v in populace]), #meanfree
        mean([maximum(v.ΔG_list) for v in populace]), #maxfree
        minimum([sum(v.μ_counts) for v in populace]), #minmut
        mean([sum(v.μ_counts) for v in populace]), #meanmut
        maximum([sum(v.μ_counts) for v in populace]), #maxmut
        ])
end
    
function plot_simulation(report)
    abscissa = range(1, size(report, 1))
    p1 = plot(abscissa, report.psiz, ylims=(0,maximum(report.psiz)),
        label="pop size", linewidth=3,
        title="A")
    p2 = plot(abscissa, [report.q1fit, report.meanfit, report.q2fit, report.maxfit], 
        label=["Q1 fitness" "mean fitness" "median fitness" "max fitness"], linewidth=3, title="B")
    p3 = plot(abscissa, [report.minfree, report.meanfree, report.maxfree], 
        label=["min ΔG" "mean ΔG" "max ΔG"], 
        linewidth=3, title="C")
    p4 = plot(abscissa, [report.minmut, report.meanmut, report.maxmut],
        label=["min # μ count" "mean # μ" "max # μ"], 
        linewidth=3, title="D")
    plot(p1, p2, p3, p4, titleloc = :left, titlefont = font(20), layout=(2,2), size=(1000, 700))
end
        
function synchronized_generation(populace)
    #@unpack R = p
    next_generation = []
    for parent in populace
        for r in 1:probabilistic_round(R * parent.fitness)
            child = reproduce(parent)
            ( child.fitness > 0  ) && ( push!(next_generation, child) ) # could be a zombie test
        end
    end
    return next_generation
end
        
function synchronized_simulation()
    #m = Params(description = "synchronized", replicate = '1')
    #@unpack K, N, sim_length = m
    population = initialize_population()
    population_size = N
    report = initialize_report()
    SIM_DURATION = sim_length
    while SIM_DURATION > 0
        SIM_DURATION -= 1
        if SIM_DURATION % 20 == 0 # can change reporting frequency here
            println("generation:", sim_length - SIM_DURATION)
        end
        report_update!(population, report)
        population = synchronized_generation(population)
        population_size = length(population)

        if population_size > K
            population = sample(population, K, replace=false)
            # why calculate when you can set?
            population_size = K

        elseif population_size == 0
            # they're all dead
            SIM_DURATION = 0
        end
    end
    return report
end

        
function run_sim(seed::Int)
	Random.seed!(seed)
    base_dir = "C:\\Users\\jade-\\Desktop\\Lethal Mutagenesis\\Simulation" 
    mkpath(base_dir)

    variable1 = "L$(L)"
    variable2 = "N$(N)"
    variable3 = "U$(U)"
    variable4 = "K$(K)"
    variable5 = "R$(R)"
    seed_str = "$seed"
    
    csv_filename = joinpath(base_dir, "$(variable1)_$(variable2)_$seed_str.csv")
    png_filename = joinpath(base_dir, "$(variable1)_$(variable2)_$seed_str.png")

    synchronized_report = synchronized_simulation()
    plot_simulation(synchronized_report)
            
    CSV.write(csv_filename, synchronized_report) 
    png(png_filename)

    return csv_filename, png_filename
end
    
num_runs = 5

for i in 1:num_runs
    seed = rand(Int) #Generate random seed per run
    println("Running simulation with seed: $seed")

    run_sim(seed)

    println("Simulation run $i complete.")
end

println("All simulation runs complete.")


