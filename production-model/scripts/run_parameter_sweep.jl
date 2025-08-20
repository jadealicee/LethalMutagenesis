"""
Lethal Mutagenesis Parameter Sweep Script

Runs viral evolution simulations across parameter combinations.
Based on notebook: Epistasis_JuliaLM_V5.0_L_Loop.ipynb

Usage:
    julia --project=.. --threads=24 scripts/run_parameter_sweep.jl
"""

# Import necessary packages
using Random, Distributions, StatsBase, DataFrames, Plots, CSV, Tables
using Base.Threads
using ProgressMeter
#using LethalMutagenesis #later

# Universal physical constants and fixed data
const kt_boltzmann = 0.001987204118 * (273.15 + 37)
const ΔΔG = Normal(1.0, 1.7)

# Epistasis configuration
Base.@kwdef struct EpistasisParams
    # Module configuration
    module_sizes::Dict{Symbol,Int} = Dict(
        :auxiliary => 3,
        :assembly => 3,
        :replication => 3,
        :host_interaction => 3
    )

    # Calculated properties
    module_names::Vector{Symbol} = collect(keys(module_sizes))
    G::Int = sum(values(module_sizes)) # number_of_genes
end

# Constructor function (required to calculate G = number of genes)
function EpistasisParams(; module_sizes=Dict(:auxiliary => 3, :assembly => 3, :replication => 3, :host_interaction => 3),
    module_names=collect(keys(module_sizes)))
    G = sum(values(module_sizes))
    return EpistasisParams(module_sizes, module_names, G)
end

# Instantiate with equal module sizes (the default)
equal_epi_config = EpistasisParams()

# Custom values reflecting phiX174 (unequal module sizes)
ΦX174_epi_config = EpistasisParams(
    module_sizes=Dict(
        :auxiliary => 2, # A* (blocks super-inf), K (optimises burst size) - overlaps, non-essential
        :assembly => 3, # B (internal scaffolding), D (external scaffolding), F (major capsid) - but H involved in stoichiometry
        :replication => 2, # A (rolling circle init), J (ssDNA binding)
        :host_interaction => 4 # C (DNA packaging), E (host cell lysis), G (major spike), H (DNA pilot/minor spike)
    )
)

# Custom values to remove epistasis and make everything a product
null_epi_config = EpistasisParams(
    module_sizes=Dict(
        :auxiliary => 1,
        :assembly => 1,
        :replication => 1,
        :host_interaction => 8
    )
)

# Variable parameters for the simulation
Base.@kwdef struct SimulationParams
    # Control parameters
    sim_length::Int = 500
    F_init::Float64 = -5.0
    # Sweep parameters
    U_default::Poisson = Poisson(5.0)
    L_default::Float64 = 0.1
    N_default::Int = 5000
    K_default::Int = 10000
    R_default::Int = 9
end

# Instantiate with default values
sim_config = SimulationParams()  # All defaults

# Core Virus struct
mutable struct ModularVirus
    μ_counts::Dict{Symbol,Vector{Int64}}
    ΔG_values::Dict{Symbol,Vector{Float64}}
    fitness::Float64
end

# constructor for clean initialization
function ModularVirus(sim_config::SimulationParams, epi_config::EpistasisParams, start_fitness::Float64)
    μ_counts = Dict(name => zeros(Int64, epi_config.module_sizes[name]) for name in epi_config.module_names)
    ΔG_values = Dict(name => fill(sim_config.F_init, epi_config.module_sizes[name]) for name in epi_config.module_names)
    return ModularVirus(μ_counts, ΔG_values, start_fitness)
end

# Epistasis Calulcation Functions
function auxiliary_epistasis(auxiliary_ΔGs::Vector{Float64})
    # Additive contributions
    stabilities = [1 / (1 + ℯ^(ΔG / kt_boltzmann)) for ΔG in auxiliary_ΔGs]
    # Arithmetic mean benefits from positives
    return mean(stabilities)
end

function assembly_epistasis(assembly_ΔGs::Vector{Float64})
    # Stoichiometric balance: variance in stabilities hurts fitness
    stabilities = [1 / (1 + ℯ^(ΔG / kt_boltzmann)) for ΔG in assembly_ΔGs]
    # Geometric mean penalizes imbalance
    return prod(stabilities)^(1 / length(stabilities))
end

function replication_epistasis(replication_ΔGs::Vector{Float64})
    # Rate-limiting: weakest link dominates
    stabilities = [1 / (1 + ℯ^(ΔG / kt_boltzmann)) for ΔG in replication_ΔGs]
    return minimum(stabilities)
end

function host_module_epistasis(host_ΔGs::Vector{Float64})
    # Independent effects but negatives propagate
    stabilities = [1 / (1 + ℯ^(ΔG / kt_boltzmann)) for ΔG in host_ΔGs]
    return prod(stabilities)
end

# Calculate fitness based on epistasis
function update_fitness!(virus::ModularVirus)
    auxiliary_contrib = auxiliary_epistasis(virus.ΔG_values[:auxiliary])
    assembly_contrib = assembly_epistasis(virus.ΔG_values[:assembly])
    replication_contrib = replication_epistasis(virus.ΔG_values[:replication])
    host_contrib = host_module_epistasis(virus.ΔG_values[:host_interaction])

    virus.fitness = auxiliary_contrib * assembly_contrib * replication_contrib * host_contrib
end

# Mutation Function
function mutate!(virus::ModularVirus, sim_config::SimulationParams, epi_config::EpistasisParams)
    number_of_mutations = only(rand(sim_config.U_default, 1))
    ΔΔG_values = rand(ΔΔG, number_of_mutations)  # note that ΔΔG is a global

    # We choose modules proportionally to their size **and ordered by** module_names vector
    # care must be taken as keys(module_sizes) != module_names
    module_probs = Weights([epi_config.module_sizes[name] for name in epi_config.module_names] ./ epi_config.G)

    for i in 1:number_of_mutations
        # Select module based on proportional targeting
        module_choice = sample(epi_config.module_names, module_probs)

        # Select gene within module
        module_size = epi_config.module_sizes[module_choice]
        gene_index = rand(1:module_size)

        # Apply mutation to appropriate dictionary [and module therein][and gene]
        virus.μ_counts[module_choice][gene_index] += 1
        virus.ΔG_values[module_choice][gene_index] += ΔΔG_values[i]

        # Lethal mutation check
        if rand() < sim_config.L_default
            virus.fitness = 0
            return virus.fitness
        end
    end

    # Update fitness if virus is still viable
    (virus.fitness > 0) && (update_fitness!(virus))

    return virus.fitness
end

# Reproduction Function
# Creates new offspring virus by copying parent + its mutations
function reproduce(parent::ModularVirus, sim_config::SimulationParams, epi_config::EpistasisParams)
    sprog = deepcopy(parent)
    mutate!(sprog, sim_config, epi_config)
    return sprog
end

# Creates initial virus population
function initialize_population(sim_config::SimulationParams, epi_config::EpistasisParams)
    start_fitness = (1 / (1 + ℯ^(sim_config.F_init / kt_boltzmann)))^epi_config.G
    return [ModularVirus(sim_config, epi_config, start_fitness) for _ in 1:sim_config.N_default]
end

# Helper Functions
# Used to implement fitness-proportional selection
function get_weights(populace)
    fitness_values = [v.fitness for v in populace]
    total = sum(fitness_values)
    fitness_values ./= total  # in-place division to avoid new array
    return Weights(fitness_values)
end

# Round a floating-point number to an integer
function probabilistic_round(number)
    base = floor(Int, number)
    return base + (rand() < number - base)
    # bracketed boolean expression evaluates to 0 or 1
end

# Helper function to get all ΔG values from a ModularVirus
function get_all_ΔG(virus::ModularVirus, epi_config::EpistasisParams)
    all_ΔG = Float64[]
    for module_name in epi_config.module_names
        append!(all_ΔG, virus.ΔG_values[module_name])
    end
    return all_ΔG
end

# Helper function to get all mutation counts from a ModularVirus  
function get_all_mutations(virus::ModularVirus, epi_config::EpistasisParams)
    return sum(sum(virus.μ_counts[module_name]) for module_name in epi_config.module_names)
end

# Helper function to get module-specific mean statistics
function get_module_stats(populace, module_name::Symbol)
    mut_counts = [sum(v.μ_counts[module_name]) for v in populace]
    ΔG_values = [mean(v.ΔG_values[module_name]) for v in populace]
    return mean(mut_counts), mean(ΔG_values)
end

# Reporting Functions
# Creates empty DataFrame with module-specific columns
function initialize_report()
    report = DataFrame(
        # Main metrics
        generation=Int[],
        psiz=Int[],
        q1fit=Union{Float64,Missing}[],
        meanfit=Union{Float64,Missing}[],
        q2fit=Union{Float64,Missing}[],
        maxfit=Union{Float64,Missing}[],
        minfree=Union{Float64,Missing}[],
        meanfree=Union{Float64,Missing}[],
        maxfree=Union{Float64,Missing}[],
        minmut=Union{Float64,Missing}[],
        meanmut=Union{Float64,Missing}[],
        maxmut=Union{Float64,Missing}[],
        # Module-specific mutation counts
        auxiliary_meanmut=Union{Float64,Missing}[],
        assembly_meanmut=Union{Float64,Missing}[],
        replication_meanmut=Union{Float64,Missing}[],
        host_meanmut=Union{Float64,Missing}[],
        # Module-specific free energies  
        auxiliary_meanfree=Union{Float64,Missing}[],
        assembly_meanfree=Union{Float64,Missing}[],
        replication_meanfree=Union{Float64,Missing}[],
        host_meanfree=Union{Float64,Missing}[]
    )
    return report
end

# Updates DataFrame with both original and module-specific metrics
function report_update!(populace, report, epi_config::EpistasisParams, generation::Int)

    if isempty(populace)
        # Handle extinction case - population size is 0, everything else is missing
        push!(report, (
            generation=generation,
            psiz=0,  # Population size is 0
            q1fit=missing,
            meanfit=missing,
            q2fit=missing,
            maxfit=missing,
            minfree=missing,
            meanfree=missing,
            maxfree=missing,
            minmut=missing,
            meanmut=missing,
            maxmut=missing,
            auxiliary_meanmut=missing,
            assembly_meanmut=missing,
            replication_meanmut=missing,
            host_meanmut=missing,
            auxiliary_meanfree=missing,
            assembly_meanfree=missing,
            replication_meanfree=missing,
            host_meanfree=missing
        ))
        return
    end

    # Original metrics
    fitness_values = [v.fitness for v in populace]
    all_ΔG_per_virus = [get_all_ΔG(v, epi_config) for v in populace]
    all_mutations_per_virus = [get_all_mutations(v, epi_config) for v in populace]

    # Module-specific metrics
    aux_mut, aux_free = get_module_stats(populace, :auxiliary)
    asm_mut, asm_free = get_module_stats(populace, :assembly)
    rep_mut, rep_free = get_module_stats(populace, :replication)
    host_mut, host_free = get_module_stats(populace, :host_interaction)

    # Precompute ΔG statistics for efficiency
    ΔG_mins = [minimum(ΔG_list) for ΔG_list in all_ΔG_per_virus]
    ΔG_means = [mean(ΔG_list) for ΔG_list in all_ΔG_per_virus]
    ΔG_maxs = [maximum(ΔG_list) for ΔG_list in all_ΔG_per_virus]

    push!(report, (
        generation=generation,
        psiz=length(populace),
        q1fit=quantile(fitness_values, 0.25),
        meanfit=mean(fitness_values),
        q2fit=median(fitness_values),
        maxfit=maximum(fitness_values),
        minfree=mean(ΔG_mins),
        meanfree=mean(ΔG_means),
        maxfree=mean(ΔG_maxs),
        minmut=minimum(all_mutations_per_virus),
        meanmut=mean(all_mutations_per_virus),
        maxmut=maximum(all_mutations_per_virus),
        auxiliary_meanmut=aux_mut,
        assembly_meanmut=asm_mut,
        replication_meanmut=rep_mut,
        host_meanmut=host_mut,
        auxiliary_meanfree=aux_free,
        assembly_meanfree=asm_free,
        replication_meanfree=rep_free,
        host_meanfree=host_free
    ))
end

#creates plots including module-specific panels
function plot_simulation(report)
    abscissa = 1:size(report, 1)

    p1 = plot(abscissa, report.psiz,
        ylims=(0, maximum(report.psiz)),
        label="pop size", linewidth=3, title="A) Population Size")

    p2 = plot(abscissa, [report.q1fit report.meanfit report.q2fit report.maxfit],
        label=["Q1 fitness" "mean fitness" "median fitness" "max fitness"],
        linewidth=3, title="B) Fitness Distribution")

    p3 = plot(abscissa, [report.minfree report.meanfree report.maxfree],
        label=["min ΔG" "mean ΔG" "max ΔG"],
        linewidth=3, title="C) Free Energy Distribution")

    p4 = plot(abscissa, [report.minmut report.meanmut report.maxmut],
        label=["min mutations" "mean mutations" "max mutations"],
        linewidth=3, title="D) Total Mutations")

    p5 = plot(abscissa, [report.replication_meanmut report.assembly_meanmut report.host_meanmut],
        label=["replication" "assembly" "host interaction"],
        linewidth=3, title="E) Mutations by Module")

    p6 = plot(abscissa, [report.replication_meanfree report.assembly_meanfree report.host_meanfree],
        label=["replication" "assembly" "host interaction"],
        linewidth=3, title="F) Mean ΔG by Module")

    plot(p1, p2, p3, p4, p5, p6,
        titleloc=:left,
        titlefont=font(12),
        layout=(3, 2),
        size=(1000, 1000))
end

# Simulation Function
# Simulates one generation of the population. fitness pushed to next generation
function synchronized_generation(populace, sim_config::SimulationParams, epi_config::EpistasisParams)
    next_generation = ModularVirus[]
    # equivalent to Array{ModularVirus,1}() and eliminates type checking later

    for parent in populace
        offspring_count = probabilistic_round(sim_config.R_default * parent.fitness)
        for r in 1:offspring_count
            child = reproduce(parent, sim_config, epi_config)
            # Only add viable offspring (short-circuit evaluation)
            (child.fitness > 0) && push!(next_generation, child)
        end
    end

    return next_generation
end

# Runs main sim and collects data into DataFrame
function synchronized_simulation(sim_config::SimulationParams, epi_config::EpistasisParams;
    report_frequency=20)
    # Initialize population and reporting
    population = initialize_population(sim_config, epi_config)
    report = initialize_report()

    generation = 0
    max_generations = sim_config.sim_length

    while generation < max_generations
        generation += 1

        # Progress reporting
        if generation % report_frequency == 0
            println("Generation: $generation, Population size: $(length(population))")
        end

        # Collect data before selection/reproduction
        report_update!(population, report, epi_config, generation)

        # Generate next generation
        population = synchronized_generation(population, sim_config, epi_config)
        population_size = length(population)

        # Apply carrying capacity constraint
        if population_size > sim_config.K_default
            population = sample(population, sim_config.K_default, replace=false)
        elseif population_size == 0
            println("Population extinct at generation $generation")
            report_update!([], report, epi_config, generation)
            break
        end
    end

    # Final report update if simulation completed
    if length(population) > 0
        report_update!(population, report, epi_config, generation)
    end

    return report, population
end

# Run Parameter Sweep in Parallel
# Runs a parameter sweep across multiple configurations in parallel
# Saves results to CSV files in a specified directory
# Uses a fixed number of repetitions for each parameter combination
# Uses threads for parallel execution
function run_parameter_sweep_parallel(base_dir::String, epi_config::EpistasisParams;
    num_reps=5,
    R_range=0:2:10,
    U_range=0:0.5:5.0,
    L_range=0.0:0.2:1.0)

    mkpath(base_dir)

    # Pre-generate all parameter combinations and seeds
    param_combinations = []
    for R in R_range, U in U_range, L in L_range, rep in 1:num_reps
        seed = rand(Int)
        push!(param_combinations, (R=R, U=U, L=L, rep=rep, seed=seed))
    end

    total_sims = length(param_combinations)
    println("Running $total_sims simulations on $(nthreads()) threads...")

    # Thread-safe storage for results
    all_results = Vector{DataFrame}(undef, total_sims)

    # Create progress bar
    progress = Progress(total_sims, desc="Simulations: ", barlen=50, color=:green)

    # Parallel execution
    @threads for i in 1:total_sims
        params = param_combinations[i]

        # Create parameter-specific directory (thread-safe)
        param_dir = joinpath(base_dir, "R$(params.R)_U$(params.U)_L$(params.L)")
        mkpath(param_dir)

        # Create simulation config
        sweep_sim_config = SimulationParams(
            R_default=params.R,
            U_default=Poisson(params.U),
            L_default=params.L
        )

        # Set seed for this specific simulation
        Random.seed!(params.seed)

        # Run simulation
        report, final_pop = synchronized_simulation(sweep_sim_config, epi_config; report_frequency=999999)

        # Add metadata
        insertcols!(report, 1,
            :R => params.R,
            :U => params.U,
            :L => params.L,
            :N => sim_config.N_default,
            :K => sim_config.K_default,
            :replicate => params.rep,
            :seed => params.seed
        )

        # Save CSV file
        csv_filename = joinpath(param_dir, "rep$(params.rep)_seed$(params.seed).csv")
        CSV.write(csv_filename, report)

        # Store result for later aggregation
        all_results[i] = report

        # Update progress bar (thread-safe)
        next!(progress)
        #println("Completed: R=$(params.R), U=$(params.U), L=$(params.L), rep=$(params.rep) on thread $(threadid())")
    end

    # Combine all results
    master_results = vcat(all_results...)

    # Save master results
    master_csv = joinpath(base_dir, "master_results.csv")

    if isfile(master_csv)
        # Append to existing file
        CSV.write(master_csv, master_results, append=true, header=false)
        println("\nParameter sweep complete! Results appended to existing: $master_csv")
        println("Added $(nrow(master_results)) new rows to master file")
    else
        # Create new file
        CSV.write(master_csv, master_results)
        println("\nParameter sweep complete! New master results saved to: $master_csv")
        println("Created master file with $(nrow(master_results)) rows")
    end

    return master_results
end

function main()
    # Run the parameter sweep with default configurations
    results = run_parameter_sweep_parallel("parameter_sweep_results", ΦX174_epi_config;
        num_reps=5,
        R_range=0:2:10,      # 6 values
        U_range=0:0.5:5.0,   # 11 values
        L_range=0.0:0.2:1.0) # 6 values

    println("Parameter sweep completed successfully!")
end
if abspath(PROGRAM_FILE) == @__FILE__
    main()
else
    println("This script is intended to be run directly, not included as a module.")
end
# This allows the script to be run directly or imported without executing main()