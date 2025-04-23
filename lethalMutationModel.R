# Basic Simulation
## parameters
S = 30 # simulation length
N = 500 # starting population size
K = 500 # carrying capacity (max pop)
U = 1.5 # lethal mutation rate (per genome, per generation)
R = 3 # baseline fecundity (before lethality)

## example - 500 individuals reproduce under lethal mutation
## expected number of offspring:
500 * R * exp(-U)

## function to give next generation
increment_pop <- function(population, fecundity) {
  next_generation <- population * fecundity * exp(-U)
  return(next_generation)
}

## use function to run simulation
generations <- vector(mode = "numeric", length = S)
generations[1] <- N
for (I in 2:S) {
  generations[I] =  min(K, increment_pop(generations[I-1], R))
}

plot(generations, type = "b")

# More Complex Simulation
## required (for graphics and paste0)
library(tidyverse)
library(RColorBrewer)

## parameters
S = 30 # simulation length
N = 500 # starting population size
K = 500 # carrying capacity (max pop)
U = 1.5 # lethal mutation rate (per genome, per generation)
R = seq(0, 5, by=0.5) # baseline fecundity (before lethality)
RCrit = 1/exp(-U)
R <- sort(c(R, RCrit))
NumPops = length(R)

## initiate pre-filled data array
generations <- array(data = NA, dim = c(S, NumPops+1)) # S rows and P+1 cols
colnames(generations) <- c("Generation", str_c("R", R)) # can use paste0 too
generations[,1] <- 1:S
generations[1, 1:NumPops+1] <- N # why is this 1: ??

## take a look
generations # S (30) generations (rows) and 10 populations 

for (fec in 1:NumPops) { # looping through R indices
  for (gen in 2:S) { # looping through generations (gen var)
    generations[gen, fec+1] <- min(K, increment_pop(generations[gen-1, fec+1], R[fec]))
  }
}

# take a look at output
generations

# reshape output for visualisation
longGen <- generations |> 
  as_tibble() |> 
  pivot_longer(cols = R0:last_col(), names_to = "Fecundity", values_to = "Population") |>
  mutate(Fecundity = parse_factor(str_remove(Fecundity, "R"), ordered = TRUE))

# examine
longGen

# graph
ggplot(longGen, aes(x = Generation, y = Population, colour = Fecundity)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = brewer.pal(NumPops, "Set3")) # "Spectral"

# save graph
ggsave("lethalMutationResult.pdf")


# More Advanced - Cycling through U Values as well

## parameters
R = c( seq(0, 3.5, by=0.5), 3.9, 3.99, 4.0)
NumPops = length(R)
U = seq(0, 2.5, by=0.25)
UCrit = -log(1/4) # this is for R=4
U <- sort(c(U, UCrit))
NumSims = length(U)

## function to give next generation, complicated
increment_pop <- function(population, fecundity, mutrate) {
  next_generation <- population * fecundity * exp(-mutrate)
  return(next_generation)
}

## initiate pre-filled data array
generations <- array(data = NA, dim = c(S, NumPops+1, NumSims)) # S rows and P+1 cols
colnames(generations) <- c("Generation", str_c("R", R)) # can use paste0 too
dimnames(generations)[[3]] <- as.character(U)
generations[,1,] <- 1:S
generations[1, 1:NumPops+1,] <- N # why is this 1: ??

## take a look
generations # S (30) generations (rows) and 10 populations 

for (mutidx in 1:NumSims) { # looping through U indices
  for (fecidx in 1:NumPops) { # looping through R indices
    for (gen in 2:S) { # looping through generations (gen var)
      generations[gen, fecidx+1, mutidx] <- min(K, increment_pop(generations[gen-1, fecidx+1, mutidx], R[fecidx], U[mutidx]))
    }
  }
}

# take a look at output
generations

# reshape output for visualisation (thanks to Claude for tidyverse solution)
longGen <- generations |> 
  # First convert array to a list of matrices, one for each U value
  array_branch(margin = 3) |>  
  # Name each matrix with its corresponding U value
  set_names(dimnames(generations)[[3]]) |>
  # Convert each matrix to a tibble, binding them together with its U value
  map_dfr(as_tibble, .id = "U") |>
  # Now use pivot_longer as before on the tibble data
  pivot_longer(
    cols = starts_with("R"), 
    names_to = "Fecundity", 
    values_to = "Population"
  ) |>
  mutate(
    Fecundity = parse_factor(str_remove(Fecundity, "R"), ordered = TRUE),
    U = as.numeric(U)  # Convert U back to numeric from character
  )

# examine
longGen

# graph
ggplot(longGen, aes(x = Generation, y = Population, colour = Fecundity)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = brewer.pal(NumPops, "Spectral")) +
  facet_wrap(~U) +
  theme_bw()

# save graph
ggsave("lethalMutationResultAll.pdf")


