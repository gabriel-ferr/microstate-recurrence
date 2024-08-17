#
#       By Gabriel Ferreira
#       Orientation: Prof. Thiago de Lima Prado
#       UFPR - 2024
#
#       Modify 'config.jl' if you need.
#       Run this script using 'julia main.jl'.
# ------------------------------------------------------------------
#       ** Required libraries:
#   - Local: RManager.jl
#   (all dependencies that you will need)
# ------------------------------------------------------------------
include("lib/rp.jl")
include("lib/run_manager.jl")
# ------------------------------------------------------------------
using .RManager
using Random
# ------------------------------------------------------------------
rng = MersenneTwister()
Random.seed!()
# ------------------------------------------------------------------
cfg = RManager.load_settings()
RManager.main(cfg)
# ------------------------------------------------------------------