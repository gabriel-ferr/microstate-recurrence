#
#       By Gabriel Ferreira
#       Orientation: Prof. Thiago de Lima Prado
#                   and Prof. Sérgio Roberto Lopes
#       UFPR - 2024
#
#       Modify 'config.jl' if you need.
#       Run this script using 'julia main.jl'.
# ----------------------------------------------------------------------------------------- #
#       ** Required libraries:                                                              #
#   - Local: RManager.jl                                                                    #
#   (all dependencies that the scripts will need)                                           #
# ----------------------------------------------------------------------------------------- #
include("lib/rp.jl")                                                                        #
include("lib/run_manager.jl")                                                               #
# ----------------------------------------------------------------------------------------- #
using .RManager                                                                             #
# ----------------------------------------------------------------------------------------- #
#       Load the Run Manager configuration and call up its main                             #
#   operation.                                                                              #
cfg = RManager.load_settings()                                                              #
RManager.main(cfg)                                                                          #
# ----------------------------------------------------------------------------------------- #