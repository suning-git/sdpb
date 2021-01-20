#include "../SDP_Solver.hxx"

// Create and initialize an SDPSolver for the given SDP and
// SDP_Solver_Parameters
SDP_Solver::SDP_Solver(const SDP_Solver_Parameters &parameters,
                       const Block_Info &block_info, const El::Grid &grid,
                       const size_t &dual_objective_b_height)
    : x(block_info.schur_block_sizes, block_info.block_indices,
        block_info.schur_block_sizes.size(), grid),
      X(block_info.psd_matrix_block_sizes, block_info.block_indices,
        block_info.schur_block_sizes.size(), grid),
      y(std::vector<size_t>(block_info.schur_block_sizes.size(),
                            dual_objective_b_height),
        block_info.block_indices, block_info.schur_block_sizes.size(), grid),
      Y(X), primal_residues(X),
      dual_residues(block_info.schur_block_sizes, block_info.block_indices,
                    block_info.schur_block_sizes.size(), grid),

	dx_backup(x), dX_backup(X), dy_backup(y), dY_backup(Y),

      current_generation(0)
{
  if(!load_checkpoint(parameters.checkpoint_in, block_info,
                      parameters.verbosity,
                      parameters.require_initial_checkpoint))
    {
      X.set_zero();
      Y.set_zero();
      for(auto &block : x.blocks)
        {
          Zero(block);
        }
      for(auto &block : y.blocks)
        {
          Zero(block);
        }

      // X = \Omega_p I
      X.add_diagonal(parameters.initial_matrix_scale_primal);
      // Y = \Omega_d I
      Y.add_diagonal(parameters.initial_matrix_scale_dual);
    }
}
