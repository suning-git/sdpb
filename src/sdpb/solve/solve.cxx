//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#include "SDP_Solver.hxx"
#include "../../Timers.hxx"
#include "../../set_stream_precision.hxx"

#include <El.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>


El::BigFloat dot(const Block_Vector &a, const Block_Vector &b);

void compute_dsdp(SDP & dsdp, const Block_Info &block_info, const El::Grid &grid, const SDP_Solver_Parameters &parameters)
{
	SDP sdp2(parameters.sdp2_directory, block_info, grid);

	dsdp.objective_const -= sdp2.objective_const;

	auto primal_objective_c_block(dsdp.primal_objective_c.blocks.begin());
	auto primal_objective_c2_block(sdp2.primal_objective_c.blocks.begin());
	for (auto &block_index : block_info.block_indices)
	{
		*primal_objective_c_block -= *primal_objective_c2_block;
		++primal_objective_c_block;
		++primal_objective_c2_block;
	}

	dsdp.dual_objective_b -= sdp2.dual_objective_b;

	auto free_var_matrix_block(dsdp.free_var_matrix.blocks.begin());
	auto free_var_matrix2_block(sdp2.free_var_matrix.blocks.begin());
	for (auto &block_index : block_info.block_indices)
	{
		*free_var_matrix_block -= *free_var_matrix2_block;
		++free_var_matrix_block;
		++free_var_matrix2_block;
	}

}


void compute_dsdp_mode3(SDP & dsdp, const Block_Info &block_info, const El::Grid &grid, const SDP_Solver_Parameters &parameters)
{
	SDP sdp2(parameters.sdp2_directory, block_info, grid);

	dsdp.objective_const = sdp2.objective_const;

	auto primal_objective_c_block(dsdp.primal_objective_c.blocks.begin());
	auto primal_objective_c2_block(sdp2.primal_objective_c.blocks.begin());
	for (auto &block_index : block_info.block_indices)
	{
		*primal_objective_c_block = *primal_objective_c2_block;
		++primal_objective_c_block;
		++primal_objective_c2_block;
	}

	dsdp.dual_objective_b = sdp2.dual_objective_b;

	auto free_var_matrix_block(dsdp.free_var_matrix.blocks.begin());
	auto free_var_matrix2_block(sdp2.free_var_matrix.blocks.begin());
	for (auto &block_index : block_info.block_indices)
	{
		*free_var_matrix_block = *free_var_matrix2_block;
		++free_var_matrix_block;
		++free_var_matrix2_block;
	}

}




Timers solve(const Block_Info &block_info, const SDP_Solver_Parameters &parameters)
{
  // Read an SDP from sdpFile and create a solver for it
  El::Grid grid(block_info.mpi_comm.value);

  Timers timers(parameters.verbosity >= Verbosity::debug);

  SDP sdp(parameters.sdp_directory, block_info, grid);

  SDP dsdp(sdp);

  if(parameters.compute_derivative_dBdbdc)
  { 
	  if (parameters.sdpd_mode_dBdbdc)
	  {
		  compute_dsdp_mode3(dsdp, block_info, grid, parameters);
	  }
	  else
		  compute_dsdp(dsdp, block_info, grid, parameters);
  }

  SDP_Solver solver(parameters, block_info, grid,
                    sdp.dual_objective_b.Height());

  SDP_Solver_Terminate_Reason reason
    = solver.run(parameters, block_info, sdp, dsdp, grid, timers);

  El::BigFloat sdpd_cdx = dot(sdp.primal_objective_c, solver.dx_backup);
  El::BigFloat sdpd_dcx = dsdp.objective_const + dot(dsdp.primal_objective_c, solver.x);
  El::BigFloat dprimalobj = sdp.objective_const + dot(sdp.primal_objective_c, solver.dx_backup);

  if(parameters.verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
    {
      set_stream_precision(std::cout);

	  if (parameters.compute_derivative_dBdbdc)
	  {
		  std::cout << El::mpi::Rank() << " dc.x= " << sdpd_dcx << "\n";
		  std::cout << El::mpi::Rank() << " c.dx= " << sdpd_cdx << "\n";
		  std::cout << El::mpi::Rank() << " d(c.x)= " << sdpd_cdx + sdpd_dcx << "\n";

		  std::cout << "[SDPDReturnBegin]" << sdpd_cdx + sdpd_dcx << "[SDPDReturnEnd]\n";
	  }
	  else
	  {
		  std::cout << El::mpi::Rank() << " dprimalobj= " << dprimalobj << "\n";

		  std::cout << "[SDPDReturnBegin]" << dprimalobj << "[SDPDReturnEnd]\n";
	  }

      std::cout //<< "-----" << reason << "-----\n"
                << '\n'
                << "primalObjective = " << solver.primal_objective << '\n'
                << "dualObjective   = " << solver.dual_objective << '\n'
                << "dualityGap      = " << solver.duality_gap << '\n'
                << "primalError     = " << solver.primal_error() << '\n'
                << "dualError       = " << solver.dual_error << '\n'
                << '\n';
    }

  if (!parameters.no_final_checkpoint)
  {
	  solver.save_solution(reason, timers.front(), parameters.out_directory,
		  parameters.write_solution,
		  block_info.block_indices,
		  parameters.verbosity);
  }

  /*
  if(!parameters.no_final_checkpoint)
    {
      solver.save_checkpoint(parameters);
    }*/

  return timers;
}
