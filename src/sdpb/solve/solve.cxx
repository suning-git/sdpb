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

void compute_dsdp(SDP & dsdp, const Block_Info &block_info, const El::Grid &grid, const boost::filesystem::path & dsdp_filepath)
{
	SDP sdp2(dsdp_filepath, block_info, grid);

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


void compute_dsdp_mode3(SDP & dsdp, const Block_Info &block_info, const El::Grid &grid, const boost::filesystem::path & dsdp_filepath)
{
	SDP sdp2(dsdp_filepath, block_info, grid);

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



El::BigFloat compute_xBy(const Block_Info &block_info, const SDP &sdp,
	const Block_Vector &x, const Block_Vector &y)
{
	Block_Vector By(x);

	auto By_block(By.blocks.begin());
	auto primal_objective_c_block(sdp.primal_objective_c.blocks.begin());
	auto y_block(y.blocks.begin());
	auto free_var_matrix_block(sdp.free_var_matrix.blocks.begin());

	for (auto &block_index : block_info.block_indices)
	{
		// By = 0
		Zero(*By_block);
		const size_t block_size(block_info.degrees[block_index] + 1);

		// By -= FreeVarMatrix * y
		Gemm(El::Orientation::NORMAL, El::Orientation::NORMAL, El::BigFloat(1),
			*free_var_matrix_block, *y_block, El::BigFloat(1),
			*By_block);

		++y_block;
		++free_var_matrix_block;
		++By_block;
	}

	return dot(x, By);
}



El::BigFloat compute_hessian_component(const SDP &dsdp, const DSDPSOLUTION & sol, const SDP_Solver & solver, const Block_Info &block_info)
{

	El::BigFloat dbdy;
	if (!sol.dy.blocks.empty())
	{
		dbdy = El::Dotu(dsdp.dual_objective_b, sol.dy.blocks.front());
	}

	El::BigFloat dcdx = dot(dsdp.primal_objective_c, sol.dx);

	El::BigFloat dxdBy = compute_xBy(block_info, dsdp, sol.dx, solver.y);

	El::BigFloat xdBdy = compute_xBy(block_info, dsdp, solver.x, sol.dy);

	/*
	if (El::mpi::Rank() == 0)
	{
		std::cout << El::mpi::Rank() << " dbdy= " << dbdy << "\n";
		std::cout << El::mpi::Rank() << " dcdx= " << dcdx << "\n";
		std::cout << El::mpi::Rank() << " dxdBy= " << dxdBy << "\n";
		std::cout << El::mpi::Rank() << " xdBdy= " << xdBdy << "\n";
	}*/

	El::BigFloat rslt = 2*dbdy + dcdx - dxdBy - xdBdy;

	return rslt;
}


int initialize_dsdp(SDP & sdp0, std::vector<SDP*> & dsdp_list, const Block_Info &block_info, const El::Grid &grid, const SDP_Solver_Parameters &parameters)
{
	if (El::mpi::Rank() == 0)
	{
		std::cout << " Total # of dsdp : " << parameters.list_sdp2_path.size() << "\n";
	}

	for(int i=0;i<parameters.list_sdp2_path.size();i++)
	{
		if (El::mpi::Rank() == 0)
		{
			std::cout << " Reading #" << i << " dsdp file from " << parameters.list_sdp2_path[i] << "\n";
		}

		SDP*psdp = new SDP(sdp0);

		if (parameters.compute_derivative_dBdbdc)
		{
			if (parameters.sdpd_mode_dBdbdc)
			{
				compute_dsdp_mode3(*psdp, block_info, grid, parameters.list_sdp2_path[i]);
			}
			else
				compute_dsdp(*psdp, block_info, grid, parameters.list_sdp2_path[i]);
		}

		dsdp_list.push_back(psdp);
	}
	return 1;
}

El::BigFloat compute_derivative_Balt_formula(SDP & dsdp, SDP_Solver & solver, const Block_Info &block_info)
{
	El::BigFloat dby;
	if (!solver.y.blocks.empty())
	{
		dby = dsdp.objective_const + El::Dotu(dsdp.dual_objective_b, solver.y.blocks.front());
	}

	El::BigFloat xdBy = compute_xBy(block_info, dsdp, solver.x, solver.y);

	El::BigFloat dprimalobj_Balt = dot(dsdp.primal_objective_c, solver.x) + dby - xdBy;

	return dprimalobj_Balt;
}

int compute_gradient(std::vector<El::BigFloat> & dobj_list, std::vector<SDP*> & dsdp_list, SDP_Solver & solver, const Block_Info &block_info)
{
	for (int i = 0; i < dsdp_list.size(); i++)
	{
		dobj_list.push_back(compute_derivative_Balt_formula(*dsdp_list[i], solver, block_info));
	}
	return 1;
}

int compute_hessian(std::vector<std::vector<El::BigFloat>> & hessian, const std::vector<SDP*> & dsdp_list, const SDP_Solver & solver, const Block_Info &block_info)
{
	for (int i = 0; i < dsdp_list.size(); i++)
	{
		std::vector<El::BigFloat> vec;
		for (int j = 0; j < dsdp_list.size(); j++)
		{
			vec.push_back(compute_hessian_component(*dsdp_list[i], *solver.dsdp_sol_list[j], solver, block_info));
		}
		hessian.push_back(vec);
	}
	return 1;
}

Timers solve(const Block_Info &block_info, const SDP_Solver_Parameters &parameters)
{
  // Read an SDP from sdpFile and create a solver for it
  El::Grid grid(block_info.mpi_comm.value);

  Timers timers(parameters.verbosity >= Verbosity::debug);

  SDP sdp(parameters.sdp_directory, block_info, grid);

  std::vector<SDP*> dsdp_list;

  initialize_dsdp(sdp, dsdp_list, block_info, grid, parameters);

  //if (El::mpi::Rank() == 0) std::cout << "test A\n";

  SDP_Solver solver(parameters, block_info, grid,
                    sdp.dual_objective_b.Height());

  std::vector<El::BigFloat> dobj_list;
  compute_gradient(dobj_list, dsdp_list, solver, block_info);

  SDP_Solver_Terminate_Reason reason
    = solver.run(parameters, block_info, sdp, dsdp_list, grid, timers);

  std::vector<std::vector<El::BigFloat>> hessian;

  compute_hessian(hessian, dsdp_list, solver, block_info);

  El::BigFloat sdpd_cdx = dot(sdp.primal_objective_c, solver.dsdp_sol_list[0]->dx);
  El::BigFloat sdpd_dcx = dsdp_list[0]->objective_const + dot(dsdp_list[0]->primal_objective_c, solver.x);
  El::BigFloat dprimalobj = sdp.objective_const + dot(sdp.primal_objective_c, solver.dsdp_sol_list[0]->dx);

  El::BigFloat sdpd_bdy;
  if (!solver.dsdp_sol_list[0]->dy.blocks.empty())
  {
	  sdpd_bdy = El::Dotu(sdp.dual_objective_b, solver.dsdp_sol_list[0]->dy.blocks.front());
  }

  El::BigFloat sdpd_dby;
  if (!solver.y.blocks.empty())
  {
	  sdpd_dby = dsdp_list[0]->objective_const + El::Dotu(dsdp_list[0]->dual_objective_b, solver.y.blocks.front());
  }

  El::BigFloat xdBy = compute_xBy(block_info, *dsdp_list[0], solver.x, solver.y);

  El::BigFloat dprimalobj_Balt = dot(dsdp_list[0]->primal_objective_c, solver.x) + sdpd_dby - xdBy;

  if(parameters.verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
    {
      set_stream_precision(std::cout);

	  if (parameters.compute_derivative_dBdbdc)
	  {
		  std::cout << "[SDPDReturnBegin.Gradient]\n";
		  std::cout << "{\n";
		  for (int i = 0; i < dobj_list.size(); i++)
		  {
			  std::cout << dobj_list[i] << "\n";
			  if (i != dobj_list.size() - 1) std::cout << ",";
		  }
		  std::cout << "}\n";
		  std::cout << "[SDPDReturnEnd.Gradient]\n";

		  std::cout << "[SDPDReturnBegin.Hessian]\n";
		  std::cout << "{\n";
		  for (int i = 0; i < dsdp_list.size(); i++)
		  {
			  std::cout << "{\n";
			  for (int j = 0; j < dsdp_list.size(); j++)
			  {
				  std::cout << hessian[i][j];
				  if (j != dsdp_list.size() - 1) std::cout << ",";
			  }
			  std::cout << "}\n";
			  if (i != dsdp_list.size() - 1) std::cout << ",";
		  }
		  std::cout << "}\n";
		  std::cout << "[SDPDReturnEnd.Hessian]\n";

		  /*
		  std::cout << El::mpi::Rank() << " dsdp.objective_const= " << dsdp_list[0]->objective_const << "\n";
		  std::cout << El::mpi::Rank() << " sdp.objective_const= " << sdp.objective_const << "\n";

		  std::cout << El::mpi::Rank() << " dc.x= " << sdpd_dcx << "\n";
		  std::cout << El::mpi::Rank() << " c.dx= " << sdpd_cdx << "\n";
		  std::cout << El::mpi::Rank() << " d(c.x)= " << sdpd_cdx + sdpd_dcx << "\n";

		  std::cout << El::mpi::Rank() << " b.dy= " << sdpd_bdy << "\n";
		  std::cout << El::mpi::Rank() << " db.y= " << sdpd_dby << "\n";
		  std::cout << El::mpi::Rank() << " d(b.y)= " << sdpd_bdy + sdpd_dby << "\n";

		  std::cout << El::mpi::Rank() << " -x.dB.y= " << -xdBy << "\n";

		  std::cout << El::mpi::Rank() << " dc.x+db.y-x.dB.y= " << dprimalobj_Balt << "\n";
		  */
		  std::cout << "[SDPDReturnBegin.Gradient.NingFormula.Primal]" << sdpd_cdx + sdpd_dcx << "[SDPDReturnEnd.Gradient.NingFormula.Primal]\n";
		  std::cout << "[SDPDReturnBegin.Gradient.NingFormula.Dual]" << sdpd_bdy + sdpd_dby << "[SDPDReturnEnd.Gradient.NingFormula.Dual]\n";

		  std::cout << "[SDPDReturnBegin]" << sdpd_cdx + sdpd_dcx << "[SDPDReturnEnd]\n";
	  }
	  else
	  {
		  std::cout << El::mpi::Rank() << " dprimalobj= " << dprimalobj << "\n";

		  std::cout << "[SDPDReturnBegin]" << dprimalobj << "[SDPDReturnEnd]\n";
	  }

	  /*
      std::cout //<< "-----" << reason << "-----\n"
                << '\n'
                << "primalObjective = " << solver.primal_objective << '\n'
                << "dualObjective   = " << solver.dual_objective << '\n'
                << "dualityGap      = " << solver.duality_gap << '\n'
                << "primalError     = " << solver.primal_error() << '\n'
                << "dualError       = " << solver.dual_error << '\n'
                << '\n';
				*/
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
