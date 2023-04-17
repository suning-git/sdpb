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


void sdp_substract(SDP&sdp1, const SDP&sdp2, const Block_Info &block_info)
{
	for (int i = 0; i < sdp1.bilinear_bases_local.size(); i++)
		sdp1.bilinear_bases_local[i] -= sdp2.bilinear_bases_local[i];
	for (int i = 0; i < sdp1.bilinear_bases_dist.size(); i++)
		sdp1.bilinear_bases_dist[i] -= sdp2.bilinear_bases_dist[i];

	sdp1.objective_const -= sdp2.objective_const;

	auto primal_objective_c_block(sdp1.primal_objective_c.blocks.begin());
	auto primal_objective_c2_block(sdp2.primal_objective_c.blocks.begin());
	for (auto &block_index : block_info.block_indices)
	{
		*primal_objective_c_block -= *primal_objective_c2_block;
		++primal_objective_c_block;
		++primal_objective_c2_block;
	}
	
	sdp1.dual_objective_b -= sdp2.dual_objective_b;

	auto free_var_matrix_block(sdp1.free_var_matrix.blocks.begin());
	auto free_var_matrix2_block(sdp2.free_var_matrix.blocks.begin());
	for (auto &block_index : block_info.block_indices)
	{
		*free_var_matrix_block -= *free_var_matrix2_block;
		++free_var_matrix_block;
		++free_var_matrix2_block;
	}
	return;
}


void sdp_add(SDP&sdp1, const SDP&sdp2, const Block_Info &block_info)
{
	for (int i = 0; i < sdp1.bilinear_bases_local.size(); i++)
		sdp1.bilinear_bases_local[i] += sdp2.bilinear_bases_local[i];
	for (int i = 0; i < sdp1.bilinear_bases_dist.size(); i++)
		sdp1.bilinear_bases_dist[i] += sdp2.bilinear_bases_dist[i];

	sdp1.objective_const += sdp2.objective_const;

	auto primal_objective_c_block(sdp1.primal_objective_c.blocks.begin());
	auto primal_objective_c2_block(sdp2.primal_objective_c.blocks.begin());
	for (auto &block_index : block_info.block_indices)
	{
		*primal_objective_c_block += *primal_objective_c2_block;
		++primal_objective_c_block;
		++primal_objective_c2_block;
	}

	sdp1.dual_objective_b += sdp2.dual_objective_b;

	auto free_var_matrix_block(sdp1.free_var_matrix.blocks.begin());
	auto free_var_matrix2_block(sdp2.free_var_matrix.blocks.begin());
	for (auto &block_index : block_info.block_indices)
	{
		*free_var_matrix_block += *free_var_matrix2_block;
		++free_var_matrix_block;
		++free_var_matrix2_block;
	}
	return;
}


void sdp_equate(SDP&sdp1, const SDP&sdp2, const Block_Info &block_info)
{
	for (int i = 0; i < sdp1.bilinear_bases_local.size(); i++)
		sdp1.bilinear_bases_local[i] = sdp2.bilinear_bases_local[i];
	for (int i = 0; i < sdp1.bilinear_bases_dist.size(); i++)
		sdp1.bilinear_bases_dist[i] = sdp2.bilinear_bases_dist[i];

	sdp1.objective_const = sdp2.objective_const;

	auto primal_objective_c_block(sdp1.primal_objective_c.blocks.begin());
	auto primal_objective_c2_block(sdp2.primal_objective_c.blocks.begin());
	for (auto &block_index : block_info.block_indices)
	{
		*primal_objective_c_block = *primal_objective_c2_block;
		++primal_objective_c_block;
		++primal_objective_c2_block;
	}

	sdp1.dual_objective_b = sdp2.dual_objective_b;

	auto free_var_matrix_block(sdp1.free_var_matrix.blocks.begin());
	auto free_var_matrix2_block(sdp2.free_var_matrix.blocks.begin());
	for (auto &block_index : block_info.block_indices)
	{
		*free_var_matrix_block = *free_var_matrix2_block;
		++free_var_matrix_block;
		++free_var_matrix2_block;
	}
	return;
}

// multiply sdp1*scalar
void sdp_multiply(SDP&sdp1, const El::BigFloat scalar, const Block_Info &block_info)
{
	for (int i = 0; i < sdp1.bilinear_bases_local.size(); i++)
		sdp1.bilinear_bases_local[i] *= scalar;
	for (int i = 0; i < sdp1.bilinear_bases_dist.size(); i++)
		sdp1.bilinear_bases_dist[i] *= scalar;

	sdp1.dual_objective_b *= scalar;

	sdp1.objective_const *= scalar;

	auto primal_objective_c_block(sdp1.primal_objective_c.blocks.begin());
	for (auto &block_index : block_info.block_indices)
	{
		*primal_objective_c_block *= scalar;
		++primal_objective_c_block;
	}

	auto free_var_matrix_block(sdp1.free_var_matrix.blocks.begin());
	for (auto &block_index : block_info.block_indices)
	{
		*free_var_matrix_block *= scalar;
		++free_var_matrix_block;
	}
	return;
}

// set sdp1 to 0
void sdp_zero(SDP&sdp1, const Block_Info &block_info)
{
	sdp_multiply(sdp1,0, block_info);
}


void compute_dsdp(SDP & dsdp, const Block_Info &block_info, const El::Grid &grid, const boost::filesystem::path & dsdp_filepath, const bool exact_dSDP_Q, const El::BigFloat & finiteDifference)
{
	SDP sdp2(dsdp_filepath, block_info, grid);

	if (exact_dSDP_Q)
	{
		sdp_equate(dsdp,sdp2, block_info);
	}
	else
	{
		sdp_substract(dsdp, sdp2, block_info);
		sdp_multiply(dsdp, 1 / finiteDifference, block_info);
	}
}

void compute_dsdp_only_b(SDP & sdp0, SDP & dsdp, const Block_Info &block_info, const El::Grid &grid, const boost::filesystem::path & dsdp_filepath, const bool exact_dSDP_Q, const El::BigFloat & finiteDifference)
{
	SDP sdp2(dsdp_filepath, block_info, grid, true);

	if (exact_dSDP_Q)
	{
		sdp_zero(dsdp, block_info);
		dsdp.dual_objective_b = sdp2.dual_objective_b;
		dsdp.objective_const = sdp2.objective_const;
	}
	else
	{
		sdp_zero(dsdp, block_info);
		dsdp.dual_objective_b = sdp0.dual_objective_b;
		dsdp.dual_objective_b -= sdp2.dual_objective_b;
		dsdp.objective_const = sdp0.objective_const;
		dsdp.objective_const -= sdp2.objective_const;

		sdp_multiply(dsdp, 1 / finiteDifference, block_info);
	}
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

		if (parameters.sdpd_mode_dSDP_only_b)
		{
			compute_dsdp_only_b(sdp0, *psdp, block_info, grid, parameters.list_sdp2_path[i], parameters.sdpd_mode_exact_dSDP, parameters.finiteDifference);
		}
		else
		{
			compute_dsdp(*psdp, block_info, grid, parameters.list_sdp2_path[i], parameters.sdpd_mode_exact_dSDP, parameters.finiteDifference);
		}

		dsdp_list.push_back(psdp);
	}
	return 1;
}

// see sdpb study.nb  Section "2nd order pertubation" for the logic of this implementation
int initialize_ddsdp_ij(const SDP &sdp00, SDP* & dsdp_i, SDP* & dsdp_j, SDP* & ddsdp_ij,
	const Block_Info &block_info, const El::Grid &grid, const SDP_Solver_Parameters &parameters)
{
	//set_stream_precision(std::cout);
	//if (El::mpi::Rank() == 0) std::cout << "sdp0.objective_const = " << sdp0.objective_const << "\n";
	if(parameters.sdpd_mode_exact_dSDP)exit(0);
	if (parameters.list_sdp2_path.size() < 3) { std::cout << "parameters.list_sdp2_path.size() < 3 \n"; exit(0); }

	std::vector<SDP*> sdp_10_01_11;

	sdp_10_01_11.push_back(new SDP(parameters.list_sdp2_path[0], block_info, grid));
	sdp_10_01_11.push_back(new SDP(parameters.list_sdp2_path[1], block_info, grid));
	sdp_10_01_11.push_back(new SDP(parameters.list_sdp2_path[2], block_info, grid));

	dsdp_i = new SDP(sdp00);
	sdp_substract(*dsdp_i, *sdp_10_01_11[0], block_info);
	sdp_multiply(*dsdp_i, -1 / parameters.finiteDifference, block_info);

	dsdp_j = new SDP(sdp00);
	sdp_substract(*dsdp_j, *sdp_10_01_11[1], block_info);
	sdp_multiply(*dsdp_j, -1 / parameters.finiteDifference, block_info);

	ddsdp_ij = new SDP(sdp00);
	sdp_substract(*ddsdp_ij, *sdp_10_01_11[0], block_info);
	sdp_substract(*ddsdp_ij, *sdp_10_01_11[1], block_info);
	sdp_add(*ddsdp_ij, *sdp_10_01_11[2], block_info);
	sdp_multiply(*ddsdp_ij, 1 / (parameters.finiteDifference*parameters.finiteDifference), block_info);

	for (int i = 0; i < sdp_10_01_11.size(); i++) delete sdp_10_01_11[i];
	return 1;
}



// In this case, the -d file should be f[x+e], f[x-e] , where f is sdp data.
// So dsdp_list should be {(f[x] - f[x+e])/e, (f[x] - f[x-e])/e}
int initialize_ddsdp(SDP & sdp0, SDP* & dsdp, SDP* & ddsdp, const Block_Info &block_info, const El::Grid &grid, const SDP_Solver_Parameters &parameters)
{
	if (El::mpi::Rank() == 0) std::cout << "sdpd_mode_exact_dSDP:" << parameters.sdpd_mode_exact_dSDP << "\n";

	if (!parameters.sdpd_mode_exact_dSDP)
	{
		std::vector<SDP*> dsdp_list;

		initialize_dsdp(sdp0, dsdp_list, block_info, grid, parameters);

		//set_stream_precision(std::cout);
		//if (El::mpi::Rank() == 0) std::cout << "sdp0.objective_const = " << sdp0.objective_const << "\n";
		//if (El::mpi::Rank() == 0) std::cout << "dsdp[0].objective_const = " << dsdp_list[0]->objective_const << "\n";
		//if (El::mpi::Rank() == 0) std::cout << "dsdp[1].objective_const = " << dsdp_list[1]->objective_const << "\n";


		dsdp = new SDP(*dsdp_list[1]);
		sdp_substract(*dsdp, *dsdp_list[0], block_info);
		sdp_multiply(*dsdp, El::BigFloat(0.5), block_info);

		ddsdp = new SDP(*dsdp_list[0]);
		sdp_add(*ddsdp, *dsdp_list[1], block_info);
		sdp_multiply(*ddsdp, -1/parameters.finiteDifference, block_info);

		for (int i = 0; i < dsdp_list.size(); i++)
			delete dsdp_list[i];
	}
	else
	{
		if (parameters.sdpd_mode_dSDP_only_b)
		{
			dsdp = new SDP(sdp0);
			ddsdp = new SDP(sdp0);

			sdp_zero(*dsdp, block_info);
			sdp_zero(*ddsdp, block_info);

			SDP*dsdp_temp = new SDP(parameters.list_sdp2_path[0], block_info, grid, true);
			SDP*ddsdp_temp = new SDP(parameters.list_sdp2_path[1], block_info, grid, true);

			dsdp->dual_objective_b = dsdp_temp->dual_objective_b;
			ddsdp->dual_objective_b = ddsdp_temp->dual_objective_b;

			dsdp->objective_const = dsdp_temp->objective_const;
			ddsdp->objective_const = ddsdp_temp->objective_const;

			delete dsdp_temp;
			delete ddsdp_temp;
		}
		else
		{
			dsdp = new SDP(parameters.list_sdp2_path[0], block_info, grid);
			ddsdp = new SDP(parameters.list_sdp2_path[1], block_info, grid);
		}
	}

	return 1;
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

		// By += FreeVarMatrix * y
		Gemm(El::Orientation::NORMAL, El::Orientation::NORMAL, El::BigFloat(1),
			*free_var_matrix_block, *y_block, El::BigFloat(1),
			*By_block);

		++y_block;
		++free_var_matrix_block;
		++By_block;
	}

	return dot(x, By);
}

El::BigFloat compute_xdAY_part1(const Block_Info &block_info, const SDP &sdp, const SDP &dsdp,
	const Block_Vector &x, const Block_Diagonal_Matrix &Y)
{
	Block_Vector AY = x;

	auto dx_block(AY.blocks.begin());

	auto Z_block(Y.blocks.begin());
	auto bilinear_bases_block(sdp.bilinear_bases_dist.begin());
	auto d_bilinear_bases_block(dsdp.bilinear_bases_dist.begin());

	for (auto &block_index : block_info.block_indices)
	{
		// dx = 0
		*dx_block *= 0;
		const size_t dx_block_size(block_info.degrees[block_index] + 1);

		// dx[p] -= Tr(A_p Z)
		// Not sure whether it is better to first loop over blocks in
		// the result or over sub-blocks in Z
		for (size_t parity = 0; parity < 2; ++parity)
		{
			const size_t Z_block_size(bilinear_bases_block->Height());
			El::DistMatrix<El::BigFloat> ones(Z_block->Grid());
			El::Ones(ones, Z_block_size, 1);

			for (size_t column_block = 0;
				column_block < block_info.dimensions[block_index];
				++column_block)
				for (size_t row_block = 0; row_block <= column_block; ++row_block)
				{
					size_t column_offset(column_block * Z_block_size),
						row_offset(row_block * Z_block_size);

					El::DistMatrix<El::BigFloat> Z_sub_block(
						El::LockedView(*Z_block, row_offset, column_offset,
							Z_block_size, Z_block_size)),
						Z_times_q(Z_block_size, dx_block_size, Z_block->Grid());
					El::Zero(Z_times_q);
					El::DistMatrix<El::BigFloat> q_Z_q(Z_times_q);
					El::Zero(q_Z_q);

					El::Gemm(El::Orientation::NORMAL, El::Orientation::NORMAL,
						El::BigFloat(1), Z_sub_block, *bilinear_bases_block,
						El::BigFloat(0), Z_times_q);

					El::Hadamard(Z_times_q, *d_bilinear_bases_block, q_Z_q);

					const size_t dx_row_offset(
						((column_block * (column_block + 1)) / 2 + row_block)
						* dx_block_size);
					El::DistMatrix<El::BigFloat> dx_sub_block(
						El::View(*dx_block, dx_row_offset, 0, dx_block_size, 1));

					El::Gemv(El::Orientation::TRANSPOSE, El::BigFloat(-1), q_Z_q,
						ones, El::BigFloat(1), dx_sub_block);
				}
			++Z_block;
			++bilinear_bases_block;
			++d_bilinear_bases_block;
		}
		++dx_block;
	}

	return -dot(x, AY);
}

El::BigFloat compute_xdAY_part2(const Block_Info &block_info, const SDP &sdp, const SDP &dsdp,
	const Block_Vector &x, const Block_Diagonal_Matrix &Y)
{
	Block_Vector AY = x;

	auto dx_block(AY.blocks.begin());

	auto Z_block(Y.blocks.begin());
	auto bilinear_bases_block(sdp.bilinear_bases_dist.begin());
	auto d_bilinear_bases_block(dsdp.bilinear_bases_dist.begin());

	for (auto &block_index : block_info.block_indices)
	{
		// dx = 0
		*dx_block *= 0;
		const size_t dx_block_size(block_info.degrees[block_index] + 1);

		// dx[p] -= Tr(A_p Z)
		// Not sure whether it is better to first loop over blocks in
		// the result or over sub-blocks in Z
		for (size_t parity = 0; parity < 2; ++parity)
		{
			const size_t Z_block_size(bilinear_bases_block->Height());
			El::DistMatrix<El::BigFloat> ones(Z_block->Grid());
			El::Ones(ones, Z_block_size, 1);

			for (size_t column_block = 0;
				column_block < block_info.dimensions[block_index];
				++column_block)
				for (size_t row_block = 0; row_block <= column_block; ++row_block)
				{
					size_t column_offset(column_block * Z_block_size),
						row_offset(row_block * Z_block_size);

					El::DistMatrix<El::BigFloat> Z_sub_block(
						El::LockedView(*Z_block, row_offset, column_offset,
							Z_block_size, Z_block_size)),
						Z_times_q(Z_block_size, dx_block_size, Z_block->Grid());
					El::Zero(Z_times_q);
					El::DistMatrix<El::BigFloat> q_Z_q(Z_times_q);
					El::Zero(q_Z_q);

					El::Gemm(El::Orientation::NORMAL, El::Orientation::NORMAL,
						El::BigFloat(1), Z_sub_block, *d_bilinear_bases_block,
						El::BigFloat(0), Z_times_q);

					El::Hadamard(Z_times_q, *bilinear_bases_block, q_Z_q);

					const size_t dx_row_offset(
						((column_block * (column_block + 1)) / 2 + row_block)
						* dx_block_size);
					El::DistMatrix<El::BigFloat> dx_sub_block(
						El::View(*dx_block, dx_row_offset, 0, dx_block_size, 1));

					El::Gemv(El::Orientation::TRANSPOSE, El::BigFloat(-1), q_Z_q,
						ones, El::BigFloat(1), dx_sub_block);
				}
			++Z_block;
			++bilinear_bases_block;
			++d_bilinear_bases_block;
		}
		++dx_block;
	}

	return -dot(x, AY);
}

El::BigFloat compute_hessian_component(const SDP &dsdp, const SDP &ddsdp, const DSDPSOLUTION & sol, const SDP_Solver & solver, const Block_Info &block_info)
{

	El::BigFloat db_dy;
	if (!sol.dy.blocks.empty())
	{
		db_dy = El::Dotu(dsdp.dual_objective_b, sol.dy.blocks.front());
	}
	El::BigFloat dc_dx = dot(dsdp.primal_objective_c, sol.dx);
	El::BigFloat dx_dB_y = compute_xBy(block_info, dsdp, sol.dx, solver.y);
	El::BigFloat x_dB_dy = compute_xBy(block_info, dsdp, solver.x, sol.dy);

	El::BigFloat ddb_y;
	if (!sol.dy.blocks.empty())
	{
		ddb_y = El::Dotu(ddsdp.dual_objective_b, solver.y.blocks.front());
	}
	El::BigFloat ddc_x = dot(ddsdp.primal_objective_c, solver.x);
	El::BigFloat x_ddB_y = compute_xBy(block_info, ddsdp, solver.x, solver.y);

	El::BigFloat rslt_1st = db_dy + dc_dx - dx_dB_y - x_dB_dy;

	El::BigFloat rslt_2nd = ddb_y + ddc_x - x_ddB_y;

	return rslt_1st + rslt_2nd + ddsdp.objective_const;
}


El::BigFloat compute_hessian_component_off_diagonal(const SDP &dsdp_i, const SDP &dsdp_j,
	const SDP &ddsdp_ij, const DSDPSOLUTION & sol_i, const DSDPSOLUTION & sol_j,
	const SDP_Solver & solver, const Block_Info &block_info)
{

	El::BigFloat db_dy;
	if (!sol_i.dy.blocks.empty())
	{
		db_dy = El::Dotu(dsdp_i.dual_objective_b, sol_j.dy.blocks.front());
	}
	El::BigFloat dc_dx = dot(dsdp_i.primal_objective_c, sol_j.dx);
	El::BigFloat dx_dB_y = compute_xBy(block_info, dsdp_i, sol_j.dx, solver.y);
	El::BigFloat x_dB_dy = compute_xBy(block_info, dsdp_i, solver.x, sol_j.dy);

	El::BigFloat ddb_y;
	if (!sol_i.dy.blocks.empty())
	{
		ddb_y = El::Dotu(ddsdp_ij.dual_objective_b, solver.y.blocks.front());
	}
	El::BigFloat ddc_x = dot(ddsdp_ij.primal_objective_c, solver.x);
	El::BigFloat x_ddB_y = compute_xBy(block_info, ddsdp_ij, solver.x, solver.y);

	El::BigFloat rslt_1st = db_dy + dc_dx - dx_dB_y - x_dB_dy;
	El::BigFloat rslt_2nd = ddb_y + ddc_x - x_ddB_y;

	return rslt_1st + rslt_2nd + ddsdp_ij.objective_const;
}




El::BigFloat approxobj_db_dy, approxobj_dc_dx, approxobj_dx_dB_y, approxobj_x_dB_dy;

El::BigFloat compute_approx_objective(const SDP &dsdp, const DSDPSOLUTION & sol, const SDP_Solver & solver, const Block_Info &block_info)
{

	El::BigFloat db_dy;
	if (!sol.dy.blocks.empty())
	{
		db_dy = El::Dotu(dsdp.dual_objective_b, sol.dy.blocks.front());
	}
	El::BigFloat dc_dx = dot(dsdp.primal_objective_c, sol.dx);
	El::BigFloat dx_dB_y = compute_xBy(block_info, dsdp, sol.dx, solver.y);
	El::BigFloat x_dB_dy = compute_xBy(block_info, dsdp, solver.x, sol.dy);

	approxobj_db_dy = db_dy;
	approxobj_dc_dx = dc_dx;
	approxobj_dx_dB_y = dx_dB_y;
	approxobj_x_dB_dy = x_dB_dy;

	El::BigFloat rslt_1st = db_dy + dc_dx - dx_dB_y - x_dB_dy;

	return rslt_1st * El::BigFloat("0.5");
}


El::BigFloat Balt_part1, Balt_part2, Balt_part3, Balt_part4;

El::BigFloat compute_derivative_Balt_formula(const SDP & sdp, const SDP & dsdp, SDP_Solver & solver, const Block_Info &block_info)
{
	El::BigFloat dby;
	if (!solver.y.blocks.empty())
	{
		dby = dsdp.objective_const + El::Dotu(dsdp.dual_objective_b, solver.y.blocks.front());
	}

	El::BigFloat xdBy = compute_xBy(block_info, dsdp, solver.x, solver.y);

	El::BigFloat xdAY = compute_xdAY_part1(block_info, sdp, dsdp, solver.x, solver.Y)+ compute_xdAY_part2(block_info, sdp, dsdp, solver.x, solver.Y);

	El::BigFloat dprimalobj_Balt = dot(dsdp.primal_objective_c, solver.x) + dby - xdBy - xdAY;

	Balt_part1 = dot(dsdp.primal_objective_c, solver.x);
	Balt_part2 = dby;
	Balt_part3 = xdBy;
	Balt_part4 = xdAY;

	return dprimalobj_Balt;
}

int compute_gradient(std::vector<El::BigFloat> & dobj_list, const SDP & sdp, std::vector<SDP*> & dsdp_list, SDP_Solver & solver, const Block_Info &block_info)
{
	for (int i = 0; i < dsdp_list.size(); i++)
	{
		dobj_list.push_back(compute_derivative_Balt_formula(sdp, *dsdp_list[i], solver, block_info));
	}
	return 1;
}


/*
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
*/


void SDPDD_gradient(SDP &sdp, SDP_Solver &solver, El::Grid& grid, Timers &timers, const Block_Info &block_info, const SDP_Solver_Parameters &parameters)
{
	std::vector<SDP*> dsdp_list;
	initialize_dsdp(sdp, dsdp_list, block_info, grid, parameters);

	std::vector<El::BigFloat> dobj_list;

	compute_gradient(dobj_list, sdp, dsdp_list, solver, block_info);

	if (El::mpi::Rank() == 0)
	{
		set_stream_precision(std::cout);

		std::cout << "dc.x = " << Balt_part1 << "\n";
		std::cout << "db.y = " << Balt_part2 << "\n";
		std::cout << "x.dB.y = " << Balt_part3 << "\n";
		std::cout << "x.dA.Y = " << Balt_part4 << "\n";

		std::cout << "[SDPDReturnBegin.Gradient]\n";
		std::cout << "{\n";
		for (int i = 0; i < dobj_list.size(); i++)
		{
			std::cout << dobj_list[i] << "\n";
			if (i != dobj_list.size() - 1) std::cout << ",";
		}
		std::cout << "}\n";
		std::cout << "[SDPDReturnEnd.Gradient]\n";

		std::cout << "[SDPDReturnBegin]" << dobj_list[0] << "[SDPDReturnEnd]\n";
	}

	return;
}


void SDPDD_hessian_ij(SDP &sdp00, SDP_Solver &solver, El::Grid& grid, Timers &timers, const Block_Info &block_info, const SDP_Solver_Parameters &parameters)
{
	SDP*dsdp_i,*dsdp_j,*ddsdp_ij;
	initialize_ddsdp_ij(sdp00, dsdp_i, dsdp_j, ddsdp_ij, block_info, grid, parameters);

	std::vector<SDP*> dsdp_list;
	dsdp_list.push_back(dsdp_i);
	dsdp_list.push_back(dsdp_j);

	SDP_Solver_Terminate_Reason reason = solver.run(parameters, block_info, sdp00, dsdp_list, grid, timers);

	El::BigFloat dobj_i = compute_derivative_Balt_formula(sdp00, *dsdp_i, solver, block_info);
	El::BigFloat dobj_j = compute_derivative_Balt_formula(sdp00, *ddsdp_ij, solver, block_info);

	El::BigFloat ddobj_ij = compute_hessian_component_off_diagonal(*dsdp_i, *dsdp_j, *ddsdp_ij,
		*solver.dsdp_sol_list[0], *solver.dsdp_sol_list[1],
		solver, block_info);

	if (parameters.verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
	{
		set_stream_precision(std::cout);
		std::cout << "[SDPDReturnBegin.Gradient]\n";
		std::cout << "{" << dobj_i << "," << dobj_j << "}";
		std::cout << "[SDPDReturnEnd.Gradient]\n";

		std::cout << "[SDPDReturnBegin.Hessian]\n";
		std::cout << ddobj_ij;
		std::cout << "[SDPDReturnEnd.Hessian]\n";
	}

	return;
}


void SDPDD_hessian_ii(SDP &sdp, SDP_Solver &solver, El::Grid& grid, Timers &timers, const Block_Info &block_info, const SDP_Solver_Parameters &parameters)
{
	SDP*dsdp;
	SDP*ddsdp;

	initialize_ddsdp(sdp, dsdp, ddsdp, block_info, grid, parameters);

	//set_stream_precision(std::cout);
	//El::BigFloat obj_const_over_FD = dsdp->objective_const / parameters.finiteDifference;
	//if (El::mpi::Rank() == 0) std::cout << "finiteDifference = " << parameters.finiteDifference << "\n1/finiteDifference = " << 1 / parameters.finiteDifference << "\n";
	//if (El::mpi::Rank() == 0) std::cout << "dsdp.objective_const = " << dsdp->objective_const << "\n";
	//if (El::mpi::Rank() == 0) std::cout << "dsdp.objective_const/finiteDifference = " << obj_const_over_FD << "\n";
	//if (El::mpi::Rank() == 0) std::cout << "dsdp.objective_const/(dsdp.objective_const/finiteDifference) = " << dsdp->objective_const/obj_const_over_FD << "\n";
	//if (El::mpi::Rank() == 0) std::cout << "finiteDifference.prec = " << parameters.finiteDifference.Precision() << "\n";
	//if (El::mpi::Rank() == 0) std::cout << "dsdp.objective_const.prec = " << dsdp->objective_const.Precision() << "\n";
	//if (El::mpi::Rank() == 0) std::cout << "El::BigFloat(1 / 2)=" << El::BigFloat(0.5) << "\nEl::BigFloat(1 / 2).prec = " << El::BigFloat(0.5).Precision() << "\n";
	//if (El::mpi::Rank() == 0) std::cout << "ddsdp->objective_const = " << ddsdp->objective_const;

	std::vector<SDP*> dsdp_list;
	dsdp_list.push_back(dsdp);

	//if (El::mpi::Rank() == 0) std::cout << "initialize_ddsdp finished \n";

	SDP_Solver_Terminate_Reason reason
		= solver.run(parameters, block_info, sdp, dsdp_list, grid, timers);

	//if (El::mpi::Rank() == 0) std::cout << "dx,dX,dy,dY finished \n";

	El::BigFloat dobj = compute_derivative_Balt_formula(sdp, *dsdp, solver, block_info);
	El::BigFloat ddobj = compute_hessian_component(*dsdp, *ddsdp, *solver.dsdp_sol_list[0], solver, block_info);

	if (parameters.verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
	{
		//std::cout << "El::mpi::Rank() : " << El::mpi::Rank() << "\n";
		//std::cout << "parameters.verbosity : " << static_cast<int>(parameters.verbosity) << "\n";

		set_stream_precision(std::cout);
		std::cout << "[SDPDReturnBegin.C1Derivative]\n";
		std::cout << dobj;
		std::cout << "[SDPDReturnEnd.C1Derivative]\n";

		std::cout << "[SDPDReturnBegin.C2Derivative]\n";
		std::cout << ddobj;
		std::cout << "[SDPDReturnEnd.C2Derivative]\n";
	}

	return;
}

void SDPDD_approxobj(SDP &sdp, SDP_Solver &solver, El::Grid& grid, Timers &timers, const Block_Info &block_info, const SDP_Solver_Parameters &parameters)
{
	std::vector<SDP*> dsdp_list;
	initialize_dsdp(sdp, dsdp_list, block_info, grid, parameters);

	SDP_Solver_Terminate_Reason reason
		= solver.run(parameters, block_info, sdp, dsdp_list, grid, timers);

	El::BigFloat approxobj = compute_approx_objective(*dsdp_list[0], *solver.dsdp_sol_list[0], solver, block_info);

	if (El::mpi::Rank() == 0)
	{
		//rslt_1st = db_dy + dc_dx - dx_dB_y - x_dB_dy;

		set_stream_precision(std::cout);

		std::cout << "db_dy = " << approxobj_db_dy << "\n";
		std::cout << "dc_dx = " << approxobj_dc_dx << "\n";
		std::cout << "dx_dB_y = " << approxobj_dx_dB_y << "\n";
		std::cout << "x_dB_dy = " << approxobj_x_dB_dy << "\n";

		std::cout << "[SDPDReturnBegin.ApproxObjective]\n";
		std::cout << approxobj << "\n";
		std::cout << "[SDPDReturnEnd.ApproxObjective]\n";

		std::cout << "[SDPDReturnBegin]" << approxobj << "[SDPDReturnEnd]\n";
	}

	return;
}

Timers solve(const Block_Info &block_info, const SDP_Solver_Parameters &parameters)
{
  // Read an SDP from sdpFile and create a solver for it
  El::Grid grid(block_info.mpi_comm.value);

  Timers timers(parameters.verbosity >= Verbosity::debug);

  SDP sdp(parameters.sdp_directory, block_info, grid);

  SDP_Solver solver(parameters, block_info, grid,
                    sdp.dual_objective_b.Height());

  if (parameters.sdpd_mode_hessian_ii)
  {
	  SDPDD_hessian_ii(sdp, solver, grid, timers, block_info, parameters);
	  return timers;
  }

  if (parameters.sdpd_mode_hessian_ij)
  {
	  SDPDD_hessian_ij(sdp, solver, grid, timers, block_info, parameters);
	  return timers;
  }

  if (parameters.sdpd_mode_approx_objective)
  {
	  SDPDD_approxobj(sdp, solver, grid, timers, block_info, parameters);
	  return timers;
  }

  // default : gradient mode
  SDPDD_gradient(sdp, solver, grid, timers, block_info, parameters);
  if (!parameters.no_final_checkpoint)
  {
	  solver.save_solution(SDP_Solver_Terminate_Reason::PrimalDualOptimal, timers.front(), parameters.out_directory,
		  parameters.write_solution, 0,
		  block_info.block_indices,
		  parameters.verbosity);
  }

  return timers;
}







//////////////////// trash ///////////////////////




/*
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
*/

/*
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

std::cout << "[SDPDReturnBegin.Gradient.NingFormula.Primal]" << sdpd_cdx + sdpd_dcx << "[SDPDReturnEnd.Gradient.NingFormula.Primal]\n";
std::cout << "[SDPDReturnBegin.Gradient.NingFormula.Dual]" << sdpd_bdy + sdpd_dby << "[SDPDReturnEnd.Gradient.NingFormula.Dual]\n";

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
*/

/*
if(!parameters.no_final_checkpoint)
{
solver.save_checkpoint(parameters);
}*/


