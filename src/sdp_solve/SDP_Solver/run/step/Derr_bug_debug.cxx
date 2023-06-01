#include "../../../SDP_Solver.hxx"


// C := alpha*A*B + beta*C
void scale_multiply_add(const El::BigFloat &alpha,
	const Block_Diagonal_Matrix &A,
	const Block_Diagonal_Matrix &B,
	const El::BigFloat &beta, Block_Diagonal_Matrix &C);

// C := A B
inline void multiply(const Block_Diagonal_Matrix &A,
	const Block_Diagonal_Matrix &B, Block_Diagonal_Matrix &C)
{
	scale_multiply_add(El::BigFloat(1), A, B, El::BigFloat(0), C);
}

// X := ACholesky^{-T} ACholesky^{-1} X = A^{-1} X
void cholesky_solve(const Block_Diagonal_Matrix &ACholesky,
	Block_Diagonal_Matrix &X);

void compute_schur_RHS(const Block_Info &block_info, const SDP &sdp,
	const Block_Vector &dual_residues,
	const Block_Diagonal_Matrix &Z,
	Block_Vector &dx);

void solve_schur_complement_equation(
	const Block_Diagonal_Matrix &schur_complement_cholesky,
	const Block_Matrix &schur_off_diagonal,
	const El::DistMatrix<El::BigFloat> &Q, Block_Vector &dx, Block_Vector &dy);

El::BigFloat dot(const Block_Vector &A, const Block_Vector &B);


// this function is based on compute_search_direction
void compute_minus_d_minus_TrApZ(
	const Block_Info &block_info, const SDP &sdp, const SDP_Solver &solver,
	const Block_Diagonal_Matrix &X_cholesky, const El::BigFloat beta,
	const El::BigFloat &mu,
	const bool &is_corrector_phase,
	Block_Vector &dx, const Block_Diagonal_Matrix &dX, const Block_Diagonal_Matrix &dY, Block_Diagonal_Matrix &Z)
{
	// R = beta mu I - X Y (predictor phase)
	// R = beta mu I - X Y - dX dY (corrector phase)
	Block_Diagonal_Matrix R(solver.X);

	scale_multiply_add(El::BigFloat(-1), solver.X, solver.Y, El::BigFloat(0), R);
	if (is_corrector_phase)
	{
		scale_multiply_add(El::BigFloat(-1), dX, dY, El::BigFloat(1), R);
	}
	R.add_diagonal(beta * mu);

	// Z = Symmetrize(X^{-1} (PrimalResidues Y - R))
	multiply(solver.primal_residues, solver.Y, Z);
	Z -= R;
	cholesky_solve(X_cholesky, Z);
	Z.symmetrize();

	// dx[p] = -dual_residues[p] - Tr(A_p Z)
	// dy[n] = dualObjective[n] - (FreeVarMatrix^T x)_n
	compute_schur_RHS(block_info, sdp, solver.dual_residues, Z, dx);
}



//////////////// compute tr(A Z) //////////////////////////
// this function is based on compute_schur_RHS

void compute_trAp_Z(const Block_Info &block_info, const SDP &sdp,
	const Block_Diagonal_Matrix &Z,
	Block_Vector &result)
{
	auto dx_block(result.blocks.begin());

	Block_Diagonal_Matrix Z_copy(Z);
	Z_copy.symmetrize();

	auto Z_block(Z_copy.blocks.begin());

	auto bilinear_bases_block(sdp.bilinear_bases.begin());

	for (auto &block_index : block_info.block_indices)
	{
		// dx = 0
		*dx_block *= 0;
		const size_t dx_block_size(block_info.num_points[block_index]);

		// dx[p] = Tr(A_p Z)
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

					El::Hadamard(Z_times_q, *bilinear_bases_block, q_Z_q);

					const size_t dx_row_offset(
						((column_block * (column_block + 1)) / 2 + row_block)
						* dx_block_size);
					El::DistMatrix<El::BigFloat> dx_sub_block(
						El::View(*dx_block, dx_row_offset, 0, dx_block_size, 1));

					El::Gemv(El::Orientation::TRANSPOSE, El::BigFloat(1), q_Z_q,
						ones, El::BigFloat(0), dx_sub_block);
				}
			++Z_block;
			++bilinear_bases_block;
		}
		++dx_block;
	}
}

/////////////////////////////////////////////////////////////////////////////

// this function is based on solve_schur_complement_equation
El::BigFloat Block_Vector_sum_dot(const Block_Vector &dy)
{
	//Block_Vector dy(dy_org);

	int64_t height = dy.blocks.front().LocalHeight();

	El::DistMatrix<El::BigFloat> dy_dist;
	Zeros(dy_dist, height, 1);
	{
		El::Matrix<El::BigFloat> dy_sum;
		Zeros(dy_sum, height, 1);

		for (size_t block = 0; block < dy.blocks.size(); ++block)
		{
			// Locally sum contributions to dy
			for (int64_t row = 0; row < dy.blocks[block].LocalHeight(); ++row)
			{
				int64_t global_row(dy.blocks[block].GlobalRow(row));
				for (int64_t column = 0; column < dy.blocks[block].LocalWidth();
					++column)
				{
					int64_t global_column(dy.blocks[block].GlobalCol(column));
					dy_sum(global_row, global_column)
						+= dy.blocks[block].GetLocal(row, column);
				}
			}
		}

		// Send out updates for dy
		El::BigFloat zero(0);
		for (int64_t row = 0; row < dy_sum.Height(); ++row)
			for (int64_t column = 0; column < dy_sum.Width(); ++column)
			{
				if (dy_sum(row, column) != zero)
				{
					dy_dist.QueueUpdate(row, column, dy_sum(row, column));
				}
			}
	}
	dy_dist.ProcessQueues();

	El::BigFloat local_primal_error(0);
	for (int64_t row = 0; row < dy_dist.LocalHeight(); ++row)
		for (int64_t column = 0; column < dy_dist.LocalWidth();
			++column)
	{
		local_primal_error += dy_dist.GetLocal(row, column)*dy_dist.GetLocal(row, column);
	}

	El::BigFloat local_sum = Dotu(dy_dist, dy_dist);

	// somehow maxvalue value are zero for some ranks (possible ranks with duplicate block), but local_sum are non-zero for those ranks
	// std::cout << "rank" << El::mpi::Rank() << " local_sum=" << local_sum << ", maxvalue=" << local_primal_error << "\n";

	El::BigFloat maxvalue = El::mpi::AllReduce(local_primal_error, El::mpi::SUM, El::mpi::COMM_WORLD);

	//El::BigFloat p2 = El::mpi::AllReduce(local_sum, El::mpi::COMM_WORLD);
	//if (El::mpi::Rank() == 0)std::cout << "primal_error^2 = " << p2 << "\n";

	return maxvalue; // it seems the value is the same for all ranks.
}

// this function is based on solve_schur_complement_equation
El::BigFloat Block_Vector_sum_max(const Block_Vector &dy)
{
	//Block_Vector dy(dy_org);

	int64_t height = dy.blocks.front().LocalHeight();

	El::DistMatrix<El::BigFloat> dy_dist;
	Zeros(dy_dist, height, 1);
	{
		El::Matrix<El::BigFloat> dy_sum;
		Zeros(dy_sum, height, 1);

		for (size_t block = 0; block < dy.blocks.size(); ++block)
		{
			// Locally sum contributions to dy
			for (int64_t row = 0; row < dy.blocks[block].LocalHeight(); ++row)
			{
				int64_t global_row(dy.blocks[block].GlobalRow(row));
				for (int64_t column = 0; column < dy.blocks[block].LocalWidth();
					++column)
				{
					int64_t global_column(dy.blocks[block].GlobalCol(column));
					dy_sum(global_row, global_column)
						+= dy.blocks[block].GetLocal(row, column);
				}
			}
		}

		// Send out updates for dy
		El::BigFloat zero(0);
		for (int64_t row = 0; row < dy_sum.Height(); ++row)
			for (int64_t column = 0; column < dy_sum.Width(); ++column)
			{
				if (dy_sum(row, column) != zero)
				{
					dy_dist.QueueUpdate(row, column, dy_sum(row, column));
				}
			}
	}
	dy_dist.ProcessQueues();


	// Get the max error.
	El::BigFloat local_primal_error(0);
	for (int64_t row = 0; row < dy_dist.LocalHeight(); ++row)
		for (int64_t column = 0; column < dy_dist.LocalWidth(); ++column)
		{
			local_primal_error = std::max(local_primal_error, El::Abs(dy_dist.GetLocal(row, column)));
		}

	El::BigFloat maxvalue = El::mpi::AllReduce(local_primal_error, El::mpi::MAX, El::mpi::COMM_WORLD);

	return maxvalue; // it seems the value is the same for all ranks.
}



void Assert_DistMatrix_Local(const El::DistMatrix<El::BigFloat> & mat, auto message)
{
	if (mat.Height() != mat.LocalHeight() || mat.Height() != mat.LocalHeight())
	{
		std::cout << "Assert_DistMatrix_local failture : DistMatrix is not local. \n";
		std::cout << message;
		El::mpi::Abort(El::mpi::COMM_WORLD, 1);
	}
}

// experimental fact : X.blocks.size()=Y.blocks.size()=2*x.blocks.size()=2*y.blocks.size()
// !x.blocks.empty() && x.blocks.at(0).Grid().Rank()  are the same as   !X.blocks.empty() && X.blocks.at(0).Grid().Rank()
void Block_Diagonal_Matrix_Visit(const Block_Info &block_info, Block_Diagonal_Matrix &mat, auto operation)
{
	bool fakeblockQ = !mat.blocks.empty() && mat.blocks.at(0).Grid().Rank() != 0;
	if (fakeblockQ) return;

	for (size_t block = 0; block != block_info.block_indices.size(); ++block)
	{
		size_t block_index(block_info.block_indices.at(block));
		for (size_t psd_block(0); psd_block < 2; ++psd_block)
		{
			int globalID = 2 * block_index + psd_block;
			int localID = 2 * block + psd_block;

			auto &curblock = mat.blocks.at(localID);

			Assert_DistMatrix_Local(curblock, "Block_Diagonal_Matrix_Visit failture.\n");
			operation(globalID, curblock);
		}
	}
}

// Visit function to avoid overcounting the duplicate blocks
void Block_Vector_Visit(const Block_Info &block_info, auto &mat, auto operation)
{
	bool fakeblockQ = !mat.blocks.empty() && mat.blocks.at(0).Grid().Rank() != 0;
	if (fakeblockQ) return;

	for (size_t localID = 0; localID != block_info.block_indices.size(); ++localID)
	{
		size_t globalID(block_info.block_indices.at(localID));

		auto &curblock = mat.blocks.at(localID);

		Assert_DistMatrix_Local(curblock, "Block_Vector_Visit failture.\n");
		operation(globalID, curblock);
	}
}

// I believe this version is correct
El::BigFloat Block_Vector_p_max_V3(const Block_Info &block_info, Block_Vector &vec)
{
	El::BigFloat local_max(0);

	Block_Vector_Visit(block_info, vec,
		[&local_max](size_t globalID, El::DistMatrix<El::BigFloat> & mat) {

		int height = mat.Height();
		int width = mat.Width();

		if (width == 1)
		{
			for (int i = 0; i < height; i++)
				local_max = El::Max(local_max, El::Abs(mat.GetLocal(i, 0)));
		}
		else if (height == 1)
		{
			for (int i = 0; i < width; i++)
				local_max = El::Max(local_max, El::Abs(mat.GetLocal(0, i)));
		}
		else
		{
			if (El::mpi::Rank() == 0)std::cout << "compute_Block_Vector_norm2 error : input is not a vector" << "\n";
			exit(0);
		}
	});

	El::BigFloat global_max = El::mpi::AllReduce(local_max, El::mpi::MAX, El::mpi::COMM_WORLD);

	return global_max;
}

El::BigFloat Block_Vector_Square(const Block_Info &block_info, const Block_Vector &vec)
{
	El::BigFloat square(0);

	Block_Vector_Visit(block_info, vec,
		[&square](size_t globalID, const El::DistMatrix<El::BigFloat> & mat) {

		int height = mat.Height();
		int width = mat.Width();

		if (width == 1)
		{
			for (int i = 0; i < height; i++)
				square += mat.GetLocal(i, 0)*mat.GetLocal(i, 0);
		}
		else if (height == 1)
		{
			for (int i = 0; i < width; i++)
				square += mat.GetLocal(0, i)*mat.GetLocal(0, i);
		}
		else
		{
			if (El::mpi::Rank() == 0)std::cout << "compute_Block_Vector_norm2 error : input is not a vector" << "\n";
			exit(0);
		}
	});

	El::BigFloat global_square = El::mpi::AllReduce(square, El::mpi::SUM, El::mpi::COMM_WORLD);
	return global_square;
}


//////////////// compute tr(A Y) //////////////////////////
// this function is based on compute_dual_residues_and_error
void compute_trA_Y_V2(
	const Block_Info &block_info, const SDP &sdp,
	const std::array<
	std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
	&A_Y,
	Block_Vector &result)
{
	auto dual_residues_block(result.blocks.begin());

	size_t Q_index(0);
	El::BigFloat local_max(0);
	for (auto &block_index : block_info.block_indices)
	{
		El::Zero(*dual_residues_block);
		const size_t block_size(block_info.num_points[block_index]),
			dim(block_info.dimensions[block_index]);

		for (auto &A_Y_parity : A_Y)
		{
			for (size_t column_block = 0; column_block < dim; ++column_block)
				for (size_t row_block = 0; row_block <= column_block; ++row_block)
				{
					El::DistMatrix<El::BigFloat> lower_diagonal(El::GetDiagonal(
						A_Y_parity[Q_index][column_block][row_block]));

					size_t residue_row_offset(
						((column_block * (column_block + 1)) / 2 + row_block)
						* block_size);

					El::DistMatrix<El::BigFloat> residue_sub_block(El::View(
						*dual_residues_block, residue_row_offset, 0, block_size, 1));
					El::Axpy(El::BigFloat(1.0), lower_diagonal,
						residue_sub_block);
				}
		}

		local_max = El::Max(local_max, El::MaxAbs(*dual_residues_block));

		++dual_residues_block;
		++Q_index;
	}
	//dual_error = El::mpi::AllReduce(local_max, El::mpi::MAX, El::mpi::COMM_WORLD);
}


//////////////// compute B.y //////////////////////////////

// this function is based on compute_dual_residues_and_error
void compute_By(
	const Block_Info &block_info, const SDP &sdp, const Block_Vector &y, Block_Vector &result)
{
	auto dual_residues_block(result.blocks.begin());
	auto y_block(y.blocks.begin());
	auto free_var_matrix_block(sdp.free_var_matrix.blocks.begin());

	El::BigFloat local_max(0);
	for (auto &block_index : block_info.block_indices)
	{
		El::Zero(*dual_residues_block);

		// dualResidues = B * y
		// TODO: Shouldn't this be Gemv since y is a vector?
		El::Gemm(El::Orientation::NORMAL, El::Orientation::NORMAL,
			El::BigFloat(1), *free_var_matrix_block, *y_block,
			El::BigFloat(0), *dual_residues_block);

		local_max = El::Max(local_max, El::MaxAbs(*dual_residues_block));

		++y_block;
		++free_var_matrix_block;
		++dual_residues_block;
	}
	El::BigFloat global_max = El::mpi::AllReduce(local_max, El::mpi::MAX, El::mpi::COMM_WORLD);
}

//////////////// Block_Vector_p_max //////////////////////////

// this function has the problem of overcounting duplicate blocks
El::BigFloat Block_Vector_p_max(const Block_Info &block_info, const Block_Vector &Vp)
{
	auto Vp_block(Vp.blocks.begin());

	El::BigFloat local_max(0);
	for (auto &block_index : block_info.block_indices)
	{
		local_max = El::Max(local_max, El::MaxAbs(*Vp_block));
		++Vp_block;
	}
	El::BigFloat global_max = El::mpi::AllReduce(local_max, El::mpi::MAX, El::mpi::COMM_WORLD);
	return global_max;
}

// this function has the problem of overcounting duplicate blocks
El::BigFloat Block_Vector_p_max_V2(const Block_Info &block_info, const Block_Vector &vec)
{
	El::BigFloat global_max(0);

	auto vec_block(vec.blocks.begin());
	for (auto &block_index : block_info.block_indices)
	{
		int blockID = vec_block - vec.blocks.begin();

		int height = vec_block->Height();
		int width = vec_block->Width();

		if (width == 1)
		{
			for (int i = 0; i < height; i++)
				global_max = El::Max(global_max, El::Abs(vec_block->Get(i, 0)));
		}
		else if (height == 1)
		{
			for (int i = 0; i < width; i++)
				global_max = El::Max(global_max, El::Abs(vec_block->Get(0, i)));
		}
		else
		{
			if (El::mpi::Rank() == 0)std::cout << "compute_Block_Vector_norm2 error : input is not a vector" << "\n";
			exit(0);
		}

		++vec_block;
	}
	return global_max;
}

void Block_Vector_plus_equal(Block_Vector &B, const Block_Vector &A)
{
	for (size_t block = 0; block < B.blocks.size(); ++block)
		El::Axpy(El::BigFloat(1), A.blocks[block], B.blocks[block]);
}

void Block_Vector_minus_equal(Block_Vector &B, const Block_Vector &A)
{
	for (size_t block = 0; block < B.blocks.size(); ++block)
		El::Axpy(El::BigFloat(-1), A.blocks[block], B.blocks[block]);
}


//////////////// compute S.x //////////////////////////////

void compute_Sx(
	const Block_Info &block_info, Block_Diagonal_Matrix &schur_complement, const Block_Vector &x, Block_Vector &result)
{
	auto x_block(x.blocks.begin());
	auto schur_complement_block(schur_complement.blocks.begin());
	auto result_block(result.blocks.begin());

	for (auto &block_index : block_info.block_indices)
	{
		El::Gemv(El::Orientation::NORMAL, El::BigFloat(1),
			*schur_complement_block, *x_block, El::BigFloat(0),
			*result_block);

		++schur_complement_block;
		++x_block;
		++result_block;
	}
}

