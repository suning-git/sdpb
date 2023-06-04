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

El::BigFloat min_eigenvalue(Block_Diagonal_Matrix &A);

El::BigFloat min_eigenvalue_safe(const Block_Diagonal_Matrix &A)
{
	Block_Diagonal_Matrix A_copy(A);
	return min_eigenvalue(A_copy);
}

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
// this function compute x_p=tr(A_p M) for a generic matrix M
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
		const size_t dx_block_size(block_info.num_points[block_index]);   //  dx_block_size=dim(k)=d_j + 1

		// dx[p] = Tr(A_p Z)
		// Not sure whether it is better to first loop over blocks in
		// the result or over sub-blocks in Z
		for (size_t parity = 0; parity < 2; ++parity)
		{
			const size_t Z_block_size(bilinear_bases_block->Height());  // delta_{j,1} +1  or  delta_{j,2} +1
			El::DistMatrix<El::BigFloat> ones(Z_block->Grid());
			El::Ones(ones, Z_block_size, 1);

			for (size_t column_block = 0;   
				column_block < block_info.dimensions[block_index];   // go through 0, ... , m_j-1
				++column_block)  
				for (size_t row_block = 0; row_block <= column_block; ++row_block)
				{
					size_t column_offset(column_block * Z_block_size),
						row_offset(row_block * Z_block_size);

					El::DistMatrix<El::BigFloat> Z_sub_block(
						El::LockedView(*Z_block, row_offset, column_offset,    //  Z_b_rs subblock
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
						ones, El::BigFloat(1), dx_sub_block);
				}
			++Z_block;
			++bilinear_bases_block;
		}
		++dx_block;
	}
}
/////////////////////////////////////////////////////////////////////////////

size_t Block_Vector_p_index_translate(const Block_Info &block_info, size_t j, size_t r, size_t s, size_t k)
{
	const size_t num_points(block_info.num_points[j]);

	if (r > s || k >= num_points)
	{
		std::cout << "Rank" << El::mpi::Rank() << " : invalide (j,r,s,k) index : (" 
			<< j << "," << r << "," << s << "," << k << ")" << ", num_points=" << num_points << "\n";
		exit(0);
	}

	return (s*(s + 1) / 2 + r)*num_points + k;
}

void Block_Vector_p_Set(const Block_Info &block_info, Block_Vector &x,
	size_t j, size_t r, size_t s, size_t k, const El::BigFloat &value)
{
	for (size_t localID = 0; localID != block_info.block_indices.size(); ++localID)
	{
		size_t globalID(block_info.block_indices.at(localID));

		auto &x_block = x.blocks.at(localID);

		bool fakeblockQ = x_block.Grid().Rank() != 0;

		size_t local_offset = 0;

		if (globalID == j)
		{
			local_offset = Block_Vector_p_index_translate(block_info, j, r, s, k);
			x_block.Set(local_offset, 0, value);
		}

		/*
		std::cout << "Rank=" << El::mpi::Rank()
			<< ", globalID=" << globalID
			<< ", global_index=" << "(" << x_block.GlobalRow(local_offset) << "," << x_block.GlobalCol(0) << ")"
			<< ", fakeblockQ=" << fakeblockQ
			<< ", grid.rank=" << x_block.Grid().Rank()
			<< ", local dim=" << "(" << x_block.LocalHeight() << "," << x_block.LocalWidth() << ")"
			<< ", global dim=" << "(" << x_block.Height() << "," << x_block.Width() << ")"
			<< ", num_points=" << block_info.num_points[globalID]
			<< "\n" << std::flush;
			*/

	}
}

void Block_Vector_p_Zero(const Block_Info &block_info, Block_Vector &x)
{
	for (size_t localID = 0; localID != block_info.block_indices.size(); ++localID)
	{
		auto &x_block = x.blocks.at(localID);
		El::Zero(x_block);
	}
}

void Block_Vector_p_Unit(const Block_Info &block_info, Block_Vector &x,
	size_t j, size_t r, size_t s, size_t k)
{
	Block_Vector_p_Zero(block_info, x);
	Block_Vector_p_Set(block_info, x, j, r, s, k, El::BigFloat(1));
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

	//if (El::mpi::Rank() == 0)std::cout << "test 1" << "\n";

	int64_t height = dy.blocks.front().LocalHeight();

	El::DistMatrix<El::BigFloat> dy_dist;
	Zeros(dy_dist, height, 1);
	{
		El::Matrix<El::BigFloat> dy_sum;
		Zeros(dy_sum, height, 1);

		//if (El::mpi::Rank() == 0)std::cout << "test 2" << "\n";

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

		//if (El::mpi::Rank() == 0)std::cout << "test 3" << "\n";

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

	//if (El::mpi::Rank() == 0)std::cout << "test 4" << "\n";

	// Get the max error.
	El::BigFloat local_primal_error(0);
	for (int64_t row = 0; row < dy_dist.LocalHeight(); ++row)
		for (int64_t column = 0; column < dy_dist.LocalWidth(); ++column)
		{
			local_primal_error = std::max(local_primal_error, El::Abs(dy_dist.GetLocal(row, column)));
		}

	//if (El::mpi::Rank() == 0)std::cout << "test 5" << "\n";

	El::BigFloat maxvalue = El::mpi::AllReduce(local_primal_error, El::mpi::MAX, El::mpi::COMM_WORLD);

	//if (El::mpi::Rank() == 0)std::cout << "test 6" << "\n";

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

// I believe this version is correct for D-err
El::BigFloat Block_Vector_p_max_V3(const Block_Info &block_info, Block_Vector &vec)
{
	El::BigFloat local_max(0);

	Block_Vector_Visit(block_info, vec,
		[&local_max](size_t globalID, El::DistMatrix<El::BigFloat> & mat) {

		int height = mat.LocalHeight();
		int width = mat.LocalWidth();

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

		int height = mat.LocalHeight();
		int width = mat.LocalWidth();

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
void compute_tr_AM_pair(
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


//////////////// compute tr(A Z) V2 //////////////////////////
// use A_Z pair
void compute_A_Y(
	const Block_Info &block_info, const Block_Diagonal_Matrix &Y,
	const std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
	std::array<std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>,
	2> &A_Y);


void compute_trAp_Z_V2(const Block_Info &block_info, const SDP &sdp,
	const Block_Diagonal_Matrix &Z,
	Block_Vector &result)
{
	std::array<std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2> A_Z;

	compute_A_Y(block_info, Z, sdp.bases_blocks, A_Z);

	compute_tr_AM_pair(block_info, sdp, A_Z, result);
	return;
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
	const Block_Info &block_info, const Block_Diagonal_Matrix &schur_complement, const Block_Vector &x, Block_Vector &result)
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

////////////////// print blocks information /////////////////////

void debug_print_Vector_blocks_info(const Block_Info &block_info, const std::vector<El::DistMatrix<El::BigFloat>> & blocks, const std::string var)
{
	bool fakeblockQ = !blocks.empty() && blocks.at(0).Grid().Rank() != 0;

	std::cout << std::flush;

	for (size_t localID = 0; localID != block_info.block_indices.size(); ++localID)
	{
		size_t globalID(block_info.block_indices.at(localID));

		const El::DistMatrix<El::BigFloat> & mat = blocks[localID];

		std::cout << var << "[Rank" << El::mpi::Rank() << "] " <<
			"localID=" << localID << ", local : (" << mat.LocalHeight() << "," << mat.LocalWidth() << ")"
			<< ", global : (" << mat.Height() << "," << mat.Width() << ")"
			<< ", fakeQ=" << fakeblockQ << ", grid.rank=" << mat.Grid().Rank()
			<< ", globalID="<< globalID << "\n" << std::flush;
	}
}

void debug_print_Matrix_blocks_info(const Block_Info &block_info, const std::vector<El::DistMatrix<El::BigFloat>> & blocks, const std::string var)
{
	bool fakeblockQ = !blocks.empty() && blocks.at(0).Grid().Rank() != 0;

	std::cout << std::flush;

	for (size_t block = 0; block != block_info.block_indices.size(); ++block)
	{
		size_t block_index(block_info.block_indices.at(block));
		for (size_t psd_block(0); psd_block < 2; ++psd_block)
		{
			int globalID = 2 * block_index + psd_block;
			int localID = 2 * block + psd_block;

			const El::DistMatrix<El::BigFloat> & mat = blocks[localID];

			std::cout << var << "[Rank" << El::mpi::Rank() << "] " <<
				"localID=" << localID << ", local : (" << mat.LocalHeight() << "," << mat.LocalWidth() << ")"
				<< ", global : (" << mat.Height() << "," << mat.Width() << ")"
				<< ", fakeQ=" << fakeblockQ << ", grid.rank=" << mat.Grid().Rank()
				<< ", globalID=" << globalID << "\n" << std::flush;
		}
	}
}

void debug_print_blocks_info(const Block_Info &block_info, const std::vector<El::DistMatrix<El::BigFloat>> & blocks, const std::string var)
{
	if (blocks.size() == block_info.block_indices.size())
		return debug_print_Vector_blocks_info(block_info, blocks, var);
	else if (blocks.size() == 2 * block_info.block_indices.size())
		return debug_print_Matrix_blocks_info(block_info, blocks, var);
	else
	{
		std::cout << var << "has size " << blocks.size() << ", while block_info has size " << block_info.block_indices.size() <<  "\n"
			<< "no available debug_print_blocks_info function for this structure." << std::flush;
	}
}


void debug_Ap_V1_vs_V2(const Block_Info &block_info, const std::vector<El::DistMatrix<El::BigFloat>> & blocks, const std::string var)
{
	bool fakeblockQ = !blocks.empty() && blocks.at(0).Grid().Rank() != 0;

	std::cout << std::flush;

	for (size_t localID = 0; localID != block_info.block_indices.size(); ++localID)
	{
		size_t globalID(block_info.block_indices.at(localID));

		const El::DistMatrix<El::BigFloat> & mat = blocks[localID];

		if (globalID == 1)
		{
			std::cout << var << " globalID=" << globalID << "\n" << std::flush;
			El::Print(mat);
		}
	}
}


// this function is based on solve_schur_complement_equation
El::BigFloat Block_Vector_sum_dot_debug(const Block_Vector &dy)
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
					std::ostringstream message;
					message.precision(100);

					message << "QueueUpdate[" << "Rank->" << El::mpi::Rank() << ", "
						<< " Entry[" << row << "," << column << "]->" << dy_sum(row, column) << "]\n";
					
					//std::cout << message.str() << std::flush;

					dy_dist.QueueUpdate(row, column, dy_sum(row, column));
				}
			}

		/*
		for (int rankID = 0; rankID < 3; rankID++)
		{
			El::mpi::Barrier(MPI_COMM_WORLD);
			usleep(10000);
			if (rankID == El::mpi::Rank())
			{
				std::cout << "dy_sum on rank" << El::mpi::Rank() << "\n" << std::flush;
				El::Print(dy_sum);
				std::cout << "\n" << std::flush;
			}
			usleep(100000);
			El::mpi::Barrier(MPI_COMM_WORLD);
		}
		dy_dist.QueueUpdate(El::mpi::Rank(), 0, El::BigFloat(El::mpi::Rank()));
		dy_dist.QueueUpdate(El::mpi::Rank(), 0, El::BigFloat(El::mpi::Rank()));
		*/
	}
	dy_dist.ProcessQueues();

	El::BigFloat local_sum_V1(0);
	for (int64_t row = 0; row < dy_dist.LocalHeight(); ++row)
		for (int64_t column = 0; column < dy_dist.LocalWidth();
			++column)
	{
		local_sum_V1 += dy_dist.GetLocal(row, column)*dy_dist.GetLocal(row, column);
	}

	El::BigFloat local_sum_V2 = Dotu(dy_dist, dy_dist);

	El::BigFloat globalsum = El::mpi::AllReduce(local_sum_V1, El::mpi::SUM, El::mpi::COMM_WORLD);

	// somehow maxvalue value are zero for some ranks (possible ranks with duplicate block), but local_sum are non-zero for those ranks

	//const El::DistMatrix<El::BigFloat> & mat = dy_dist;
	//bool fakeblockQ = mat.Grid().Rank() != 0;

	//El::Print(dy_dist, "dy_dist");

	/*
	std::cout << "p" << "[Rank" << El::mpi::Rank() << "] "
		<< ", local : (" << mat.LocalHeight() << "," << mat.LocalWidth() << ")"
		<< ", global : (" << mat.Height() << "," << mat.Width() << ")"
		<< ", fakeQ=" << fakeblockQ << ", grid.rank=" << mat.Grid().Rank()
		<< ", local_sum_V1=" << local_sum_V1
		<< ", local_sum_V2=" << local_sum_V2
		<< ", globalsum=" << globalsum 
		<< "\n" << std::flush;
		*/

	return globalsum; // it seems the value is the same for all ranks.
}

void constraint_matrix_weighted_sum(const Block_Info &block_info,
	const SDP &sdp, const Block_Vector &a,
	Block_Diagonal_Matrix &Result);

void compute_minus_InvX_Apdx_Y(const Block_Info &block_info, const SDP &sdp, const Block_Vector &dx,
	const Block_Diagonal_Matrix &X_cholesky, const Block_Diagonal_Matrix &Y, Block_Diagonal_Matrix &result)
{
	Block_Diagonal_Matrix Apdx(Y);
	constraint_matrix_weighted_sum(block_info, sdp, dx, Apdx);

	scale_multiply_add(El::BigFloat(-1), Apdx, Y, El::BigFloat(0), result);
	cholesky_solve(X_cholesky, result);
	return;
}



////////////////////// check tr(A_i X^-1 (A.dx) Y)  vs  S.dx //////////////////

void compute_tr_A_InvX_Adx_Y(const Block_Info &block_info, const SDP &sdp, const Block_Vector &dx,
	const Block_Diagonal_Matrix &X_cholesky, const Block_Diagonal_Matrix &Y, Block_Vector &result)
{
	Block_Diagonal_Matrix Apdx(Y);
	constraint_matrix_weighted_sum(block_info, sdp, dx, Apdx);

	Block_Diagonal_Matrix InvX_Adx_Y(Y);

	scale_multiply_add(El::BigFloat(1), Apdx, Y, El::BigFloat(0), InvX_Adx_Y);
	cholesky_solve(X_cholesky, InvX_Adx_Y);

	compute_trAp_Z(block_info, sdp, InvX_Adx_Y, result);
}


void compute_tr_A_InvX_Adx_Y_vs_Sdx(const Block_Info &block_info, const SDP &sdp, const Block_Vector &dx,
	const Block_Diagonal_Matrix &X_cholesky, const Block_Diagonal_Matrix &Y, const Block_Diagonal_Matrix &schur_complement, Block_Vector &result)
{
	compute_tr_A_InvX_Adx_Y(block_info, sdp, dx, X_cholesky, Y, result);

	Block_Vector Sdx(dx);
	compute_Sx(block_info, schur_complement, dx, Sdx);

	Block_Vector_minus_equal(result, Sdx);
}


void print_Matrix_block(const Block_Info &block_info, const Block_Diagonal_Matrix &mat, size_t globalID_to_print, const std::string file_suffix);



// X := ACholesky^{-T} ACholesky^{-1} X = A^{-1} X
void cholesky_solve_fix(const Block_Diagonal_Matrix &ACholesky,
	Block_Diagonal_Matrix &X)
{
	for (size_t b = 0; b < X.blocks.size(); b++)
	{
		El::Trsm(El::LeftOrRight::LEFT, El::UpperOrLowerNS::LOWER,
			El::Orientation::NORMAL, El::UnitOrNonUnit::NON_UNIT,
			El::BigFloat(1), ACholesky.blocks[b], X.blocks[b]);
	}
	for (size_t b = 0; b < X.blocks.size(); b++)
	{
		El::Trsm(El::LeftOrRight::LEFT, El::UpperOrLowerNS::LOWER,
			El::Orientation::ADJOINT, El::UnitOrNonUnit::NON_UNIT,
			El::BigFloat(1), ACholesky.blocks[b], X.blocks[b]);
	}
}



void check_A_InvX(const Block_Info &block_info, const SDP &sdp, 
	const std::array<
	std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
	&A_X_inv, const Block_Diagonal_Matrix &X_cholesky, const Block_Diagonal_Matrix &X, const Block_Vector &x)
{
	Block_Vector tr_A_InvX(x);
	compute_tr_AM_pair(block_info, sdp, A_X_inv, tr_A_InvX);

	Block_Diagonal_Matrix InvX(X);
	InvX *= 0;
	InvX.add_diagonal(1);
	Block_Diagonal_Matrix InvX_org(InvX);

	cholesky_solve(X_cholesky, InvX);

	Block_Vector tr_A_InvX_V2(x);
	compute_trAp_Z_V2(block_info, sdp, InvX, tr_A_InvX_V2);

	Block_Vector check_tr_A_InvX_V1_vs_V2(tr_A_InvX);
	Block_Vector_minus_equal(check_tr_A_InvX_V1_vs_V2, tr_A_InvX_V2);
	El::BigFloat check_tr_A_InvX_V1_vs_V2_max = Block_Vector_p_max_V3(block_info, check_tr_A_InvX_V1_vs_V2);
	if (El::mpi::Rank() == 0)std::cout << "[check A_InvX] : V1-V2=" << check_tr_A_InvX_V1_vs_V2_max << "\n";

	El::BigFloat XInvmin = min_eigenvalue_safe(InvX);
	if (El::mpi::Rank() == 0)std::cout << "[check A_InvX] : min(InvX)=" << XInvmin << "\n";

	Block_Diagonal_Matrix check_InvX_sym(InvX);
	check_InvX_sym.symmetrize();
	check_InvX_sym -= InvX;
	El::BigFloat check_InvX_sym_max = check_InvX_sym.max_abs_mpi();
	if (El::mpi::Rank() == 0)std::cout << "[check A_InvX] : check_InvX_sym_max=" << check_InvX_sym_max << "\n";

	/**/
	print_Matrix_block(block_info, InvX, 3, "pd35_InvX");
	print_Matrix_block(block_info, X, 3, "pd35_X");
	print_Matrix_block(block_info, X_cholesky, 3, "pd35_choleskyX");
	

}



void check_S_identity(const Block_Info &block_info, const SDP &sdp,
	const Block_Diagonal_Matrix &schur_complement, const Block_Diagonal_Matrix &X_cholesky,
	const Block_Diagonal_Matrix &X, const Block_Diagonal_Matrix &Y, 
	const Block_Vector &x, const Block_Vector &dx)
{
	Block_Vector vec1(x), vec2(x);
	Block_Vector_p_Unit(block_info, vec1, 3, 0, 0, 3);
	Block_Vector_p_Unit(block_info, vec2, 3, 0, 0, 4);

	El::BigFloat v1_v1, v1_v2, v2_v2;
	v1_v1 = dot(vec1, vec1);
	v1_v2 = dot(vec1, vec2);
	v2_v2 = dot(vec2, vec2);

	Block_Vector S_v1(x), S_v2(x);
	compute_tr_A_InvX_Adx_Y_vs_Sdx(block_info, sdp, vec1, X_cholesky, Y, schur_complement, S_v1);
	compute_tr_A_InvX_Adx_Y_vs_Sdx(block_info, sdp, vec2, X_cholesky, Y, schur_complement, S_v2);

	El::BigFloat S_11 = dot(vec1, S_v1);
	El::BigFloat S_12 = dot(vec1, S_v2);
	El::BigFloat S_21 = dot(vec2, S_v1);
	El::BigFloat S_22 = dot(vec2, S_v2);

	if (El::mpi::Rank() == 0)std::cout << "[S identity] : S11 V1-V2 =" << S_11 << "\n";
	if (El::mpi::Rank() == 0)std::cout << "[S identity] : S12 V1-V2 =" << S_12 << "\n";
	if (El::mpi::Rank() == 0)std::cout << "[S identity] : S21 V1-V2 =" << S_21 << "\n";
	if (El::mpi::Rank() == 0)std::cout << "[S identity] : S22 V1-V2 =" << S_22 << "\n";

	El::BigFloat SvsA_p1_max = Block_Vector_p_max_V3(block_info, S_v1);
	El::BigFloat SvsA_p2_max = Block_Vector_p_max_V3(block_info, S_v2);

	if (El::mpi::Rank() == 0)std::cout << "[S identity] : S V1-V2 at p1.max=" << SvsA_p1_max
		<< ", p2.max=" << SvsA_p2_max  << "\n";

}

void check_S_identity_V2(const Block_Info &block_info, const SDP &sdp,
	const Block_Diagonal_Matrix &schur_complement, const Block_Diagonal_Matrix &X_cholesky,
	const Block_Diagonal_Matrix &X, const Block_Diagonal_Matrix &Y,
	const Block_Vector &x, const Block_Vector &dx)
{
	Block_Vector vec1(x), vec2(x);
	Block_Vector_p_Unit(block_info, vec1, 3, 0, 0, 3);
	Block_Vector_p_Unit(block_info, vec2, 3, 0, 0, 4);

	El::BigFloat v1_v1, v1_v2, v2_v2;
	v1_v1 = dot(vec1, vec1);
	v1_v2 = dot(vec1, vec2);
	v2_v2 = dot(vec2, vec2);

	Block_Vector AInvXAY_v1(x), AInvXAY_v2(x), S_v1(x), S_v2(x);
	compute_tr_A_InvX_Adx_Y(block_info, sdp, vec1, X_cholesky, Y, AInvXAY_v1);
	compute_tr_A_InvX_Adx_Y(block_info, sdp, vec2, X_cholesky, Y, AInvXAY_v2);
	compute_Sx(block_info, schur_complement, vec1, S_v1);
	compute_Sx(block_info, schur_complement, vec2, S_v2);

	El::BigFloat S_11 = dot(vec1, S_v1);
	El::BigFloat S_12 = dot(vec1, S_v2);
	El::BigFloat S_21 = dot(vec2, S_v1);
	El::BigFloat S_22 = dot(vec2, S_v2);

	El::BigFloat AInvXAY_11 = dot(vec1, AInvXAY_v1);
	El::BigFloat AInvXAY_12 = dot(vec1, AInvXAY_v2);
	El::BigFloat AInvXAY_21 = dot(vec2, AInvXAY_v1);
	El::BigFloat AInvXAY_22 = dot(vec2, AInvXAY_v2);

	if (El::mpi::Rank() == 0)std::cout << 
		"[S identity] : S11=" << S_11 << ", S12=" << S_12 << ", S21=" << S_21 << ", S22=" << S_22
		<< ", S12-S21=" << (S_12 - S_21) << "\n";

	if (El::mpi::Rank() == 0)std::cout <<
		"[S identity] : A11=" << AInvXAY_11 << ", A12=" << AInvXAY_12 << ", A21=" << AInvXAY_21 << ", A22=" << AInvXAY_22 
		<< ", A12-A21=" << (AInvXAY_12 - AInvXAY_21) << "\n";

	if (El::mpi::Rank() == 0)std::cout <<
		"[S identity] : A12-S12=" << (S_12 - AInvXAY_12) << ", A21-S21=" << (AInvXAY_21 - S_21) << "\n";
}





/*

Block_Diagonal_Matrix InvXCholesky(X);
InvXCholesky *= 0;
InvXCholesky.add_diagonal(1);
for (size_t b = 0; b < X.blocks.size(); b++)
{
El::Trsm(El::LeftOrRight::LEFT, El::UpperOrLowerNS::LOWER,
El::Orientation::NORMAL, El::UnitOrNonUnit::NON_UNIT,
El::BigFloat(1), X_cholesky.blocks[b], InvXCholesky.blocks[b]);
}
print_Matrix_block(block_info, InvXCholesky, 3, "pd11_InvCholeskyX");

Block_Diagonal_Matrix InvTXCholesky_InvXCholesky(InvXCholesky);
for (size_t b = 0; b < X.blocks.size(); b++)
{
El::Trsm(El::LeftOrRight::LEFT, El::UpperOrLowerNS::LOWER,
El::Orientation::ADJOINT, El::UnitOrNonUnit::NON_UNIT,
El::BigFloat(1), X_cholesky.blocks[b], InvTXCholesky_InvXCholesky.blocks[b]);
}
print_Matrix_block(block_info, InvTXCholesky_InvXCholesky, 3, "pd11_InvTXCholesky_InvCholeskyX");
*/


void write_psd_block(const boost::filesystem::path &outfile,
	const El::DistMatrix<El::BigFloat> &block)
{
	boost::filesystem::ofstream stream;
	if (block.DistRank() == block.Root())
	{
		stream.open(outfile);
	}
	El::Print(block,
		std::to_string(block.Height()) + " "
		+ std::to_string(block.Width()),
		"\n", stream);
	if (block.DistRank() == block.Root())
	{
		stream << "\n";
		if (!stream.good())
		{
			throw std::runtime_error("Error when writing to: "
				+ outfile.string());
		}
	}
}

void print_Matrix_block(const Block_Info &block_info, const Block_Diagonal_Matrix &mat, size_t globalID_to_print, const std::string file_suffix)
{
	for (size_t block = 0; block != block_info.block_indices.size(); ++block)
	{
		size_t block_index(block_info.block_indices.at(block));
		for (size_t psd_block(0); psd_block < 2; ++psd_block)
		{
			int globalID = 2 * block_index + psd_block;
			int localID = 2 * block + psd_block;

			auto &curblock = mat.blocks.at(localID);

			if (globalID == globalID_to_print)
			{
				write_psd_block("./testdata/"+ file_suffix + "_" + std::to_string(globalID) + ".txt", curblock);
			}
		}
	}
}