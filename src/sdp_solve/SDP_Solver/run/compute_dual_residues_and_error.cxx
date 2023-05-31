#include "../../SDP_Solver.hxx"

// dual_residues[p] = c[p] - A[p,a,b] Y[a,b] - B[p,a] y[a]
//
// A[p,a,c] Y[c,b] = (1/2) A_Y[parity,r,s,p,a,b] + swap (r <-> s)

void compute_dual_residues_and_error(
  const Block_Info &block_info, const SDP &sdp, const Block_Vector &y,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_Y,
  Block_Vector &dual_residues, El::BigFloat &dual_error, Timers &timers)
{
  auto &dual_residues_timer(timers.add_and_start("run.computeDualResidues"));

  auto dual_residues_block(dual_residues.blocks.begin());
  auto primal_objective_c_block(sdp.primal_objective_c.blocks.begin());
  auto y_block(y.blocks.begin());
  auto free_var_matrix_block(sdp.free_var_matrix.blocks.begin());

  size_t Q_index(0);
  El::BigFloat local_max(0);
  for(auto &block_index : block_info.block_indices)
    {
      El::Zero(*dual_residues_block);
      const size_t block_size(block_info.num_points[block_index]),
        dim(block_info.dimensions[block_index]);

      for(auto &A_Y_parity : A_Y)
        {
          for(size_t column_block = 0; column_block < dim; ++column_block)
            for(size_t row_block = 0; row_block <= column_block; ++row_block)
              {
                El::DistMatrix<El::BigFloat> lower_diagonal(El::GetDiagonal(
                  A_Y_parity[Q_index][column_block][row_block]));

                size_t residue_row_offset(
                  ((column_block * (column_block + 1)) / 2 + row_block)
                  * block_size);

                El::DistMatrix<El::BigFloat> residue_sub_block(El::View(
                  *dual_residues_block, residue_row_offset, 0, block_size, 1));
                El::Axpy(El::BigFloat(-1.0), lower_diagonal,
                         residue_sub_block);
              }
        }
      // dualResidues -= B * y
      // TODO: Shouldn't this be Gemv since y is a vector?
      El::Gemm(El::Orientation::NORMAL, El::Orientation::NORMAL,
               El::BigFloat(-1), *free_var_matrix_block, *y_block,
               El::BigFloat(1), *dual_residues_block);
      // dualResidues += c
      El::Axpy(El::BigFloat(1), *primal_objective_c_block,
               *dual_residues_block);

      local_max = El::Max(local_max, El::MaxAbs(*dual_residues_block));

      ++primal_objective_c_block;
      ++y_block;
      ++free_var_matrix_block;
      ++dual_residues_block;
      ++Q_index;
    }
  dual_error
    = El::mpi::AllReduce(local_max, El::mpi::MAX, El::mpi::COMM_WORLD);
  dual_residues_timer.stop();
}

//////////////// compute tr(A Z) //////////////////////////

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

///////////////////////// Block_Vector ////////////////////////////////

El::BigFloat compute_Block_Vector_norm2(const Block_Info &block_info, const Block_Vector &vec)
{
	El::BigFloat norm2(0);

	auto vec_block(vec.blocks.begin());
	for (auto &block_index : block_info.block_indices)
	{
		int blockID = vec_block - vec.blocks.begin();

		int height = vec_block->Height();
		int width = vec_block->Width();

		if (width == 1)
		{
			for (int i = 0; i < height; i++)
				norm2 += (vec_block->Get(i, 0))*(vec_block->Get(i, 0));
		}
		else if (height == 1)
		{
			for (int i = 0; i < width; i++)
				norm2 += (vec_block->Get(0, i))*(vec_block->Get(0, i));
		}
		else
		{
			if (El::mpi::Rank() == 0)std::cout << "compute_Block_Vector_norm2 error : input is not a vector" << "\n";
			exit(0);
		}

		++vec_block;
	}
	return norm2;
}
