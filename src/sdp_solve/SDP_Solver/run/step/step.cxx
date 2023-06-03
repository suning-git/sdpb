#include "../../../SDP_Solver.hxx"

// Tr(A B), where A and B are symmetric
El::BigFloat frobenius_product_symmetric(const Block_Diagonal_Matrix &A,
                                         const Block_Diagonal_Matrix &B);

void initialize_schur_complement_solver(
  const Block_Info &block_info, const SDP &sdp,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_X_inv,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_Y,
  const El::Grid &block_grid, Block_Diagonal_Matrix &schur_complement_cholesky,
  Block_Matrix &schur_off_diagonal, El::DistMatrix<El::BigFloat> &Q,
  Timers &timers);

void compute_search_direction(
  const Block_Info &block_info, const SDP &sdp, const SDP_Solver &solver,
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const Block_Diagonal_Matrix &X_cholesky, const El::BigFloat beta,
  const El::BigFloat &mu, const Block_Vector &primal_residue_p,
  const bool &is_corrector_phase, const El::DistMatrix<El::BigFloat> &Q,
  Block_Vector &dx, Block_Diagonal_Matrix &dX, Block_Vector &dy,
  Block_Diagonal_Matrix &dY);

El::BigFloat predictor_centering_parameter(const Solver_Parameters &parameters,
                                           const bool is_primal_dual_feasible);

El::BigFloat corrector_centering_parameter(
  const Solver_Parameters &parameters, const Block_Diagonal_Matrix &X,
  const Block_Diagonal_Matrix &dX, const Block_Diagonal_Matrix &Y,
  const Block_Diagonal_Matrix &dY, const El::BigFloat &mu,
  const bool is_primal_dual_feasible, const size_t &total_num_rows);

El::BigFloat
step_length(const Block_Diagonal_Matrix &MCholesky,
            const Block_Diagonal_Matrix &dM, const El::BigFloat &gamma,
            const std::string &timer_name, Timers &timers);

void compute_R_error(const std::size_t &total_psd_rows,
	const Block_Diagonal_Matrix &X, const Block_Diagonal_Matrix &Y, El::BigFloat & R_error);

El::BigFloat max_step_length(const Block_Diagonal_Matrix &MCholesky,
	const Block_Diagonal_Matrix &dM,
	const std::string &timer_name,
	Timers &timers);

El::BigFloat dot(const Block_Vector &A, const Block_Vector &B);

El::BigFloat min_eigenvalue(Block_Diagonal_Matrix &A);

El::BigFloat Block_Vector_sum_dot(const Block_Vector &dy);
El::BigFloat Block_Vector_sum_max(const Block_Vector &dy);

void compute_schur_complement(
	const Block_Info &block_info,
	const std::array<
	std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
	&A_X_inv,
	const std::array<
	std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
	&A_Y,
	Block_Diagonal_Matrix &schur_complement, Timers &timers);


void compute_Sx(
	const Block_Info &block_info, Block_Diagonal_Matrix &schur_complement, const Block_Vector &x, Block_Vector &result);
void compute_By(
	const Block_Info &block_info, const SDP &sdp, const Block_Vector &y, Block_Vector &result);
void compute_minus_d_minus_TrApZ(
	const Block_Info &block_info, const SDP &sdp, const SDP_Solver &solver,
	const Block_Diagonal_Matrix &X_cholesky, const El::BigFloat beta,
	const El::BigFloat &mu,
	const bool &is_corrector_phase,
	Block_Vector &dx, const Block_Diagonal_Matrix &dX, const Block_Diagonal_Matrix &dY, Block_Diagonal_Matrix &Z);

void Block_Vector_plus_equal(Block_Vector &B, const Block_Vector &A);
void Block_Vector_minus_equal(Block_Vector &B, const Block_Vector &A);
El::BigFloat Block_Vector_p_max(const Block_Info &block_info, const Block_Vector &Vp);

void compute_trAp_Z(const Block_Info &block_info, const SDP &sdp,
	const Block_Diagonal_Matrix &Z,
	Block_Vector &result);

void compute_trA_Y_V2(
	const Block_Info &block_info, const SDP &sdp,
	const std::array<
	std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
	&A_Y,
	Block_Vector &result);

El::BigFloat Block_Vector_p_max_V2(const Block_Info &block_info, const Block_Vector &vec);
El::BigFloat Block_Vector_p_max_V3(const Block_Info &block_info, Block_Vector &vec);

El::BigFloat Block_Vector_sum_dot_debug(const Block_Vector &dy);

void Assert_DistMatrix_Local(const El::DistMatrix<El::BigFloat> & mat, auto message);
El::BigFloat Block_Vector_Square(const Block_Info &block_info, const Block_Vector &vec);

void debug_print_Matrix_blocks_info(const Block_Info &block_info, const std::vector<El::DistMatrix<El::BigFloat>> & blocks, const std::string var);
void debug_print_Vector_blocks_info(const Block_Info &block_info, const std::vector<El::DistMatrix<El::BigFloat>> & blocks, const std::string var);
void debug_print_blocks_info(const Block_Info &block_info, const std::vector<El::DistMatrix<El::BigFloat>> & blocks, const std::string var);
void debug_Ap_V1_vs_V2(const Block_Info &block_info, const std::vector<El::DistMatrix<El::BigFloat>> & blocks, const std::string var);

void compute_trAp_Z_V2(const Block_Info &block_info, const SDP &sdp, const Block_Diagonal_Matrix &Z, Block_Vector &result);

void SDP_Solver::step(
  const Solver_Parameters &parameters, const std::size_t &total_psd_rows,
  const bool &is_primal_and_dual_feasible, const Block_Info &block_info,
  const SDP &sdp, const El::Grid &grid,
  const Block_Diagonal_Matrix &X_cholesky,
  const Block_Diagonal_Matrix &Y_cholesky,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_X_inv,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_Y,
  const Block_Vector &primal_residue_p, El::BigFloat &mu,
  El::BigFloat &beta_corrector, El::BigFloat &primal_step_length,
  El::BigFloat &dual_step_length, bool &terminate_now, Timers &timers)
{
  El::BigFloat R_error;


  auto &step_timer(timers.add_and_start("run.step"));
  El::BigFloat beta_predictor;

  // Search direction: These quantities have the same structure
  // as (x, X, y, Y). They are computed twice each iteration:
  // once in the predictor step, and once in the corrector step.
  Block_Vector dx(x), dy(y);
  Block_Diagonal_Matrix dX(X), dY(Y);
  {
    // SchurComplementCholesky = L', the Cholesky decomposition of the
    // Schur complement matrix S.
    Block_Diagonal_Matrix schur_complement_cholesky(
      block_info.schur_block_sizes(), block_info.block_indices,
      block_info.num_points.size(), grid);

    // SchurOffDiagonal = L'^{-1} FreeVarMatrix, needed in solving the
    // Schur complement equation.
    Block_Matrix schur_off_diagonal;

    // Q = B' L'^{-T} L'^{-1} B' - {{0, 0}, {0, 1}}, where B' =
    // (FreeVarMatrix U).  Q is needed in the factorization of the Schur
    // complement equation.  Q has dimension N'xN', where
    //
    //   N' = cols(B) + cols(U) = N + cols(U)
    //
    // where N is the dimension of the dual objective function.  Note
    // that N' could change with each iteration.
    El::DistMatrix<El::BigFloat> Q(sdp.dual_objective_b.Height(),
                                   sdp.dual_objective_b.Height());

    // Compute SchurComplement and prepare to solve the Schur
    // complement equation for dx, dy
    initialize_schur_complement_solver(block_info, sdp, A_X_inv, A_Y, grid,
                                       schur_complement_cholesky,
                                       schur_off_diagonal, Q, timers);

	///////// compute S_ij for debug purpose ///////////////////
	Block_Diagonal_Matrix schur_complement(
		block_info.schur_block_sizes(), block_info.block_indices,
		block_info.num_points.size(), grid);
	compute_schur_complement(block_info, A_X_inv, A_Y, schur_complement,
		timers);
	////////////////////////////////////////////////////////////

    // Compute the complementarity mu = Tr(X Y)/X.dim
    auto &frobenius_timer(
      timers.add_and_start("run.step.frobenius_product_symmetric"));
    mu = frobenius_product_symmetric(X, Y) / total_psd_rows;
    frobenius_timer.stop();
    if(mu > parameters.max_complementarity)
      {
        terminate_now = true;
        return;
      }

    auto &predictor_timer(
      timers.add_and_start("run.step.computeSearchDirection(betaPredictor)"));

	{
 		El::BigFloat primal_residue_p_2 = Block_Vector_sum_max(primal_residue_p);
		if (El::mpi::Rank() == 0)std::cout << "Block_Vector_sum_max(primal_residue_p) =" << primal_residue_p_2 << "\n";

		//El::mpi::Barrier(block_info.mpi_comm.value);
		//if (El::mpi::Rank() == 0)std::cout << "----------- p info ---------------" << "\n" << std::flush;
		El::mpi::Barrier(MPI_COMM_WORLD); usleep(100000);
		//debug_print_Vector_blocks_info(block_info, primal_residue_p.blocks, "p");
		
		//usleep(1000000);

		//if (El::mpi::Rank() == 0)std::cout << "----------- dy_dist info ---------------" << "\n" << std::flush;

		//El::BigFloat primal_residue_p_dot = Block_Vector_sum_dot_debug(primal_residue_p);
		//if (El::mpi::Rank() == 0)std::cout << "Block_Vector_sum_dot_debug(primal_residue_p) =" << primal_residue_p_dot << "\n";
		//usleep(1000000);
	}

    // Compute the predictor solution for (dx, dX, dy, dY)
    beta_predictor
      = predictor_centering_parameter(parameters, is_primal_and_dual_feasible);
    compute_search_direction(block_info, sdp, *this, schur_complement_cholesky,
                             schur_off_diagonal, X_cholesky, beta_predictor,
                             mu, primal_residue_p, false, Q, dx, dX, dy, dY);
    predictor_timer.stop();

	// try to check Sdx
	{
		Block_Vector Bdy(x);
		compute_By(block_info, sdp, dy, Bdy);

		Block_Vector minus_d_minus_TrApZBdy(x);
		Block_Diagonal_Matrix Z(X);
		compute_minus_d_minus_TrApZ(block_info, sdp, *this, X_cholesky, beta_predictor,
			mu, false, minus_d_minus_TrApZBdy, dX, dY, Z);

		Block_Vector Sdx(x);
		compute_Sx(block_info, schur_complement, dx, Sdx);

		Block_Vector result_test(Sdx);
		Block_Vector_minus_equal(result_test, Bdy);
		Block_Vector_minus_equal(result_test, minus_d_minus_TrApZBdy);

		El::BigFloat result_max = Block_Vector_p_max_V3(block_info, result_test);

		if (El::mpi::Rank() == 0)std::cout << "[Predictor] Sdx-Bdy+d+tr(AZ)=" << result_max << "\n";

		Block_Diagonal_Matrix dY_plus_Z(dY);
		dY_plus_Z += Z;
		Block_Vector trAp_dY_plus_Z(x);

		compute_trAp_Z(block_info, sdp, dY_plus_Z, trAp_dY_plus_Z);

		Block_Vector result_test2(trAp_dY_plus_Z);
		Block_Vector_plus_equal(result_test2, Sdx);

		El::BigFloat result2_max = Block_Vector_p_max_V3(block_info, result_test2);

		if (El::mpi::Rank() == 0)std::cout << "[Predictor] tr(A(dY+Z)+Sdx)=" << result2_max << "\n";
	}

	{   
   /*
    Block_Vector dy_plus_p(dy);
    for(size_t block = 0; block < dy.blocks.size(); ++block)
    {
      El::Axpy(1, primal_residue_p.blocks[block], dy.blocks[block]);
    }
    El::BigFloat test=dot(dy_plus_p, dy_plus_p);
    if (El::mpi::Rank() == 0)std::cout << "compute dy+p :" << test << "\n";
    */

		El::BigFloat dx2 = dot(dx, dx);
		El::BigFloat dy2 = dot(dy, dy);
		El::BigFloat dx2crt = 0; // Block_Vector_sum_dot(dx);
		El::BigFloat dy2crt = 0; // Block_Vector_sum_dot(dy);

		if (El::mpi::Rank() == 0)std::cout << "Predictor : |dx.dx|=" << dx2 <<
			", sum_dot(dx)=" << dx2crt <<
			", |dy.dy|=" << dy2 << ", sum_dot(dy)=" << dy2crt << "\n";
	}

    // Compute the corrector solution for (dx, dX, dy, dY)
    auto &corrector_timer(
      timers.add_and_start("run.step.computeSearchDirection(betaCorrector)"));
    beta_corrector = corrector_centering_parameter(
      parameters, X, dX, Y, dY, mu, is_primal_and_dual_feasible,
      total_psd_rows);

    //compute_search_direction(block_info, sdp, *this, schur_complement_cholesky,
    //                         schur_off_diagonal, X_cholesky, beta_corrector,
    //                         mu, primal_residue_p, true, Q, dx, dX, dy, dY);

	{
		Block_Vector minus_d_minus_TrApZBdy(x);
		Block_Diagonal_Matrix Z(X);
		compute_minus_d_minus_TrApZ(block_info, sdp, *this, X_cholesky, beta_corrector,
			mu, true, minus_d_minus_TrApZBdy, dX, dY, Z);

		compute_search_direction(block_info, sdp, *this, schur_complement_cholesky,
		                         schur_off_diagonal, X_cholesky, beta_corrector,
		                         mu, primal_residue_p, true, Q, dx, dX, dy, dY);

		Block_Vector Bdy(x);
		compute_By(block_info, sdp, dy, Bdy);

		Block_Vector Sdx(x);
		compute_Sx(block_info, schur_complement, dx, Sdx);

		Block_Vector result_test(Sdx);
		Block_Vector_minus_equal(result_test, Bdy);
		Block_Vector_minus_equal(result_test, minus_d_minus_TrApZBdy);

		El::BigFloat result_max = Block_Vector_p_max_V3(block_info, result_test);
		if (El::mpi::Rank() == 0)std::cout << "[Corrector] Sdx-Bdy+d+tr(AZ)=" << result_max << "\n";
		//result_max = Block_Vector_sum_max(result_test);
		//if (El::mpi::Rank() == 0)std::cout << "[Corrector]#Sdx-Bdy+d+tr(AZ)=" << result_max << "\n";

		Block_Diagonal_Matrix dY_plus_Z(dY);
		dY_plus_Z += Z;
		Block_Vector trAp_dY_plus_Z(x);

		compute_trAp_Z(block_info, sdp, dY_plus_Z, trAp_dY_plus_Z);

		Block_Vector result_test2(trAp_dY_plus_Z);
		Block_Vector_plus_equal(result_test2, Sdx);

		El::BigFloat result2_max = Block_Vector_p_max_V3(block_info, result_test2);
		if (El::mpi::Rank() == 0)std::cout << "[Corrector] tr(A(dY+Z)+Sdx)=" << result2_max << "\n";
		//result2_max = Block_Vector_sum_max(result_test2);
		//if (El::mpi::Rank() == 0)std::cout << "[Corrector]#tr(A(dY+Z)+Sdx)=" << result2_max << "\n";

		Block_Vector trAp_dY_plus_Z_V2(x);
		compute_trAp_Z_V2(block_info, sdp, dY_plus_Z, trAp_dY_plus_Z_V2);
		Block_Vector trAp_dY_plus_Z_plus_Sdx(trAp_dY_plus_Z_V2);
		Block_Vector_plus_equal(trAp_dY_plus_Z_plus_Sdx, Sdx);
		El::BigFloat trAp_dY_plus_Z_plus_Sdx_V2_max = Block_Vector_p_max_V3(block_info, trAp_dY_plus_Z_plus_Sdx);
		if (El::mpi::Rank() == 0)std::cout << "[Corrector] tr(A(dY+Z)+Sdx)_V2=" << trAp_dY_plus_Z_plus_Sdx_V2_max << "\n";


		Block_Vector trApY(x);
		compute_trAp_Z(block_info, sdp, Y, trApY);
		Block_Vector By(x);
		compute_By(block_info, sdp, y, By);
		Block_Vector result_old_d(trApY);
		Block_Vector_plus_equal(result_old_d, By);
		Block_Vector_minus_equal(result_old_d, sdp.primal_objective_c);

		El::BigFloat result_old_d_max = Block_Vector_p_max_V3(block_info, result_old_d);
		if (El::mpi::Rank() == 0)std::cout << "[Corrector] tr(Ap Y)+By-c=" << result_old_d_max << "\n";
		//result_old_d_max = Block_Vector_sum_max(result_old_d);
		//if (El::mpi::Rank() == 0)std::cout << "[Corrector]#tr(Ap Y)+By-c=" << result_old_d_max << "\n";


		Block_Vector trApY_V2(x);

		compute_trA_Y_V2(block_info, sdp, A_Y, trApY_V2);

		/*
		if (El::mpi::Rank() == 0)std::cout << "----------- B info ---------------" << "\n" << std::flush;
		El::mpi::Barrier(MPI_COMM_WORLD); usleep(100000);
		debug_print_blocks_info(block_info, sdp.bilinear_bases, "[B]");
		El::mpi::Barrier(MPI_COMM_WORLD); usleep(100000);

		if (El::mpi::Rank() == 0)std::cout << "----------- Z info ---------------" << "\n" << std::flush;
		El::mpi::Barrier(MPI_COMM_WORLD); usleep(100000);
		debug_print_blocks_info(block_info, Z.blocks, "[Z]");
		El::mpi::Barrier(MPI_COMM_WORLD); usleep(100000);
		*/

		/*
		for (auto &block_index : block_info.block_indices)
		{
			std::cout << "Rank=" << El::mpi::Rank() << ", globalID=" << block_index
				<< ", num_points=" << block_info.num_points[block_index] 
				<< ", dimensions=" << block_info.dimensions[block_index] << "\n";
		}
		*/

		/*
		if (El::mpi::Rank() == 0)std::cout << "----------- tr(A Y) V1 info ---------------" << "\n" << std::flush;
		El::mpi::Barrier(MPI_COMM_WORLD); usleep(100000);
		debug_print_blocks_info(block_info, trApY.blocks, "[trAY V1]");
		El::mpi::Barrier(MPI_COMM_WORLD); usleep(100000);

		if (El::mpi::Rank() == 0)std::cout << "----------- tr(A Y) V2 info ---------------" << "\n" << std::flush;
		El::mpi::Barrier(MPI_COMM_WORLD); usleep(100000);
		debug_print_blocks_info(block_info, trApY_V2.blocks, "[trAY V2]");
		El::mpi::Barrier(MPI_COMM_WORLD); usleep(100000);

		if (El::mpi::Rank() == 0)std::cout << std::flush << "----------- tr(A Y) V1 data ---------------" << "\n" << std::flush;
		El::mpi::Barrier(MPI_COMM_WORLD); usleep(100000);
		debug_Ap_V1_vs_V2(block_info, trApY.blocks, "[trAY V1]");
		El::mpi::Barrier(MPI_COMM_WORLD); usleep(1000000);

		if (El::mpi::Rank() == 0)std::cout << std::flush << "----------- tr(A Y) V2 data ---------------" << "\n" << std::flush;
		El::mpi::Barrier(MPI_COMM_WORLD); usleep(100000);
		debug_Ap_V1_vs_V2(block_info, trApY_V2.blocks, "[trAY V2]");
		El::mpi::Barrier(MPI_COMM_WORLD); usleep(1000000);
		*/

		Block_Vector result_old_d_V2(trApY_V2);
		Block_Vector_plus_equal(result_old_d_V2, By);
		Block_Vector_minus_equal(result_old_d_V2, sdp.primal_objective_c);

		El::BigFloat result_old_d_V2_max = Block_Vector_p_max_V3(block_info, result_old_d_V2);
		if (El::mpi::Rank() == 0)std::cout << "[Corrector] V2 : tr(Ap Y)+By-c=" << result_old_d_V2_max << "\n";
		//result_old_d_V2_max = Block_Vector_sum_max(result_old_d_V2);
		//if (El::mpi::Rank() == 0)std::cout << "[Corrector] V2 :#tr(Ap Y)+By-c=" << result_old_d_V2_max << "\n";

		Block_Vector trApY_V2_minus_V1(trApY_V2); 
		Block_Vector_minus_equal(trApY_V2_minus_V1, trApY);

		El::BigFloat trApY_V2_minus_V1_max = Block_Vector_p_max_V3(block_info, trApY_V2_minus_V1);
		if (El::mpi::Rank() == 0)std::cout << "[Corrector] : tr(Ap Y)_V2 - tr(Ap Y)_V1=" << trApY_V2_minus_V1_max << "\n";
		//trApY_V2_minus_V1_max = Block_Vector_sum_max(trApY_V2_minus_V1);
		//if (El::mpi::Rank() == 0)std::cout << "[Corrector] :#tr(Ap Y)_V2 - tr(Ap Y)_V1=" << trApY_V2_minus_V1_max << "\n";
		

		/*
		/////////////// print x,y,X,Y information ///////////////////
		//El::mpi::Barrier(block_info.mpi_comm.value);
		if (El::mpi::Rank() == 0)std::cout << "----------- x info ---------------" << "\n" << std::flush;
		El::mpi::Barrier(MPI_COMM_WORLD); usleep(100000);
		debug_print_Vector_blocks_info(block_info, x.blocks, "x");
		if (El::mpi::Rank() == 0)std::cout << "----------- y info ---------------" << "\n" << std::flush;
		El::mpi::Barrier(MPI_COMM_WORLD); usleep(100000);
		debug_print_Vector_blocks_info(block_info, y.blocks, "y");
		if (El::mpi::Rank() == 0)std::cout << "----------- X info ---------------" << "\n" << std::flush;
		El::mpi::Barrier(MPI_COMM_WORLD); usleep(100000);
		debug_print_Matrix_blocks_info(block_info, X.blocks, "X");
		if (El::mpi::Rank() == 0)std::cout << "----------- Y info ---------------" << "\n" << std::flush;
		El::mpi::Barrier(MPI_COMM_WORLD); usleep(100000);
		debug_print_Matrix_blocks_info(block_info, Y.blocks, "Y");
		El::mpi::Barrier(MPI_COMM_WORLD); usleep(100000);
		/////////////////////////////////////////////////////////////
		*/


		//if (El::mpi::Rank() == 0)std::cout << "----------- debug tr(Ap Y) V2 ---------------" << "\n";
		//if (El::mpi::Rank() == 0 || El::mpi::Rank() == 1)
		/*{
			std::cout << "[Rank" << El::mpi::Rank() << "] blocks.size=" << trApY_V2.blocks.size() << "\n";

			auto trApY_block(trApY_V2.blocks.begin());
			for (auto &block_index : block_info.block_indices)
			{
				int blockID = trApY_block - trApY_V2.blocks.begin();

				std::cout << "[Rank" << El::mpi::Rank() << "]" << "block " << blockID
					<< " height="<< trApY_block->Height() << " width=" << trApY_block->Width() << "\n";

				for (int i = 0; i < trApY_block->Height() && i < 4; i++)
				{
					std::cout << "[Rank" << El::mpi::Rank() << "]" << "block " << blockID
						<< "[" << i << "," << 0 << "]=" << trApY_block->Get(i, 0) << "\n";
				}

				++trApY_block;
			}
		}
		*/

		/*

		if (El::mpi::Rank() == 0)std::cout << "----------- debug tr(Ap Y) V1 ---------------" << "\n";

		if (El::mpi::Rank() == 0 || El::mpi::Rank() == 1)
		{
			std::cout << "[Rank" << El::mpi::Rank() << "] blocks.size=" << trApY.blocks.size() << "\n";

			auto trApY_block(trApY.blocks.begin());
			for (auto &block_index : block_info.block_indices)
			{
				int blockID = trApY_block - trApY.blocks.begin();

				std::cout << "[Rank" << El::mpi::Rank() << "]" << "block " << blockID
					<< " height=" << trApY_block->Height() << " width=" << trApY_block->Width() << "\n";

				for (int i = 0; i < trApY_block->Height() && i < 4; i++)
				{
					std::cout << "[Rank" << El::mpi::Rank() << "]" << "block " << blockID
						<< "[" << 0 << "," << i << "]=" << trApY_block->Get(i, 0) << "\n";
				}
				++trApY_block;
			}
		}

		if (El::mpi::Rank() == 0)std::cout << "---------------------------------------------" << "\n";
		*/

	}


	{
		El::BigFloat dx2 = dot(dx, dx);
		El::BigFloat dy2 = dot(dy, dy);
		El::BigFloat dx2crt = Block_Vector_sum_dot(dx);
		El::BigFloat dy2crt = Block_Vector_sum_dot(dy);

		if (El::mpi::Rank() == 0)std::cout << "Corrector : |dx.dx|=" << dx2 <<
			", sum_dot(dx)=" << dx2crt <<
			", |dy.dy|=" << dy2 << ", sum_dot(dy)=" << dy2crt << "\n";
	}

    corrector_timer.stop();
  }
  // Compute step-lengths that preserve positive definiteness of X, Y
  primal_step_length
    = step_length(X_cholesky, dX, parameters.step_length_reduction,
                  "run.step.stepLength(XCholesky)", timers);

  dual_step_length
    = step_length(Y_cholesky, dY, parameters.step_length_reduction,
                  "run.step.stepLength(YCholesky)", timers);


  El::BigFloat primal_step_maxlength, dual_step_maxlength;

  primal_step_maxlength = max_step_length(X_cholesky, dX, "run.step.stepLength(XCholesky)", timers);
  dual_step_maxlength = max_step_length(Y_cholesky, dY, "run.step.stepLength(YCholesky)", timers);

  // If our problem is both dual-feasible and primal-feasible,
  // ensure we're following the true Newton direction.
  if(is_primal_and_dual_feasible)
    {
      primal_step_length = El::Min(primal_step_length, dual_step_length);
      dual_step_length = primal_step_length;
    }

  // Update the primal point (x, X) += primalStepLength*(dx, dX)
  for(size_t block = 0; block < x.blocks.size(); ++block)
    {
      El::Axpy(primal_step_length, dx.blocks[block], x.blocks[block]);
    }
  dX *= primal_step_length;

  X += dX;

  // Update the dual point (y, Y) += dualStepLength*(dy, dY)
  for(size_t block = 0; block < dy.blocks.size(); ++block)
    {
      El::Axpy(dual_step_length, dy.blocks[block], y.blocks[block]);
    }
  dY *= dual_step_length;

  Y += dY;


  //////////////// debug Slater /////////////

  compute_R_error(total_psd_rows, X, Y, R_error);

  if (El::mpi::Rank() == 0)std::cout << "Rerr=" << R_error << "\n";

  if (El::mpi::Rank() == 0)std::cout << "max_p_step=" << primal_step_maxlength <<
	  ", max_d_step=" << dual_step_maxlength << "\n";
     
     
     Block_Diagonal_Matrix X_copy(X), Y_copy(Y);
     
     El::BigFloat Xmin=min_eigenvalue(X_copy);
     El::BigFloat Ymin=min_eigenvalue(Y_copy);
     
    if (El::mpi::Rank() == 0)std::cout << "minX=" << Xmin <<
	  ", Ymin=" << Ymin << "\n";

	////////////////// check the structure of x /////////////////////

	/*
	bool fakeblock = !x.blocks.empty() && x.blocks.at(0).Grid().Rank() != 0;
	bool fakeblock2 = !X.blocks.empty() && X.blocks.at(0).Grid().Rank() != 0;

	for (size_t block = 0; block != x.blocks.size(); ++block)
	{
		size_t block_index(block_info.block_indices.at(block));
		for (size_t psd_block(0); psd_block < 2; ++psd_block)
		{
			int blockID_global = 2 * block_index + psd_block;
			int blockID = 2 * block + psd_block;

			std::cout << "mpi::Rank=" << El::mpi::Rank() <<
				", fakeQ=" << fakeblock <<
				", fake2Q=" << fakeblock2 <<
				", globalID=" << blockID_global <<
				", localID=" << blockID <<
				", Grid.Rank=" << x.blocks.at(block).Grid().Rank() << 
				"\n";
		}
	}
	*/

	/*
	std::cout << "mpi::Rank=" << El::mpi::Rank() <<
		"len(block_indices)=" << block_info.block_indices.size() <<
		", len(x.blocks)=" << x.blocks.size() <<
		", len(y.blocks)=" << y.blocks.size() <<
		", len(X.blocks)=" << X.blocks.size() <<
		", len(Y.blocks)=" << Y.blocks.size() <<
		"\n";
		*/

  step_timer.stop();
}

