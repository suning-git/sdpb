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
		El::BigFloat primal_residue_p_2;
		primal_residue_p_2 = dot(primal_residue_p, primal_residue_p);
		if (El::mpi::Rank() == 0)std::cout << "Before Predictor p.p : " << primal_residue_p_2 << "\n";
	}

    // Compute the predictor solution for (dx, dX, dy, dY)
    beta_predictor
      = predictor_centering_parameter(parameters, is_primal_and_dual_feasible);
    compute_search_direction(block_info, sdp, *this, schur_complement_cholesky,
                             schur_off_diagonal, X_cholesky, beta_predictor,
                             mu, primal_residue_p, false, Q, dx, dX, dy, dY);
    predictor_timer.stop();

	{
		El::BigFloat primal_residue_p_2;
		primal_residue_p_2 = dot(primal_residue_p, primal_residue_p);
		if (El::mpi::Rank() == 0)std::cout << "After Predictor p.p : " << primal_residue_p_2 << "\n";

		El::BigFloat dx2, dy2;
		dx2 = dot(dx, dx);
		dy2 = dot(dy, dy);
		if (El::mpi::Rank() == 0)std::cout << "Predictor beta : " << beta_predictor << "\n";
		if (El::mpi::Rank() == 0)std::cout << "Predictor : |dx.dx|=" << dx2 <<
			", |dy.dy|=" << dy2 << "\n";
	}

    // Compute the corrector solution for (dx, dX, dy, dY)
    auto &corrector_timer(
      timers.add_and_start("run.step.computeSearchDirection(betaCorrector)"));
    beta_corrector = corrector_centering_parameter(
      parameters, X, dX, Y, dY, mu, is_primal_and_dual_feasible,
      total_psd_rows);

    compute_search_direction(block_info, sdp, *this, schur_complement_cholesky,
                             schur_off_diagonal, X_cholesky, beta_corrector,
                             mu, primal_residue_p, true, Q, dx, dX, dy, dY);

	{
		El::BigFloat primal_residue_p_2;
		primal_residue_p_2 = dot(primal_residue_p, primal_residue_p);
		if (El::mpi::Rank() == 0)std::cout << "After Corrector p.p : " << primal_residue_p_2 << "\n";

		El::BigFloat dx2, dy2;
		dx2 = dot(dx, dx);
		dy2 = dot(dy, dy);
		if (El::mpi::Rank() == 0)std::cout << "Corrector beta : " << beta_corrector << "\n";
		if (El::mpi::Rank() == 0)std::cout << "Corrector : |dx.dx|=" << dx2 <<
			", |dy.dy|=" << dy2 << "\n";
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

  step_timer.stop();
}
