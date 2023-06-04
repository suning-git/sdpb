#include "../../../Block_Diagonal_Matrix.hxx"
#include "../../../Block_Info.hxx"
#include "../../../../Timers.hxx"

void compute_A_X_inv(
  const Block_Info &block_info, const Block_Diagonal_Matrix &X_cholesky,
  const std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  std::array<std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>,
             2> &A_X_inv);

void compute_A_Y(
  const Block_Info &block_info, const Block_Diagonal_Matrix &Y,
  const std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  std::array<std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>,
             2> &A_Y);


// X := ACholesky^{-T} ACholesky^{-1} X = A^{-1} X
void cholesky_solve(const Block_Diagonal_Matrix &ACholesky,
	Block_Diagonal_Matrix &X);

void compute_bilinear_pairings(
  const Block_Info &block_info, const Block_Diagonal_Matrix &X_cholesky,
  const Block_Diagonal_Matrix &Y,
  const std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  std::array<std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>,
             2> &A_X_inv,
  std::array<std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>,
             2> &A_Y,
  Timers &timers)
{
  auto &congruence_timer(timers.add_and_start("run.bilinear_pairings"));
  
 
  compute_A_X_inv(block_info, X_cholesky, bases_blocks, A_X_inv);

  /*
  Block_Diagonal_Matrix InvX(Y);
  InvX *= 0;
  InvX.add_diagonal(1);
  cholesky_solve(X_cholesky, InvX);
  //InvX.symmetrize();
  compute_A_Y(block_info, InvX, bases_blocks, A_X_inv);
  */

  

  compute_A_Y(block_info, Y, bases_blocks, A_Y);
  congruence_timer.stop();
}
