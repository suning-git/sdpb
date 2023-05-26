#include "../../../../sdp_solve/SDP_Solver.hxx"


#include "../../../Block_Diagonal_Matrix.hxx"

// C := alpha*A*B + beta*C
void scale_multiply_add(const El::BigFloat &alpha,
	const Block_Diagonal_Matrix &A,
	const Block_Diagonal_Matrix &B,
	const El::BigFloat &beta, Block_Diagonal_Matrix &C)
{
	for (size_t block = 0; block < A.blocks.size(); ++block)
	{
		El::Gemm(El::OrientationNS::NORMAL, El::OrientationNS::NORMAL, alpha,
			A.blocks[block], B.blocks[block], beta, C.blocks[block]);
	}
}



// Tr(A B), where A and B are symmetric
El::BigFloat frobenius_product_symmetric(const Block_Diagonal_Matrix &A,
	const Block_Diagonal_Matrix &B);

// R_error= tr(XY)/X.dim * I - XY

void compute_R_error(const std::size_t &total_psd_rows, const Block_Diagonal_Matrix &X, const Block_Diagonal_Matrix &Y, 
	Block_Diagonal_Matrix &R, El::BigFloat & R_error, El::BigFloat & mu)
{
	mu = frobenius_product_symmetric(X, Y) / total_psd_rows;

	R = X;
	scale_multiply_add(El::BigFloat(-1), X, Y, El::BigFloat(0), R);
	R.add_diagonal(mu);

	R_error = R.max_abs_mpi();
}

void compute_R_error(const std::size_t &total_psd_rows, const Block_Diagonal_Matrix &X, const Block_Diagonal_Matrix &Y
	, El::BigFloat & R_error, El::BigFloat & mu)
{
	Block_Diagonal_Matrix R(X);
	return compute_R_error(total_psd_rows, X, Y, R, R_error, mu);
}

void compute_R_error(const std::size_t &total_psd_rows,
	const Block_Diagonal_Matrix &X, const Block_Diagonal_Matrix &Y, El::BigFloat & R_error)
{
	El::BigFloat mu;
	return compute_R_error(total_psd_rows, X, Y, R_error, mu);
}
