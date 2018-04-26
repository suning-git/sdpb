#include "../Matrix.hxx"

// y := alpha A^T x + beta y
void vector_scale_matrix_multiply_transpose_add(const Real alpha,
                                                const Matrix &A,
                                                const Vector &x,
                                                const Real beta, Vector &y)
{
  assert(A.cols <= y.size());
  assert(A.rows <= x.size());

  for(size_t n = 0; n < A.cols; n++)
    {
      Real tmp = 0;
      for(size_t p = 0; p < A.rows; p++)
        {
          tmp += A.elt(p, n) * x[p];
        }
      y[n] = alpha * tmp + beta * y[n];
    }
}