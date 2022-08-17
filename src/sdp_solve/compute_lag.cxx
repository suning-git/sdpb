#include "SDP_Solver.hxx"
#include <El.hpp>
#include <math.h>
#include <iostream>
//#include <gmpxx.h>
//#include <boost/multiprecision/gmp.hpp>
//#include <mpfr.h>
//#include <mpf2mpfr.h>
//using namespace boost::multiprecision;

//Does the same thing as Approx_Objective computes the linear difference 

// Tr(A B), where A and B are symmetric
El::BigFloat frobenius_product_symmetric(const Block_Diagonal_Matrix &A,
                                         const Block_Diagonal_Matrix &B);


El::Base<El::BigFloat> det_log_cholesky(const Block_Diagonal_Matrix &L){
   El::Base<El::BigFloat> local_det_log(0);
   for(size_t b = 0; b < L.blocks.size(); b++)
     {
       El::SafeProduct<El::Base<El::BigFloat>> safeDet = El::hpd_det::AfterCholesky(El::UpperOrLowerNS::LOWER,L.blocks[b]);
       local_det_log += safeDet.kappa*safeDet.n;
       //Notice that El::hpd_det::AfterCholesky multiply the det by 2 so the result is det(X = LL^T) instead of det(L);
     }
   //std::cout << "before all reduce: " << local_det_log << std::endl;
   El::Base<El::BigFloat> final_result = El::mpi::AllReduce(local_det_log, El::mpi::SUM, El::mpi::COMM_WORLD);
   //std::cout << "after all reduce: " << final_result   << std::endl;

   return final_result;
  
}


// solver 
// dual_residues[p] = c[p] - A[p,a,b] Y[a,b] - B[p,a] y[a]
// Lagrangian = dual_residues.x + b.y + Tr(X,Y) - mu log det (X) 

El::BigFloat compute_lag(const El::BigFloat mu, const Block_Diagonal_Matrix &X_cholesky, const SDP_Solver &solver){
  El::BigFloat lag(0);
  // dual_objective = f + b . y
  lag = solver.dual_objective;
  
  // Tr(XY)
  lag += frobenius_product_symmetric(solver.X, solver.Y);

  
  // mu log det X = mu log det LL^T = 2 mu log det L; 
  // The El::HPDDeterminant routine takes log of the diagonal entries and then exponentiate the sum
  // depending on if we want to use the whole routine or just use the sum of the logs  
  //Notice that El::hpd_det::AfterCholesky already accounts for the factor of 2 
  lag -= mu * (det_log_cholesky(X_cholesky)); 

  El::BigFloat local_residues(0);
  for(size_t block(0); block != solver.x.blocks.size(); ++block)
    {
      // dual_residues[p] = c[p] - A[p,a,b] Y[a,b] - B[p,a] y[a]
      // Lagrangian += dual_residues.x 
      local_residues += Dotu(solver.dual_residues.blocks.at(block), solver.x.blocks.at(block));
   }

   lag 
    += El::mpi::AllReduce(local_residues, El::mpi::SUM, El::mpi::COMM_WORLD);


  return lag;
}
