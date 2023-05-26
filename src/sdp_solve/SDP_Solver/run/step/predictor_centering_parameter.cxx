#include "../../../SDP_Solver.hxx"

// Centering parameter \beta_p for the predictor step
El::BigFloat predictor_centering_parameter(const Solver_Parameters &parameters,
                                           const bool is_primal_dual_feasible)
{
  if (parameters.infeasible_centering_parameter == El::BigFloat(1) && parameters.feasible_centering_parameter == El::BigFloat(1))
	  return El::BigFloat(1);
  return is_primal_dual_feasible ? El::BigFloat(0)
                                 : parameters.infeasible_centering_parameter;
}
