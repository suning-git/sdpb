import os, subprocess

def options(opt):
    opt.load(['compiler_cxx','gnu_dirs','cxx14','boost','gmpxx','mpfr',
              'elemental','libxml2', 'rapidjson'])

def configure(conf):
    if not 'CXX' in os.environ or os.environ['CXX']=='g++' or os.environ['CXX']=='icpc':
        conf.environ['CXX']='mpicxx'

    conf.load(['compiler_cxx','gnu_dirs','cxx14','boost','gmpxx','mpfr',
               'elemental','libxml2', 'rapidjson'])

    conf.env.git_version=subprocess.check_output('git describe --dirty', universal_newlines=True, shell=True).rstrip()
    
def build(bld):
    default_flags=['-Wall', '-Wextra', '-O3', '-D SDPB_VERSION_STRING="' + bld.env.git_version + '"']
    # default_flags=['-Wall', '-Wextra', '-g', '-D SDPB_VERSION_STRING="' + bld.env.git_version + '"']
    use_packages=['cxx14','boost','gmpxx','mpfr','elemental','libxml2', 'rapidjson']
    
    # Main executable
    bld.program(source=['src/sdpb/main.cxx',
                        'src/sdpb/Write_Solution.cxx',
                        'src/sdpb/SDP_Solver_Parameters/SDP_Solver_Parameters.cxx',
                        'src/sdpb/SDP_Solver_Parameters/ostream.cxx',
                        'src/sdpb/SDP_Solver_Parameters/to_property_tree.cxx',
                        'src/sdpb/solve/solve.cxx',
                        'src/compute_block_grid_mapping.cxx',
                        'src/sdpb/Block_Info/Block_Info.cxx',
                        'src/sdpb/Block_Info/read_block_info.cxx',
                        'src/sdpb/Block_Info/read_block_costs.cxx',
                        'src/sdpb/Block_Info/allocate_blocks.cxx',
                        'src/sdpb/write_timing.cxx',
                        'src/sdpb/solve/SDP/SDP/SDP.cxx',
                        'src/sdpb/solve/SDP/SDP/read_objectives.cxx',
                        'src/sdpb/solve/SDP/SDP/read_bilinear_bases.cxx',
                        'src/sdpb/solve/SDP/SDP/read_primal_objective_c.cxx',
                        'src/sdpb/solve/SDP/SDP/read_free_var_matrix.cxx',
                        'src/sdpb/solve/SDP_Solver/save_solution.cxx',
                        'src/sdpb/solve/SDP_Solver/save_checkpoint.cxx',
                        'src/sdpb/solve/SDP_Solver/load_checkpoint/load_checkpoint.cxx',
                        'src/sdpb/solve/SDP_Solver/load_checkpoint/load_binary_checkpoint.cxx',
                        'src/sdpb/solve/SDP_Solver/load_checkpoint/load_text_checkpoint.cxx',
                        'src/sdpb/solve/SDP_Solver/SDP_Solver.cxx',
                        'src/sdpb/solve/SDP_Solver/run/run.cxx',
                        'src/sdpb/solve/SDP_Solver/run/cholesky_decomposition.cxx',
                        'src/sdpb/solve/SDP_Solver/run/constraint_matrix_weighted_sum.cxx',
                        'src/sdpb/solve/SDP_Solver/run/compute_dual_residues_and_error.cxx',
                        'src/sdpb/solve/SDP_Solver/run/compute_primal_residues_and_error_P_Ax_X.cxx',
                        'src/sdpb/solve/SDP_Solver/run/compute_primal_residues_and_error_p_b_Bx.cxx',
                        'src/sdpb/solve/SDP_Solver/run/compute_objectives/compute_objectives.cxx',
                        'src/sdpb/solve/SDP_Solver/run/compute_objectives/dot.cxx',
                        'src/sdpb/solve/SDP_Solver/run/compute_bilinear_pairings/compute_bilinear_pairings.cxx',
                        'src/sdpb/solve/SDP_Solver/run/compute_bilinear_pairings/compute_bilinear_pairings_X_inv.cxx',
                        'src/sdpb/solve/SDP_Solver/run/compute_bilinear_pairings/compute_bilinear_pairings_Y.cxx',
                        'src/sdpb/solve/SDP_Solver/run/compute_feasible_and_termination.cxx',
                        'src/sdpb/solve/SDP_Solver/run/print_header.cxx',
                        'src/sdpb/solve/SDP_Solver/run/print_iteration.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/step.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/initialize_schur_complement_solver/initialize_schur_complement_solver.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/initialize_schur_complement_solver/compute_schur_complement.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/initialize_schur_complement_solver/initialize_Q_group.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/initialize_schur_complement_solver/synchronize_Q.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/compute_search_direction/compute_search_direction.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/compute_search_direction/cholesky_solve.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/compute_search_direction/compute_schur_RHS.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/compute_search_direction/scale_multiply_add.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/compute_search_direction/solve_schur_complement_equation.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/predictor_centering_parameter.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/corrector_centering_parameter/corrector_centering_parameter.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/corrector_centering_parameter/frobenius_product_of_sums.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/frobenius_product_symmetric.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/step_length/step_length.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/step_length/min_eigenvalue.cxx',
                        'src/sdpb/solve/SDP_Solver/run/step/step_length/lower_triangular_inverse_congruence.cxx',
                        'src/sdpb/solve/SDP_Solver_Terminate_Reason/ostream.cxx',
                        'src/sdpb/solve/lower_triangular_transpose_solve.cxx',
                        'src/sdpb/solve/Block_Diagonal_Matrix/ostream.cxx'],
                target='sdpd',
                cxxflags=default_flags,
                use=use_packages
                )
    
