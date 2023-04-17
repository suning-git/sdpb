import os, subprocess

def options(opt):
    opt.load(['compiler_cxx','gnu_dirs','cxx17','boost','gmpxx','mpfr',
              'elemental','libxml2', 'rapidjson'])

def configure(conf):
    if not 'CXX' in os.environ or os.environ['CXX']=='g++' or os.environ['CXX']=='icpc':
        conf.environ['CXX']='mpicxx'

    conf.load(['compiler_cxx','gnu_dirs','cxx17','boost','gmpxx','mpfr',
               'elemental','libxml2', 'rapidjson'])

    conf.env.git_version=subprocess.check_output('git describe --dirty', universal_newlines=True, shell=True).rstrip()
    
def build(bld):
    #default_flags=['-Wall', '-Wextra', '-O3', '-D SDPB_VERSION_STRING="' + bld.env.git_version + '"']
    default_flags=['-Wall', '-Wextra', '-g', '-D SDPB_VERSION_STRING="' + bld.env.git_version + '"']
    use_packages=['cxx17','boost','gmpxx','mpfr','elemental','libxml2', 'rapidjson']
    
    library_sources=['src/sdp_convert/Dual_Constraint_Group/Dual_Constraint_Group/Dual_Constraint_Group.cxx',
                     'src/sdp_convert/Dual_Constraint_Group/Dual_Constraint_Group/sample_bilinear_basis.cxx',
                     'src/sdp_convert/write_objectives.cxx',
                     'src/sdp_convert/write_bilinear_bases.cxx',
                     'src/sdp_convert/write_blocks.cxx',
                     'src/sdp_convert/write_primal_objective_c.cxx',
                     'src/sdp_convert/write_free_var_matrix.cxx',
                     'src/sdp_convert/write_sdpb_input_files.cxx',
                     'src/sdp_convert/read_file_list.cxx']

    bld.stlib(source=library_sources,
              target='sdp_convert',
              cxxflags=default_flags,
              use=use_packages)

    bld.program(source=['src/sdp2input/main.cxx',
                        'src/sdp2input/parse_MMA_expr.cxx',
                        'src/sdp2input/simpleboot.cxx',
                        'src/sdp2input/read_input/read_input.cxx',
                        'src/sdp2input/read_input/read_json/read_json.cxx',
                        'src/sdp2input/read_input/read_json/Positive_Matrix_With_Prefactor_State/json_key.cxx',
                        'src/sdp2input/read_input/read_json/Positive_Matrix_With_Prefactor_State/json_string.cxx',
                        'src/sdp2input/read_input/read_json/Positive_Matrix_With_Prefactor_State/json_start_array.cxx',
                        'src/sdp2input/read_input/read_json/Positive_Matrix_With_Prefactor_State/json_end_array.cxx',
                        'src/sdp2input/read_input/read_json/Positive_Matrix_With_Prefactor_State/json_start_object.cxx',
                        'src/sdp2input/read_input/read_json/Positive_Matrix_With_Prefactor_State/json_end_object.cxx',
                        'src/sdp2input/read_input/read_json/Damped_Rational_State/json_key.cxx',
                        'src/sdp2input/read_input/read_json/Damped_Rational_State/json_string.cxx',
                        'src/sdp2input/read_input/read_json/Damped_Rational_State/json_start_array.cxx',
                        'src/sdp2input/read_input/read_json/Damped_Rational_State/json_end_array.cxx',
                        'src/sdp2input/read_input/read_json/Damped_Rational_State/json_start_object.cxx',
                        'src/sdp2input/read_input/read_json/Damped_Rational_State/json_end_object.cxx',
                        'src/sdp2input/read_input/read_json/JSON_Parser/Key.cxx',
                        'src/sdp2input/read_input/read_json/JSON_Parser/String.cxx',
                        'src/sdp2input/read_input/read_json/JSON_Parser/StartArray.cxx',
                        'src/sdp2input/read_input/read_json/JSON_Parser/EndArray.cxx',
                        'src/sdp2input/read_input/read_json/JSON_Parser/StartObject.cxx',
                        'src/sdp2input/read_input/read_json/JSON_Parser/EndObject.cxx',
                        'src/sdp2input/read_input/read_mathematica/read_mathematica.cxx',
                        'src/sdp2input/read_input/read_mathematica/parse_SDP/parse_SDP.cxx',
                        'src/sdp2input/read_input/read_mathematica/parse_SDP/parse_matrices.cxx',
                        'src/sdp2input/read_input/read_mathematica/parse_SDP/parse_number.cxx',
                        'src/sdp2input/read_input/read_mathematica/parse_SDP/parse_polynomial.cxx',
                        'src/sdp2input/read_input/read_mathematica/parse_SDP/parse_matrix/parse_matrix.cxx',
                        'src/sdp2input/read_input/read_mathematica/parse_SDP/parse_matrix/parse_damped_rational.cxx',
                        'src/sdp2input/write_output/write_output.cxx',
                        'src/sdp2input/write_output/sample_points.cxx',
                        'src/sdp2input/write_output/bilinear_basis/bilinear_basis.cxx',
                        'src/sdp2input/write_output/bilinear_basis/precompute/precompute.cxx',
                        'src/sdp2input/write_output/bilinear_basis/precompute/integral.cxx',
                        'src/sdp2input/write_output/bilinear_basis/bilinear_form/bilinear_form.cxx',
                        'src/sdp2input/write_output/bilinear_basis/bilinear_form/rest.cxx',
                        'src/sdp2input/write_output/bilinear_basis/bilinear_form/dExp.cxx',
                        'src/sdp2input/write_output/bilinear_basis/bilinear_form/derivative.cxx',
                        'src/sdp2input/write_output/bilinear_basis/bilinear_form/operator_plus_set_Derivative_Term.cxx'],
                target='sdp2input_mod_2.4.0',
                cxxflags=default_flags,
                use=use_packages + ['sdp_convert']
                )
                
                        
    
