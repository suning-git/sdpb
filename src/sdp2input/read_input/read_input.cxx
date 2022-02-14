#include "../Positive_Matrix_With_Prefactor.hxx"
#include "../../sdp_convert.hxx"

#include <boost/filesystem.hpp>

void read_json(const boost::filesystem::path &input_path,
               std::vector<El::BigFloat> &objectives,
               std::vector<El::BigFloat> &normalization,
               std::vector<Positive_Matrix_With_Prefactor> &matrices);

void read_mathematica(const boost::filesystem::path &input_path,
                      std::vector<El::BigFloat> &objectives,
                      std::vector<El::BigFloat> &normalization,
                      std::vector<Positive_Matrix_With_Prefactor> &matrices,
					  std::vector<int> &matrices_valid_indices);

void read_input(const boost::filesystem::path &input_file,
                std::vector<El::BigFloat> &objectives,
                std::vector<El::BigFloat> &normalization,
                std::vector<Positive_Matrix_With_Prefactor> &matrices,
				std::vector<int> &matrices_valid_indices)
{
	read_mathematica(input_file, objectives, normalization, matrices, matrices_valid_indices);
	return;
}
