#include "Boost_Float.hxx"
#include "Positive_Matrix_With_Prefactor.hxx"
#include "../Timers.hxx"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <csignal>

namespace po = boost::program_options;

void read_input(const boost::filesystem::path &input_file,
	std::vector<El::BigFloat> &objectives,
	std::vector<El::BigFloat> &normalization,
	std::vector<Positive_Matrix_With_Prefactor> &matrices,
	std::vector<int> &matrices_valid_indices);

void write_output(const boost::filesystem::path &output_dir,
	const std::vector<El::BigFloat> &objectives,
	const std::vector<El::BigFloat> &normalization,
	const std::vector<Positive_Matrix_With_Prefactor> &matrices,
	std::vector<int> &matrices_valid_indices,
	Timers &timers);

void test_load();

int parse_parameter_file(boost::filesystem::path & param_file);

void sb_parse_input(std::vector<El::BigFloat> &objectives,
	std::vector<El::BigFloat> &normalization,
	std::vector<Positive_Matrix_With_Prefactor> &matrices,
	std::vector<int> &matrices_valid_indices
);

int main(int argc, char **argv)
{
  El::Environment env(argc, argv);

  try
    {
      int precision;
      boost::filesystem::path input_file, output_dir;
      bool debug(false);

	  boost::filesystem::path param_file;

      po::options_description options("Basic options");
      options.add_options()("help,h", "Show this helpful message.");
      options.add_options()(
        "input,i", po::value<boost::filesystem::path>(&input_file),
        "Mathematica, JSON, or NSV file with SDP definition");
      options.add_options()(
        "output,o",
        po::value<boost::filesystem::path>(&output_dir)->required(),
        "Directory to place output");
      options.add_options()(
        "precision", po::value<int>(&precision)->required(),
        "The precision, in the number of bits, for numbers in the "
        "computation. ");
      options.add_options()("debug",
                            po::value<bool>(&debug)->default_value(false),
                            "Write out debugging output.");


	  options.add_options()(
		  "parameter", po::value<boost::filesystem::path>(&param_file),
		  "*.m file that defined varialble used for sdp2input_mod");

      po::positional_options_description positional;
      positional.add("precision", 1);
      positional.add("input", 1);
      positional.add("output", 1);

      po::variables_map variables_map;
      po::store(po::parse_command_line(argc, argv, options), variables_map);

      if(variables_map.count("help") != 0)
        {
          std::cout << options << '\n';
          return 0;
        }

      po::notify(variables_map);


	  if (variables_map.count("parameter") == 0)
	  {
		  if (variables_map.count("input") == 0)
			  throw std::runtime_error("Either --parameter or --input is required.");

		  if (!boost::filesystem::exists(input_file))
		  {
			  throw std::runtime_error("Input file '" + input_file.string()
				  + "' does not exist");
		  }
		  if (boost::filesystem::is_directory(input_file))
		  {
			  throw std::runtime_error("Input file '" + input_file.string()
				  + "' is a directory, not a file");
		  }
	  }

      if(boost::filesystem::exists(output_dir)
         && !boost::filesystem::is_directory(output_dir))
        {
          throw std::runtime_error("Output directory '" + output_dir.string()
                                   + "' exists and is not a directory");
        }

      El::gmp::SetPrecision(precision);
      // El::gmp wants base-2 bits, but boost::multiprecision want
      // base-10 digits.
      Boost_Float::default_precision(precision * log(2) / log(10));

	  std::cout << "test 1 : Boost_Float default precision=" << mpfr_get_default_prec() << "\n";
	  std::cout << "test 2 : precision * log(2) / log(10)=" << precision * log(2) / log(10) << "\n";

	  if (variables_map.count("parameter") != 0)
		  parse_parameter_file(param_file); // I should do this after set the default precision

      std::vector<El::BigFloat> objectives, normalization;
      std::vector<Positive_Matrix_With_Prefactor> matrices;
	  std::vector<int> matrices_valid_indices;
      Timers timers(debug);
      auto &read_input_timer(timers.add_and_start("read_input"));

	  if (variables_map.count("parameter") == 0)
	  {
		  read_input(input_file, objectives, normalization, matrices, matrices_valid_indices);
	  }
	  else
	  {
		  auto begin_clock = std::chrono::system_clock::now();
		  std::time_t begin_time = std::chrono::system_clock::to_time_t(begin_clock);

		  std::cout << "[Rank " << El::mpi::Rank() << "] begin sb_parse_input at " << std::ctime(&begin_time);

		  sb_parse_input(objectives, normalization, matrices, matrices_valid_indices);

		  auto end_clock = std::chrono::system_clock::now();
		  std::time_t end_time = std::chrono::system_clock::to_time_t(end_clock);

		  std::chrono::duration<double> elapsed_seconds = end_clock - begin_clock;

		  std::cout << "[Rank " << El::mpi::Rank() << "] finished sb_parse_input at " << std::ctime(&end_time) << "elapsed time: " << elapsed_seconds.count() << "s\n";
	  }
      read_input_timer.stop();

	  {
		  auto begin_clock = std::chrono::system_clock::now();
		  std::time_t begin_time = std::chrono::system_clock::to_time_t(begin_clock);

		  std::cout << "[Rank " << El::mpi::Rank() << "] begin write_output at " << std::ctime(&begin_time);

		  auto &write_output_timer(timers.add_and_start("write_output"));
		  write_output(output_dir, objectives, normalization, matrices, matrices_valid_indices, timers);
		  write_output_timer.stop();

		  auto end_clock = std::chrono::system_clock::now();
		  std::time_t end_time = std::chrono::system_clock::to_time_t(end_clock);

		  std::chrono::duration<double> elapsed_seconds = end_clock - begin_clock;

		  std::cout << "[Rank " << El::mpi::Rank() << "] finished write_output at " << std::ctime(&end_time) << "elapsed time: " << elapsed_seconds.count() << "s\n";
	  }


      if(debug)
        {
          timers.write_profile(output_dir.string() + "/profiling."
                               + std::to_string(El::mpi::Rank()));
        }
    }
  catch(std::exception &e)
    {
      std::cerr << "Error: " << e.what() << "\n" << std::flush;
      El::mpi::Abort(El::mpi::COMM_WORLD, 1);
    }
  catch(...)
    {
      std::cerr << "Unknown Error\n" << std::flush;
      El::mpi::Abort(El::mpi::COMM_WORLD, 1);
    }
}
