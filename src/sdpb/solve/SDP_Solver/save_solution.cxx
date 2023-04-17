#include "../SDP_Solver.hxx"
#include "../../../set_stream_precision.hxx"

#include <boost/filesystem/fstream.hpp>

#include <iomanip>

namespace
{
  void write_psd_block(const boost::filesystem::path &outfile,
                       const El::DistMatrix<El::BigFloat> &block)
  {
    boost::filesystem::ofstream stream;
    if(block.DistRank() == block.Root())
      {
        stream.open(outfile);
      }
    El::Print(block,
              std::to_string(block.Height()) + " "
                + std::to_string(block.Width()),
              "\n", stream);
    if(block.DistRank() == block.Root())
      {
        stream << "\n";
        if(!stream.good())
          {
            throw std::runtime_error("Error when writing to: "
                                     + outfile.string());
          }
      }
  }
}



void SDP_Solver::save_solution(
	const SDP_Solver_Terminate_Reason terminate_reason,
	const std::pair<std::string, Timer> &timer_pair,
	const boost::filesystem::path &out_directory,
	const Write_Solution &write_solution, int index,
	const std::vector<size_t> &block_indices, const Verbosity &verbosity) const
{
	// Internally, El::Print() sync's everything to the root core and
	// outputs it from there.  So do not actually open the file on
	// anything but the root node.


	Block_Vector & dx = dsdp_sol_list[index]->dx;
	Block_Diagonal_Matrix & dX = dsdp_sol_list[index]->dX;
	Block_Vector & dy = dsdp_sol_list[index]->dy;
	Block_Diagonal_Matrix & dY = dsdp_sol_list[index]->dY;


	// y is duplicated among cores, so only need to print out copy on
	// the root node.
	if (write_solution.vector_y && !dy.blocks.empty())
	{
		const boost::filesystem::path y_path(out_directory / "y.txt");
		boost::filesystem::ofstream y_stream;
		if (El::mpi::Rank() == 0)
		{
			y_stream.open(y_path);
		}
		El::Print(dy.blocks.at(0),
			std::to_string(dy.blocks.at(0).Height()) + " "
			+ std::to_string(dy.blocks.at(0).Width()),
			"\n", y_stream);
		if (El::mpi::Rank() == 0)
		{
			y_stream << "\n";
			if (!y_stream.good())
			{
				throw std::runtime_error("Error when writing to: "
					+ y_path.string());
			}
		}
	}

	for (size_t block = 0; block != dx.blocks.size(); ++block)
	{
		size_t block_index(block_indices.at(block));
		if (write_solution.vector_x)
		{
			const boost::filesystem::path x_path(
				out_directory / ("x_" + std::to_string(block_index) + ".txt"));
			boost::filesystem::ofstream x_stream;
			if (dx.blocks.at(block).DistRank() == dx.blocks.at(block).Root())
			{
				x_stream.open(x_path);
			}
			El::Print(dx.blocks.at(block),
				std::to_string(dx.blocks.at(block).Height()) + " "
				+ std::to_string(dx.blocks.at(block).Width()),
				"\n", x_stream);
			if (dx.blocks.at(block).DistRank() == dx.blocks.at(block).Root())
			{
				x_stream << "\n";
				if (!x_stream.good())
				{
					throw std::runtime_error("Error when writing to: "
						+ x_path.string());
				}
			}
		}
		for (size_t psd_block(0); psd_block < 2; ++psd_block)
		{
			std::string suffix(std::to_string(2 * block_index + psd_block)
				+ ".txt");

			if (write_solution.matrix_X
				&& dX.blocks.at(2 * block + psd_block).Height() != 0)
			{
				write_psd_block(out_directory / ("X_matrix_" + suffix),
					dX.blocks.at(2 * block + psd_block));
			}
			if (write_solution.matrix_Y
				&& dY.blocks.at(2 * block + psd_block).Height() != 0)
			{
				write_psd_block(out_directory / ("Y_matrix_" + suffix),
					dY.blocks.at(2 * block + psd_block));
			}
		}
	}
}




/*
void SDP_Solver::save_solution(
  const SDP_Solver_Terminate_Reason terminate_reason,
  const std::pair<std::string, Timer> &timer_pair,
  const boost::filesystem::path &out_directory,
  const Write_Solution &write_solution,
  const std::vector<size_t> &block_indices, const Verbosity &verbosity) const
{
  // Internally, El::Print() sync's everything to the root core and
  // outputs it from there.  So do not actually open the file on
  // anything but the root node.

  boost::filesystem::ofstream out_stream;
  if(El::mpi::Rank() == 0)
    {
      if(verbosity >= Verbosity::regular)
        {
          std::cout << "Saving solution to      : " << out_directory << '\n';
        }
      const boost::filesystem::path output_path(out_directory / "out.txt");
      out_stream.open(output_path);
      set_stream_precision(out_stream);
      out_stream << "terminateReason = \"" << terminate_reason << "\";\n"
                 << "primalObjective = " << primal_objective << ";\n"
                 << "dualObjective   = " << dual_objective << ";\n"
                 << "dualityGap      = " << duality_gap << ";\n"
                 << "primalError     = " << primal_error() << ";\n"
                 << "dualError       = " << dual_error << ";\n"
                 << std::setw(16) << std::left << timer_pair.first << "= "
                 << timer_pair.second.elapsed_seconds() << ";\n";
      if(!out_stream.good())
        {
          throw std::runtime_error("Error when writing to: "
                                   + output_path.string());
        }
    }
  // y is duplicated among cores, so only need to print out copy on
  // the root node.
  if(write_solution.vector_y && !dy_backup.blocks.empty())
    {
      const boost::filesystem::path y_path(out_directory / "y.txt");
      boost::filesystem::ofstream y_stream;
      if(El::mpi::Rank() == 0)
        {
          y_stream.open(y_path);
        }
      El::Print(dy_backup.blocks.at(0),
                std::to_string(dy_backup.blocks.at(0).Height()) + " "
                  + std::to_string(dy_backup.blocks.at(0).Width()),
                "\n", y_stream);
      if(El::mpi::Rank() == 0)
        {
          y_stream << "\n";
          if(!y_stream.good())
            {
              throw std::runtime_error("Error when writing to: "
                                       + y_path.string());
            }
        }
    }

  for(size_t block = 0; block != dx_backup.blocks.size(); ++block)
    {
      size_t block_index(block_indices.at(block));
      if(write_solution.vector_x)
        {
          const boost::filesystem::path x_path(
            out_directory / ("x_" + std::to_string(block_index) + ".txt"));
          boost::filesystem::ofstream x_stream;
          if(dx_backup.blocks.at(block).DistRank() == dx_backup.blocks.at(block).Root())
            {
              x_stream.open(x_path);
            }
          El::Print(dx_backup.blocks.at(block),
                    std::to_string(dx_backup.blocks.at(block).Height()) + " "
                      + std::to_string(dx_backup.blocks.at(block).Width()),
                    "\n", x_stream);
          if(dx_backup.blocks.at(block).DistRank() == dx_backup.blocks.at(block).Root())
            {
              x_stream << "\n";
              if(!x_stream.good())
                {
                  throw std::runtime_error("Error when writing to: "
                                           + x_path.string());
                }
            }
        }
      for(size_t psd_block(0); psd_block < 2; ++psd_block)
        {
          std::string suffix(std::to_string(2 * block_index + psd_block)
                             + ".txt");

          if(write_solution.matrix_X
             && dX_backup.blocks.at(2 * block + psd_block).Height() != 0)
            {
              write_psd_block(out_directory / ("X_matrix_" + suffix),
				  dX_backup.blocks.at(2 * block + psd_block));
            }
          if(write_solution.matrix_Y
             && dY_backup.blocks.at(2 * block + psd_block).Height() != 0)
            {
              write_psd_block(out_directory / ("Y_matrix_" + suffix),
				  dY_backup.blocks.at(2 * block + psd_block));
            }
        }
    }
}
*/