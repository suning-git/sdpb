#include "../Block_Diagonal_Matrix.hxx"
#include "../Block_Info.hxx"
#include "../../Block_Vector.hxx"


void write_psd_block(const boost::filesystem::path &outfile,
	const El::DistMatrix<El::BigFloat> &block)
{
	boost::filesystem::ofstream stream;
	if (block.DistRank() == block.Root())
	{
		stream.open(outfile);
	}
	El::Print(block,
		std::to_string(block.Height()) + " "
		+ std::to_string(block.Width()),
		"\n", stream);
	if (block.DistRank() == block.Root())
	{
		stream << "\n";
		if (!stream.good())
		{
			throw std::runtime_error("Error when writing to: "
				+ outfile.string());
		}
	}
}


std::ostream &operator<<(std::ostream &os, const Block_Diagonal_Matrix &A)
{
  os << "{";
  for(auto block(A.blocks.begin()); block != A.blocks.end();)
    {
      El::Print(*block, "", ",", os);
      ++block;
      if(block != A.blocks.end())
        {
          os << ", ";
        }
    }
  os << "}";
  return os;
}


void Block_Diagonal_Matrix::print(const Block_Info &block_info, const std::string file_path_no_suffix)
{
	for (size_t block = 0; block != block_info.block_indices.size(); ++block)
	{
		size_t block_index(block_info.block_indices.at(block));
		for (size_t psd_block(0); psd_block < 2; ++psd_block)
		{
			int globalID = 2 * block_index + psd_block;
			int localID = 2 * block + psd_block;

			auto &curblock = this->blocks.at(localID);

			write_psd_block(file_path_no_suffix + "_" + std::to_string(globalID) + ".txt", curblock);
		}
	}
}

void Block_Diagonal_Matrix::print(const Block_Info &block_info, size_t globalID_to_print, const std::string file_path_no_suffix)
{
	for (size_t block = 0; block != block_info.block_indices.size(); ++block)
	{
		size_t block_index(block_info.block_indices.at(block));
		for (size_t psd_block(0); psd_block < 2; ++psd_block)
		{
			int globalID = 2 * block_index + psd_block;
			int localID = 2 * block + psd_block;

			auto &curblock = this->blocks.at(localID);

			if (globalID == globalID_to_print)
			{
				write_psd_block(file_path_no_suffix + "_" + std::to_string(globalID) + ".txt", curblock);
			}
		}
	}
}


void Block_Vector::print(const Block_Info &block_info, size_t globalID_to_print, const std::string file_path_no_suffix)
{
	for (size_t localID = 0; localID != block_info.block_indices.size(); ++localID)
	{
		size_t globalID(block_info.block_indices.at(localID));

		auto &curblock = this->blocks.at(localID);

		if (globalID == globalID_to_print)
		{
			write_psd_block(file_path_no_suffix + "_" + std::to_string(globalID) + ".txt", curblock);
		}
	}
}


void Block_Vector::print(const Block_Info &block_info, const std::string file_path_no_suffix)
{
	for (size_t localID = 0; localID != block_info.block_indices.size(); ++localID)
	{
		size_t globalID(block_info.block_indices.at(localID));

		auto &curblock = this->blocks.at(localID);

		write_psd_block(file_path_no_suffix + "_" + std::to_string(globalID) + ".txt", curblock);
	}
}

