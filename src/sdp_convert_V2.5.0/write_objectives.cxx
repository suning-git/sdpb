#include "write_vector.hxx"
#include "byte_counter.hxx"
#include "../set_stream_precision.hxx"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

size_t write_objectives_V250(const boost::filesystem::path &output_dir,
                        const El::BigFloat &objective_const,
                        const std::vector<El::BigFloat> &dual_objective_b)
{
  const boost::filesystem::path output_path(output_dir / "objectives.json");
  byte_counter counter;
  {
    boost::iostreams::filtering_ostream output_stream;
    output_stream.push(boost::ref(counter));
    // Use gzip with no compression to get a CRC
    output_stream.push(
      boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(0)));
    output_stream.push(boost::iostreams::file_sink(output_path.string()));
    set_stream_precision(output_stream);
    output_stream << "{\n  \"constant\": \"" << objective_const << "\",\n";
    write_vector(output_stream, dual_objective_b, "  ", "b");
    output_stream << "\n}\n";
    if(!output_stream.good())
      {
        throw std::runtime_error("Error when writing to: "
                                 + output_path.string());
      }
  }
  return counter.num_bytes;
}
