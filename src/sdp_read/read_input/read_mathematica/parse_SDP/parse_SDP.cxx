#include "parse_vector.hxx"
#include "parse_generic.hxx"
#include "../../../Positive_Matrix_With_Prefactor.hxx"

#include <algorithm>
#include <iterator>
#include <string>

const char *
parse_matrices(const char *begin, const char *end, const int &rank,
               const int &num_procs, const size_t &num_matrices,
               std::vector<Positive_Matrix_With_Prefactor> &matrices);

const char *parse_SDP(const char *begin, const char *end,
                      std::vector<El::BigFloat> &objectives,
                      std::vector<El::BigFloat> &normalization,
                      std::vector<Positive_Matrix_With_Prefactor> &matrices,
                      size_t &num_processed)
{
  const std::string SDP_literal("SDP[");
  auto SDP_start(
    std::search(begin, end, SDP_literal.begin(), SDP_literal.end()));
  if(SDP_start == end)
    {
      throw std::runtime_error("Could not find 'SDP['");
    }
  if(SDP_start != begin)
    {
      const char previous_char(*(SDP_start - 1));
      if(previous_char != ' ' && previous_char != '\t' && previous_char != '\n'
         && previous_char != '\r' && previous_char != ')')
        {
          throw std::runtime_error("Invalid sequence: '"
                                   + (previous_char + SDP_literal)
                                   + "'.  Is this an SDP file?");
        }
    }

  std::vector<El::BigFloat> temp_vector;
  const char *end_objective(parse_vector(SDP_start, end, temp_vector));
  if(!temp_vector.empty())
    {
      std::swap(objectives, temp_vector);
      temp_vector.clear();
    }

  const char *comma(std::find(end_objective, end, ','));
  if(comma == end)
    {
      throw std::runtime_error("Missing comma after objective");
    }

  parse_vector(std::next(comma), end, temp_vector);
  if(!temp_vector.empty())
    {
      std::swap(normalization, temp_vector);
      temp_vector.clear();
    }

  comma = std::find(end_objective, end, ',');
  if(comma == end)
    {
      throw std::runtime_error("Missing comma after normalization");
    }

  std::vector<Positive_Matrix_With_Prefactor> temp_matrices;
  const int rank(El::mpi::Rank()), num_procs(El::mpi::Size());
  const char *end_matrices(parse_matrices(
    std::next(comma), end, rank, num_procs, matrices.size(), temp_matrices));
  {
    size_t offset(matrices.size());
    matrices.resize(matrices.size() + temp_matrices.size());
    for(size_t index = 0; index < temp_matrices.size(); ++index)
      {
        std::swap(matrices[offset + index], temp_matrices[index]);
        ++num_processed;
      }
  }
  const auto close_bracket(std::find(end_matrices, end, ']'));
  if(close_bracket == end)
    {
      throw std::runtime_error("Missing ']' at end of SDP");
    }
  return std::next(close_bracket);
}
