#include "parse_vector.hxx"
#include "parse_generic.hxx"
#include "../../../Positive_Matrix_With_Prefactor.hxx"

#include <algorithm>
#include <iterator>
#include <string>

#include <vector>
#include <chrono>
#include <thread>


void mpi_counter_init(int init_value = 0);
int mpi_counter_get();


const char * findClosingParen(const char *paren_begin)
{
	char char_closing, char_opening = *paren_begin;
	switch (char_opening)
	{
	case '(':
		char_closing = ')';
		break;
	case '[':
		char_closing = ']';
		break;
	case '{':
		char_closing = '}';
		break;
	default:
		throw std::runtime_error("findClosingParen incorrect paren_begin");
	}

	const char * paren_end = paren_begin + 1;
	int counter = 1;
	while (1)
	{
		char c = *paren_end;
		if (c == char_opening) {
			counter++;
		}
		else if (c == char_closing) {
			if ((--counter) == 0)break;
		}

		if (c == 0)
			throw std::runtime_error("findClosingParen unbalanced parentheses");

		paren_end++;
	}
	return paren_end;
}


const char *parse_matrix_bypass(const char *begin, const char *end, Positive_Matrix_With_Prefactor &matrix)
{
	const std::string matrix_literal("PositiveMatrixWithPrefactor[");
	auto matrix_start(
		std::search(begin, end, matrix_literal.begin(), matrix_literal.end()));
	if (matrix_start == end)
	{
		throw std::runtime_error("Could not find '" + matrix_literal + "'");
	}

	const char *paren_begin = matrix_start + matrix_literal.size() - 1;

	//std::cout << "[Rank=" << El::mpi::Rank() << "] " << "findClosingParen paren_begin=" << std::string(paren_begin, 20) << "\n";
	const char *paren_end = findClosingParen(paren_begin);
	//std::cout << "[Rank=" << El::mpi::Rank() << "] " << "findClosingParen end\n";
	return paren_end + 1;
}


int get_next_matrix_index(const int & num_procs, const int & rank, int current_matrix_index = -1)
{
	if (current_matrix_index == -1)
	{
		return rank % num_procs;
	}
	else
	{
		return current_matrix_index + num_procs;
	}
}


const char *
parse_matrices(const char *begin, const char *end, const int &rank,
	const int &num_procs, const size_t &num_matrices,
	std::vector<Positive_Matrix_With_Prefactor> &matrices,
	std::vector<int> &matrices_valid_indices)
{
	const auto open_brace(std::find(begin, end, '{'));
	if (open_brace == end)
	{
		throw std::runtime_error("Could not find '{' to start array");
	}

	auto delimiter(open_brace);
	const std::vector<char> delimiters({ ',', '}' });
	int matrix_index(num_matrices);

	//int current_matrix_index = mpi_counter_get();
	int current_matrix_index = get_next_matrix_index(num_procs, rank);
	do
	{
		//std::cout << "[Rank=" << El::mpi::Rank() << " BEGIN] " << "looking for current_matrix_index=" << current_matrix_index << ", while matrix_index=" << matrix_index << "\n";

		auto start_matrix(std::next(delimiter));
		Positive_Matrix_With_Prefactor matrix;

		const char *end_matrix;
		if (matrix_index == current_matrix_index)
		{
			end_matrix = parse_generic(start_matrix, end, matrix);
		}
		else
		{
			end_matrix = parse_matrix_bypass(start_matrix, end, matrix);
		}
		//std::cout << "[Rank=" << El::mpi::Rank() << " parsed] " << "matrix_index=" << matrix_index << "\n";

		matrices.emplace_back();
		//if(matrix_index % num_procs == rank) swap(matrices.back(), matrix);
		if (matrix_index == current_matrix_index)// && matrix_index != 60
		{
			swap(matrices.back(), matrix);
			matrices_valid_indices.push_back(current_matrix_index);
			//current_matrix_index = mpi_counter_get();
			current_matrix_index = get_next_matrix_index(num_procs, rank, current_matrix_index);
		}
		//std::cout << "[Rank=" << El::mpi::Rank() << "   END] " << "end_matrix=" << std::string(end_matrix, 20) << "\n";

		++matrix_index;

		delimiter = std::find_first_of(end_matrix, end, delimiters.begin(),
			delimiters.end());
		if (delimiter == end)
		{
			throw std::runtime_error(
				"Missing '}' at end of array of PositiveMatrixWithPrefactor");
		}
	} while (*delimiter != '}');

	std::cout << "[Rank=" << El::mpi::Rank() << "] " << "matrices_valid_indices=";
	for (auto i : matrices_valid_indices)std::cout << i << " ";
	std::cout << "\n";

	//std::cout << "[Rank=" << El::mpi::Rank() << "] " << "enter MPI_Barrier\n";
	//El::mpi::Barrier(MPI_COMM_WORLD);
	//std::cout << "[Rank=" << El::mpi::Rank() << "] " << "finished\n";

	return std::next(delimiter);
}
