



#include <boost/filesystem.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/algorithm/string.hpp>

#include <boost/math/tools/polynomial.hpp>

#include <set>

#include <vector>
#include <string>

#include <utility>

#include <variant>

#include <El.hpp>

#include <deque>
#include <list>
#include <map>

#include "read_input/read_mathematica/parse_SDP/parse_vector.hxx"
#include "Boost_Float.hxx"
#include "../Polynomial.hxx"




using boost::lexical_cast;
using std::string;

inline El::BigFloat abs(const El::BigFloat & a)
{
	return (a > 0) ? a : -a;
}

template <typename T>
string sign_str(T const &x)
{
	return x < 0 ? "-" : "+";
}

template <typename T>
string inner_coefficient(T const &x)
{
	string result(" " + sign_str(x) + " ");
	if (abs(x) != T(1))
		result += lexical_cast<string>(abs(x));
	return result;
}

template <typename T>
string formula_format(std::vector<T> const &a)
{
	string result;
	if (a.size() == 0)
		result += lexical_cast<string>(T(0));
	else
	{
		// First one is a special case as it may need unary negate.
		unsigned i = a.size() - 1;
		if (a[i] < 0)
			result += "-";
		if (abs(a[i]) != T(1))
			result += lexical_cast<string>(abs(a[i]));

		if (i > 0)
		{
			result += "x";
			if (i > 1)
			{
				result += "^" + lexical_cast<string>(i);
				i--;
				for (; i != 1; i--)
					if (a[i])
						result += inner_coefficient(a[i]) + "x^" + lexical_cast<string>(i);

				if (a[i])
					result += inner_coefficient(a[i]) + "x";
			}
			i--;

			if (a[i])
				result += " " + sign_str(a[i]) + " " + lexical_cast<string>(abs(a[i]));
		}
	}
	return result;
} // string formula_format(polynomial<T> const &a)


namespace boost {
	namespace serialization {
		template<class Archive>
		void save(Archive& ar, El::BigFloat const & f, const boost::serialization::version_type&) {
			std::vector<uint8_t> local_array(f.SerializedSize());
			f.Serialize(local_array.data());
			ar & local_array;
		}

		template<class Archive>
		void load(Archive& ar, El::BigFloat & f, const boost::serialization::version_type&) {
			std::vector<uint8_t> local_array(f.SerializedSize());
			ar & local_array;
			f.Deserialize(local_array.data());
		}

	}  // namespace serialization
}  // namespace boost

BOOST_SERIALIZATION_SPLIT_FREE(El::BigFloat)

static auto const boost_archive_flags = boost::archive::no_header | boost::archive::no_tracking;

void test_load()
{
	return;

	boost::filesystem::ifstream ifs("saveblock_test2/A-L0.block");

	std::vector<std::vector<std::vector<El::BigFloat>>> zzb_derivs_conv_El;

	boost::archive::binary_iarchive ia(ifs, boost_archive_flags);
	ia & zzb_derivs_conv_El;

	std::cout << "CBTab[0][0]=" << formula_format(zzb_derivs_conv_El[0][0]) << "\n";

	El::BigFloat El_f = zzb_derivs_conv_El[0][0][0];
	std::cout << "test El_f = " << El_f << " ; prec = " << El_f.Precision() << "; data_size=" << El_f.SerializedSize() << "\n";

	std::cout << "_mp_size = " << El_f.gmp_float.get_mpf_t()[0]._mp_size << " ; sizeof(mp_exp_t) = " << sizeof(mp_exp_t) << "; sizeof(mp_limb_t)=" << sizeof(mp_limb_t)
		<< " ; sizeof(int)=" << sizeof(int) << " ; El::gmp::num_limbs=" << El::gmp::num_limbs << "\n";

	std::cout.precision(200);
	std::cout << std::fixed;
	std::cout << "test El_f = " << El_f << "\n";

	return;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////

#define MMA_PARSER_ERROR(flow) {std::stringstream ss; ss << "file :" << __FILE__ << " line : " << __LINE__ << " :\n" << flow; throw std::runtime_error(ss.str());}


void load_block_folder(const std::string & block_folder, std::map<std::pair<std::string, int>, std::vector<std::vector<std::vector<El::BigFloat>>> > & blockF)
{
	namespace fs = boost::filesystem;

	std::cout << "scan blocks...\n";

	std::vector<std::vector<std::vector<El::BigFloat>>> zzb_derivs_conv_El;

	for (auto const & file : fs::recursive_directory_iterator(block_folder))
	{
		if (fs::is_regular_file(file) && file.path().extension() == std::string(".block"))
		{
			const std::string filename = file.path().filename().string();
			size_t barL = filename.find("-L");
			if (barL == std::string::npos) MMA_PARSER_ERROR("Load block error : invalid block file name : " << filename << "\n");
			const std::string stamp = filename.substr(0, barL);
			barL += 2;
			size_t dot = filename.find(".", barL);
			if (dot == std::string::npos) MMA_PARSER_ERROR("Load block error : invalid block file name : " << filename << "\n");
			int spin = std::stoi(filename.substr(barL, dot));

			std::cout << "load block with stamp = " << stamp << " spin = " << spin << "\n";

			boost::filesystem::ifstream ifs(file);
			boost::archive::binary_iarchive ia(ifs, boost_archive_flags);
			ia & zzb_derivs_conv_El;

			blockF.emplace(std::make_pair(stamp, spin), zzb_derivs_conv_El);
		}
	}

	std::cout << "blocks loaded\n";

	//for (const auto&[key, value] : blockF)
	//	std::cout << "find block with stamp=" << key.first << " spin=" << key.second << " max_m=" << value.size() - 1 << "\n";

}


void load_block_folder(const std::string & block_folder, const std::string & stamp, int spin,
	std::map<std::pair<std::string, int>, std::vector<std::vector<std::vector<El::BigFloat>>> > & blockF)
{
	namespace fs = boost::filesystem;

	std::vector<std::vector<std::vector<El::BigFloat>>> zzb_derivs_conv_El;

	auto pblock = blockF.find(std::make_pair(stamp, spin));
	if (pblock != blockF.end())
		MMA_PARSER_ERROR("load_block_folder error : " << stamp << "-L" << spin << ".block" << " already exist.\n");

	std::stringstream path_str;
	path_str << block_folder << "/" << stamp << "-L" << spin << ".block";

	const boost::filesystem::path file(path_str.str());

	if(!fs::is_regular_file(file)) //exists(file) && 
		MMA_PARSER_ERROR("load_block_folder error : " << stamp << "-L" << spin << ".block" << " is missing.\n");

	boost::filesystem::ifstream ifs(file);
	boost::archive::binary_iarchive ia(ifs, boost_archive_flags);
	ia & zzb_derivs_conv_El;

	blockF.emplace(std::make_pair(stamp, spin), zzb_derivs_conv_El);

	std::cout << "Rank=" << El::mpi::Rank() << " : load block with stamp = " << stamp << " spin = " << spin << "\n";
}



//////////////////////////////////////////////////////////////////////////////////////////////////////

//std::string test_str = " 1.2+3 (2-7) + 2 3 - (2.4 * 7)/(-3+2/3*4) + Fs[stamp, 1, 1, 2, 0.134, 3, 17] , 2+3  ";
//std::string test_str = " -1.2+3/2*3-1.1*1.2 3.4 + 4 , 1+3";

//std::string test_str = "-2.3 + 2.3 Fs[stamp, 1, 1, 2, 0.134, 3, 17] , 2+3  ";
//std::string test_str = "-2.3 + 2.3";



/*

template <typename T>
const char *
parse_MMA_expr(const char *begin, const char *end, T & result)
{
	const std::vector<char> delimiters({ ',','}' });

	const auto delimiter(
		std::find_first_of(begin, end, delimiters.begin(), delimiters.end()));
	if (delimiter == end)
	{
		throw std::runtime_error("Missing '}' at end of array of polynomials");
	}
}




{
	const char * pstr = begin;
	int i = 1;
	while (pstr != end && i <= 5000)
	{
		MMA_ELEMENT element;

		pstr = parse_MMA_element(pstr, end, element);
		std::cout << "Find element #" << i << " : " << element << "\n";

		i++;
	}
	exit(0);
	return pstr;
}

{
	const char * pstr = begin;
	int i = 1;
	while (pstr != end && i <= 5000)
	{
		MMA_TOKEN word;

		pstr = parse_get_token(pstr, end, word);
		std::cout << "Find word #" << i << " : " << word << "\n";

		i++;
	}
	exit(0);
	return pstr;
}

exit(0);
return nullptr;
*/