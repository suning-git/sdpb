

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/filesystem/fstream.hpp>
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

#include "Positive_Matrix_With_Prefactor.hxx"
#include "read_input/read_mathematica/parse_SDP/parse_vector.hxx"
#include "Boost_Float.hxx"
#include "../Polynomial.hxx"


////////////////////////////////////////  infarestructure //////////////////////////////////////

inline El::BigFloat as_BigFloat(const Boost_Float & bf)
{
	mpf_t mpf_v;
	mpf_init(mpf_v);

	mpfr_get_f(mpf_v, bf.backend().data(), MPFR_RNDN);

	El::BigFloat result;
	result.gmp_float = mpf_class(mpf_v);

	return result;
}

inline Boost_Float as_BoostFloat(const El::BigFloat & ef)
{
	return Boost_Float(ef.gmp_float.get_mpf_t());
}



////////////////////////////////////////   parse Mathematica expression //////////////////////////////////////


std::string parse_number(const char *begin, const char *end);

//////////////

template <class InputIterator>
inline InputIterator find_delimiters(InputIterator begin, InputIterator end, const std::string delimiters)
{
	return std::find_first_of(begin, end, delimiters.begin(), delimiters.end());
};

const std::string MMA_expr_delimiters("()[]{}+-*/^, \t\n\v\f\r");
inline const char * find_next_MMA_expr_delimiters(const char *begin, const char *end)
{
	const char * pstr = find_delimiters(begin, end, MMA_expr_delimiters);
	if (*pstr == '\n' && pstr > begin && pstr < end && *(pstr - 1) == '\\') // special case where a line end with "\\\n"
		return find_next_MMA_expr_delimiters(pstr + 1, end);
	else if (*pstr == '\r' && pstr > begin && pstr < end - 1 && *(pstr - 1) == '\\' && *(pstr + 1) == '\n') // special case where a line end with "\\\r\n"
		return find_next_MMA_expr_delimiters(pstr + 2, end);
	else
		return pstr;
}

template <typename T>
inline bool string_containQ(const std::string & str, const T substr)
{
	return str.find(substr) != std::string::npos;
};

inline bool string_contain_delimiters_Q(const std::string & str, const std::string delimiters)
{
	return find_delimiters(str.begin(), str.end(), delimiters) != str.end();
};

inline bool string_contain_MMA_delimiters_Q(const std::string & str)
{
	return find_delimiters(str.begin(), str.end(), MMA_expr_delimiters) != str.end();
};

using MMA_TOKEN = std::variant<std::monostate, int, El::BigFloat, char, std::string, std::string>;
#define MMA_TOKEN_Invalid	0
#define MMA_TOKEN_Integer	1
#define MMA_TOKEN_Real		2
#define MMA_TOKEN_Operator	3
#define MMA_TOKEN_Symbol	4
#define MMA_TOKEN_String	5

#define AS_MMA_TOKEN(var,T) std::get<MMA_TOKEN_##T>(var)
#define SET_MMA_TOKEN(var, value,T) var.emplace<MMA_TOKEN_##T>(value)


const char * ptr_MMA_begin;
const char * ptr_MMA_current;

#define MMA_PARSER_ERROR(flow) {std::stringstream ss; ss << "current ptr in MMA file : " << ptr_MMA_current - ptr_MMA_begin << "\n" << "file :" << __FILE__ << " line : " << __LINE__ << " :\n" << flow; throw std::runtime_error(ss.str());}


std::ostream& operator<<(std::ostream& os, const MMA_TOKEN & v)
{
	switch (v.index())
	{
	case MMA_TOKEN_Invalid:
		os << "[invalid]";
		break;
	case MMA_TOKEN_Integer:
		os << "[Integer " << AS_MMA_TOKEN(v, Integer) << "]";
		break;
	case MMA_TOKEN_Real:
		os << "[Float " << AS_MMA_TOKEN(v, Real) << "]";
		break;
	case MMA_TOKEN_Operator:
		os << "[Operator " << AS_MMA_TOKEN(v, Operator) << "]";
		break;
	case MMA_TOKEN_Symbol:
		os << "[Symbol " << AS_MMA_TOKEN(v, Symbol) << "]";
		break;
	case MMA_TOKEN_String:
		os << "[String " << AS_MMA_TOKEN(v, String) << "]";
		break;
	default:
		break;
	}

	/*
	std::visit([&os](auto&& arg) {
		os << arg;
	}, v);*/
	return os;
}

const char *skip_space_from_left(const char *b, const char *e)
{
	const char * rslt = std::find_if_not(b, e, [](const char pc) {return std::isspace(pc); });
	if (*rslt == '\\' && *(rslt + 1) == '\n')
		return skip_space_from_left(rslt + 2, e);
	else
		return rslt;
};


// example : 1.`200.*^-50  or 1.`200.*^50
// p is a digit
const char * find_end_of_number(const char *p, const char *end)
{
	auto q = find_next_MMA_expr_delimiters(p + 1, end);

	if (*q != '*') return q;
	auto q2 = skip_space_from_left(q + 1, end);
	if (*q2 != '^') return q;
	q = q2;

	// the text has to be in scientific notation.
	q = skip_space_from_left(q + 1, end);
	if (*q == '+' || *q == '-') q = skip_space_from_left(q + 1, end);

	q = find_next_MMA_expr_delimiters(q, end);
	return q;
}

const char * parse_get_token(const char *begin, const char *end, MMA_TOKEN & token);
const char * parse_get_token_func(const char *begin, const char *end, MMA_TOKEN & token)
{
	ptr_MMA_current = begin;

	auto skip_space_l = [](const char *b, const char *e)
	{
		return std::find_if_not(b, e,[](const char pc) {return std::isspace(pc); });
	};

	const char * p = skip_space_from_left(begin, end); //skip_space_l(begin,end);
	
	if (p == end) { token = std::monostate{}; return p; } 	// " 2+(3+1)   "

	// parse basic operators
	const auto ops_ch = std::string("()[]{}+-*/^,");
	if (string_containQ(ops_ch,*p))
	{
		SET_MMA_TOKEN(token, *p, Operator);
		return p + 1;
	}

	// parse numbers
	if (isdigit(*p))
	{
		// this is old code with a bug : it can't handle the following case : 1.`200.*^\\\n-5
		/* 
		auto q = find_next_MMA_expr_delimiters(p + 1, end);

		if (q + 3 < end && *q == '*' && *(q + 1) == '^')  // 1.`200.*^-50  or 1.`200.*^50
		{
			if (*(q + 2) == '+' || *(q + 2) == '-')
				q = find_next_MMA_expr_delimiters(q + 3, end);
			else
				q = find_next_MMA_expr_delimiters(q + 2, end);
		}
		*/

		auto q = find_end_of_number(p, end);

		if (std::find_if_not(p, end, isdigit) == q) // integer
		{
			SET_MMA_TOKEN(token, std::stoi(std::string(p, q)), Integer);
		}
		else
		{
			//std::cout << p - end << " : token read : real number raw_string =" << std::string(p, q) << "\ncleaned_string = " << parse_number(p, q) << "\n";
			SET_MMA_TOKEN(token, El::BigFloat(parse_number(p, q)), Real);
		}
		
		return q;
	}

	// parse symbol
	if (isalpha(*p))
	{
		auto q = find_next_MMA_expr_delimiters(p + 1, end);
		SET_MMA_TOKEN(token, std::string(p, q), Symbol);
		return q;
	}

	// parse string
	if (*p == '\"')
	{
		auto q = std::find(p + 1, end, '\"');
		if (q == end)
		{
			std::cout << "Unrecognizable expression :" << std::string(p, q) << "\n";
			exit(0);
		}
		SET_MMA_TOKEN(token, std::string(p + 1, q), String);
		return q + 1;
	}

	// handle "\\\n" case
	if (*p == '\\' && *(p + 1) == '\n')
	{
		std::cout << "the code shouldn't reach here :" << std::string(p, p + 20) << "\n";
		exit(0);
		return parse_get_token(p + 2, end, token);
	}

	std::cout << "Unrecognizable expression :" << std::string(p, p + 10) << "\n";
	exit(0);
	return "";
}

// this is a wrapper to keep track of the debug information
const char * parse_get_token(const char *begin, const char *end, MMA_TOKEN & token)
{
	const char * parse_begin = begin;
	const char * parse_end = parse_get_token_func(begin, end, token);
	//std::cout << "parse_get_token parse " << std::string(parse_begin, parse_end) << " to be " << token << "\n";
	return parse_end;
}



using MMA_EXPR = std::variant<El::BigFloat, Polynomial>;
#define MMA_EXPR_Number		0
#define MMA_EXPR_Polynomial	1
#define AS_MMA_EXPR(var,T) std::get<MMA_EXPR_##T>(var)
#define SET_MMA_EXPR(var, value,T) var.emplace<MMA_EXPR_##T>(value)

using MMA_ELEMENT = std::variant<std::monostate, MMA_EXPR, char>;
#define MMA_ELEMENT_Invalid		0
#define MMA_ELEMENT_Expression	1
#define MMA_ELEMENT_Operator	2
#define AS_MMA_ELEMENT(var,T) std::get<MMA_ELEMENT_##T>(var)
#define AS_MMA_ELEMENT_Number(var) std::get<MMA_EXPR_Number>(std::get<MMA_ELEMENT_Expression>(var))
#define AS_MMA_ELEMENT_Polynomial(var) std::get<MMA_EXPR_Polynomial>(std::get<MMA_ELEMENT_Expression>(var))
#define SET_MMA_ELEMENT(var, value,T) var.emplace<MMA_ELEMENT_##T>(value)


std::ostream& operator<<(std::ostream& os, const MMA_ELEMENT & v)
{
	switch (v.index())
	{
	case MMA_ELEMENT_Invalid:
		os << "[invalid]";
		break;
	case MMA_ELEMENT_Expression:
	{
		const MMA_EXPR & elmt = AS_MMA_ELEMENT(v, Expression);
		switch (elmt.index())
		{
			//case MMA_EXPR_Invalid:
			//	os << "[invalid]";
			//	break;
		case MMA_EXPR_Number:
		{
			int prev_prec = std::cout.precision();
			os << "[Number " << std::setprecision(15) << AS_MMA_EXPR(elmt, Number) << std::setprecision(prev_prec) << "]";
		}
		break;
		case MMA_EXPR_Polynomial:
			os << "[Polynomial " << AS_MMA_EXPR(elmt, Polynomial) << "]";
			break;
		}
		break;
	}
	case MMA_ELEMENT_Operator:
		os << "[Operator " << AS_MMA_ELEMENT(v, Operator) << "]";
		break;
	}
	return os;
}

const char * parse_MMA_expr(const char *begin, const char *end, MMA_ELEMENT & result);
void parse_MMA_symbol(const std::string & name, MMA_ELEMENT & result);
const char * parse_MMA_function(const std::string & name, const char * begin, const char * end, MMA_ELEMENT & result);

const char * parse_MMA_element(const char *begin, const char *end, MMA_ELEMENT & result)
{

	MMA_TOKEN token;
	const char * pstr = parse_get_token(begin, end, token);

	switch (token.index())
	{
	case MMA_TOKEN_Invalid:
		result = std::monostate{};
		return pstr;

	case MMA_TOKEN_Integer:
		SET_MMA_ELEMENT(result, AS_MMA_TOKEN(token, Integer), Expression);
		return pstr;

	case MMA_TOKEN_Real:
		SET_MMA_ELEMENT(result, AS_MMA_TOKEN(token, Real), Expression);
		return pstr;

	case MMA_TOKEN_Operator:
	{
		char op = AS_MMA_TOKEN(token, Operator);
		switch (op)
		{
		case '+':
		case '-':
		case '*':
		case '/':
		case '^':
			SET_MMA_ELEMENT(result, op, Operator);
			return pstr;

		case '(':
		{
			pstr = parse_MMA_expr(pstr, end, result);

			MMA_TOKEN token_next;
			pstr = parse_get_token(pstr, end, token_next);

			if (token_next.index() != MMA_TOKEN_Operator || AS_MMA_TOKEN(token_next, Operator) != ')')
				MMA_PARSER_ERROR("Expecting ')' at " << std::string(pstr-1, pstr + 3) << "\nBut I got " << token_next << "\n");

			return pstr;
		}

		case ',': // end of expression
		case ')':
		case ']':
		case '}':
			//SET_MMA_ELEMENT(result, op, Operator);
			//return pstr;
			result = std::monostate{};
			pstr--;
			return pstr;

		default:
			MMA_PARSER_ERROR("Un expected operator : " << op << " before " << std::string(pstr, pstr + 10) << "\n");
			break;
		}
	}
	break;
	case MMA_TOKEN_Symbol:
		if (*pstr == '[')
		{
			pstr = parse_MMA_function(AS_MMA_TOKEN(token, Symbol), pstr + 1, end, result);
		}
		else
		{
			parse_MMA_symbol(AS_MMA_TOKEN(token, Symbol), result);
		}
		return pstr;
		break;

	case MMA_TOKEN_String:
	default:
		MMA_PARSER_ERROR("Un expected token : " << token << " before " << std::string(pstr, pstr + 10) << "\n");
		break;
	}
	return begin;
}

inline const char * parse_MMA_element_as_expression(const char *begin, const char *end, MMA_ELEMENT & result)
{
	const char * pstr = parse_MMA_element(begin, end, result);

	if (result.index() != MMA_ELEMENT_Expression)
		MMA_PARSER_ERROR("Expecting expression at : " << std::string(begin, 10) << ".\nBut I got " << result << "\n");
	return pstr;
}

const char * parse_MMA_expr_as_number(const char *begin, const char *end, El::BigFloat & result)
{
	MMA_ELEMENT element;
	const char * pstr = parse_MMA_expr(begin, end, element);

	if (element.index() != MMA_ELEMENT_Expression || AS_MMA_ELEMENT(element, Expression).index() != MMA_EXPR_Number)
		MMA_PARSER_ERROR("Expecting number at : " << std::string(begin, 10) << ".\nBut I got " << element << "\n");

	result = AS_MMA_ELEMENT_Number(element);
	return pstr;
}


inline bool parse_MMA_expr_delimitersQ(MMA_ELEMENT & elmt)
{
	if (elmt.index() == MMA_ELEMENT_Invalid) return true;
	if (elmt.index() != MMA_ELEMENT_Operator) return false;
	const char op = AS_MMA_ELEMENT(elmt, Operator);
	switch (op)
	{
	case ']':
	case ',':
		return true;
	default:
		return false;
	}
}

const char * parse_MMA_expr_list(const char *begin, const char *end, std::list<MMA_ELEMENT> & chain)
{
	MMA_ELEMENT elmt;
	const char * pstr = parse_MMA_element(begin, end, elmt);
	if (elmt.index() == MMA_ELEMENT_Operator)
	{
		char op = AS_MMA_ELEMENT(elmt, Operator);
		switch (op)
		{
		case '+':
			pstr = parse_MMA_element_as_expression(pstr, end, elmt);
			chain.push_back(std::move(elmt));
			break;
		case '-':
		{
			pstr = parse_MMA_element_as_expression(pstr, end, elmt);
			MMA_EXPR & expr = AS_MMA_ELEMENT(elmt, Expression);
			if (expr.index()==MMA_EXPR_Number)
				AS_MMA_ELEMENT_Number(elmt) *= -1;
			else
				AS_MMA_ELEMENT_Polynomial(elmt) *= -1;
			chain.push_back(std::move(elmt));
			break;
		}
		default:
			MMA_PARSER_ERROR("parse_MMA_expr error 1st op = " << op << "\n");
			break;
		}
	}
	else
	{
		chain.push_back(std::move(elmt));
	}

	while (pstr != end)
	{
		//const char * pstr_next = parse_MMA_element(pstr, end, elmt);
		//if (parse_MMA_expr_delimitersQ(elmt)) break;
		//pstr = pstr_next;

		pstr = parse_MMA_element(pstr, end, elmt);

		if (parse_MMA_expr_delimitersQ(elmt)) break;

		if (chain.back().index() == MMA_ELEMENT_Operator && elmt.index() == MMA_ELEMENT_Expression)
		{
			chain.emplace_back(std::move(elmt));
			continue;
		}

		if (chain.back().index() == MMA_ELEMENT_Expression && elmt.index() == MMA_ELEMENT_Operator)
		{
			chain.emplace_back(std::move(elmt));
			continue;
		}

		if (chain.back().index() == MMA_ELEMENT_Expression && elmt.index() == MMA_ELEMENT_Expression)
		{
			chain.emplace_back('*');
			chain.emplace_back(std::move(elmt));
			continue;
		}

		MMA_PARSER_ERROR("parse_MMA_expr_Times error : illegal element : " << elmt << " after " << chain.back() << " before " << std::string(pstr, 10) << "\n");
	}

	return pstr;
}

inline int parse_MMA_precedence(const char op)
{
	switch (op)
	{
	case '+':
	case '-':
		return 3;
	case '*':
	case '/':
		return 2;
	case '^':
		return 1;
	default:
		MMA_PARSER_ERROR("parse_MMA_precedence unexpected error : illegal op : " << op << "\n");
		break;
	}
}
inline bool parse_MMA_precedence_orderedQ(MMA_ELEMENT & op1, MMA_ELEMENT & op2)
{
	return parse_MMA_precedence(AS_MMA_ELEMENT(op1, Operator)) <= parse_MMA_precedence(AS_MMA_ELEMENT(op2, Operator));
}

inline void parse_MMA_expr_add(MMA_EXPR & e1, MMA_EXPR & e2)
{
	if (e1.index() == MMA_EXPR_Number && e2.index() == MMA_EXPR_Number)
	{
		AS_MMA_EXPR(e1, Number) += AS_MMA_EXPR(e2, Number);
		return;
	}

	if (e1.index() == MMA_EXPR_Polynomial && e2.index() == MMA_EXPR_Polynomial)
	{
		AS_MMA_EXPR(e1, Polynomial) += AS_MMA_EXPR(e2, Polynomial);
		return;
	}

	if (e1.index() == MMA_EXPR_Polynomial && e2.index() == MMA_EXPR_Number)
	{
		AS_MMA_EXPR(e1, Polynomial) += AS_MMA_EXPR(e2, Number);
		return;
	}

	if (e1.index() == MMA_EXPR_Number && e2.index() == MMA_EXPR_Polynomial)
	{
		AS_MMA_EXPR(e2, Polynomial) += AS_MMA_EXPR(e1, Number);
		SET_MMA_EXPR(e1, AS_MMA_EXPR(e2, Polynomial), Polynomial);
		return;
	}

	MMA_PARSER_ERROR("parse_MMA_expr_add unexpected error : e1.index() = " << e1.index() << " e2.index() = " << e2.index() << "\n");
}


inline void parse_MMA_expr_substract(MMA_EXPR & e1, MMA_EXPR & e2)
{
	if (e1.index() == MMA_EXPR_Number && e2.index() == MMA_EXPR_Number)
	{
		AS_MMA_EXPR(e1, Number) -= AS_MMA_EXPR(e2, Number);
		return;
	}

	if (e1.index() == MMA_EXPR_Polynomial && e2.index() == MMA_EXPR_Polynomial)
	{
		AS_MMA_EXPR(e1, Polynomial) -= AS_MMA_EXPR(e2, Polynomial);
		return;
	}

	if (e1.index() == MMA_EXPR_Polynomial && e2.index() == MMA_EXPR_Number)
	{
		AS_MMA_EXPR(e1, Polynomial) -= AS_MMA_EXPR(e2, Number);
		return;
	}

	if (e1.index() == MMA_EXPR_Number && e2.index() == MMA_EXPR_Polynomial)
	{
		-AS_MMA_EXPR(e2, Polynomial);
		AS_MMA_EXPR(e2, Polynomial) += AS_MMA_EXPR(e1, Number);
		SET_MMA_EXPR(e1, AS_MMA_EXPR(e2, Polynomial), Polynomial);
		return;
	}

	MMA_PARSER_ERROR("parse_MMA_expr_substract unexpected error : e1.index() = " << e1.index() << " e2.index() = " << e2.index() << "\n");
}
inline void parse_MMA_expr_multiply(MMA_EXPR & e1, MMA_EXPR & e2)
{
	if (e1.index() == MMA_EXPR_Number && e2.index() == MMA_EXPR_Number)
	{
		AS_MMA_EXPR(e1, Number) *= AS_MMA_EXPR(e2, Number);
		return;
	}

	if (e1.index() == MMA_EXPR_Polynomial && e2.index() == MMA_EXPR_Number)
	{
		AS_MMA_EXPR(e1, Polynomial) *= AS_MMA_EXPR(e2, Number);
		return;
	}

	if (e1.index() == MMA_EXPR_Number && e2.index() == MMA_EXPR_Polynomial)
	{
		AS_MMA_EXPR(e2, Polynomial) *= AS_MMA_EXPR(e1, Number);
		SET_MMA_EXPR(e1, AS_MMA_EXPR(e2, Polynomial), Polynomial);
		return;
	}

	MMA_PARSER_ERROR("parse_MMA_expr_multiply unexpected error : e1.index() = " << e1.index() << " e2.index() = " << e2.index() << "\n");
}

inline void parse_MMA_expr_divide(MMA_EXPR & e1, MMA_EXPR & e2)
{
	if (e1.index() == MMA_EXPR_Number && e2.index() == MMA_EXPR_Number)
	{
		AS_MMA_EXPR(e1, Number) /= AS_MMA_EXPR(e2, Number);
		return;
	}

	if (e1.index() == MMA_EXPR_Polynomial && e2.index() == MMA_EXPR_Number)
	{
		AS_MMA_EXPR(e1, Polynomial) /= AS_MMA_EXPR(e2, Number);
		return;
	}

	MMA_PARSER_ERROR("parse_MMA_expr_divide unexpected error : e1.index() = " << e1.index() << " e2.index() = " << e2.index() << "\n");
}

inline void parse_MMA_expr_power(MMA_EXPR & e1, MMA_EXPR & e2)
{
	if (e1.index() == MMA_EXPR_Number && e2.index() == MMA_EXPR_Number)
	{
		AS_MMA_EXPR(e1, Number) = as_BigFloat(pow(as_BoostFloat(AS_MMA_EXPR(e1, Number)), as_BoostFloat(AS_MMA_EXPR(e2, Number))));
		return;
	}

	MMA_PARSER_ERROR("parse_MMA_expr_power unexpected error : e1.index() = " << e1.index() << " e2.index() = " << e2.index() << "\n");
}

inline void parse_MMA_expr_single_operate(std::list<MMA_ELEMENT> & chain, std::list<MMA_ELEMENT>::iterator it)
{
	auto it_l = std::prev(it);
	auto it_r = std::next(it);
	switch (AS_MMA_ELEMENT(*it, Operator))
	{
	case '+':
		parse_MMA_expr_add(AS_MMA_ELEMENT(*it_l, Expression), AS_MMA_ELEMENT(*it_r, Expression));
		break;
	case '-':
		parse_MMA_expr_substract(AS_MMA_ELEMENT(*it_l, Expression), AS_MMA_ELEMENT(*it_r, Expression));
		break;
	case '*':
		parse_MMA_expr_multiply(AS_MMA_ELEMENT(*it_l, Expression), AS_MMA_ELEMENT(*it_r, Expression));
		break;
	case '/':
		parse_MMA_expr_divide(AS_MMA_ELEMENT(*it_l, Expression), AS_MMA_ELEMENT(*it_r, Expression));
		break;
	case '^':
		parse_MMA_expr_power(AS_MMA_ELEMENT(*it_l, Expression), AS_MMA_ELEMENT(*it_r, Expression));
		break;
	default:
		MMA_PARSER_ERROR("parse_MMA_precedence unexpected error : illegal op : " << AS_MMA_ELEMENT(*it, Operator) << "\n");
		break;
	}

	chain.erase(it);
	chain.erase(it_r);
}

const char * parse_MMA_expr(const char *begin, const char *end, MMA_ELEMENT & result)
{
	std::list<MMA_ELEMENT> chain;
	const char * pstr = parse_MMA_expr_list(begin, end, chain);

	//for (auto ele : chain) std::cout << "Find element : " << ele << "\n\n";

	while (chain.size()>1)
	{
		for (auto op = std::next(chain.begin()); op != chain.end(); std::advance(op, 2))
		{
			auto next_op = std::next(op, 2);
			if (next_op == chain.end() || parse_MMA_precedence_orderedQ(*op, *next_op))
			{
				//std::cout << "operate : " << *std::prev(op) << *op << *std::next(op) << "\n";
				parse_MMA_expr_single_operate(chain, op);
				//std::cout << "result  : " << *std::prev(next_op) << "\n";
				break;
			}
		}
	}

	result = chain.back();
	return pstr;
}



////////////  Symbol ///////////////////////////////


inline const char * parse_MMA_check_op(const char * begin, const char * end, MMA_TOKEN & token, char op)
{
	const char * pstr = parse_get_token(begin, end, token);
	if (token.index() != MMA_TOKEN_Operator || AS_MMA_TOKEN(token, Operator) != op)
		MMA_PARSER_ERROR("parse_MMA_check_op error : expect " << op << ", but I got " << token << " from text " << std::string(begin, 20) << "\n");
	return pstr;
}

inline const char * parse_MMA_get_op(const char * begin, const char * end, MMA_TOKEN & token, char & op)
{
	const char * pstr = parse_get_token(begin, end, token);
	if (token.index() != MMA_TOKEN_Operator)
		MMA_PARSER_ERROR("parse_MMA_get_op error : expect Operator, but I got " << token << " from text " << std::string(begin, 20) << "\n");
	op = AS_MMA_TOKEN(token, Operator);
	return pstr;
}


inline const char * parse_MMA_token_as_int(const char * begin, const char * end, MMA_TOKEN & token, int & intnum)
{
	const char * pstr = parse_get_token(begin, end, token);
	if (token.index() != MMA_TOKEN_Integer)
		MMA_PARSER_ERROR("parse_MMA_get_token_as_int error : expect integer, but I got " << token << " from text " << std::string(begin, 20) << "\n");
	intnum = AS_MMA_TOKEN(token, Integer);
	return pstr;
}

inline const char * parse_MMA_token_as_float(const char * begin, const char * end, MMA_TOKEN & token, El::BigFloat & f)
{
	const char * pstr = parse_get_token(begin, end, token);

	if (token.index() == MMA_TOKEN_Operator && AS_MMA_TOKEN(token, Operator) == '-')
	{
		const char *pstr2 = parse_MMA_token_as_float(pstr, end, token, f);
		f = -f;
		return pstr2;
	}

	if (token.index() != MMA_TOKEN_Real && token.index() != MMA_TOKEN_Integer)
		MMA_PARSER_ERROR("parse_MMA_get_token_as_float error : expect Real, but I got " << token << " from text " << std::string(begin, 20) << "\n");

	if (token.index() == MMA_TOKEN_Integer)
		f = AS_MMA_TOKEN(token, Integer);
	else
		f = AS_MMA_TOKEN(token, Real);
	return pstr;
}

inline const char * parse_MMA_token_as_string(const char * begin, const char * end, MMA_TOKEN & token, std::string & str)
{
	const char * pstr = parse_get_token(begin, end, token);
	if (token.index() != MMA_TOKEN_String)
		MMA_PARSER_ERROR("parse_MMA_get_token_as_float error : expect String, but I got " << token << " from text " << std::string(begin, 20) << "\n");
	str = std::move(AS_MMA_TOKEN(token, String));
	return pstr;
}

////////////////// build-in function and symbols ////////////////////////



namespace param
{
	El::BigFloat dim, nu;
	int kappa;
	int maxderivs; // only used for interval positivity
	std::string block_folder;
	std::vector<std::string> input_files;
	std::map<std::string, El::BigFloat> var_map;

	Boost_Float r_crossing_4;
	bool MPI_F_FS_parallelQ = false;
}

bool internal_print_Q = false;

std::map<std::pair<std::string, int>, std::vector<std::vector<std::vector<El::BigFloat>>> > blockF;
std::map<std::pair<std::string, int>, int > blockF_key2index;

int MPI_stamp_spin_to_rank(const std::string & stamp, int L)
{
	if (blockF_key2index.find(std::make_pair(stamp, L)) == blockF_key2index.end())
		MMA_PARSER_ERROR("MPI_stamp_spin_to_rank error : can not find stamp=" << stamp << ", L=" << L << "\n");

	int index = blockF_key2index[std::make_pair(stamp, L)];

	return index % El::mpi::Size();
}


int sb_parse_vector_mpi_allreduce(std::vector<El::BigFloat> & vec) // add vec from all MPI process
{
	for (auto & v : vec)
		v = El::mpi::AllReduce(v, El::mpi::COMM_WORLD);

	return 1;
}

void generate_blockF_key2index(const std::string & block_folder, std::map<std::pair<std::string, int>, int > & blockF_key2index)
{
	namespace fs = boost::filesystem;

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

			blockF_key2index.emplace(std::make_pair(stamp, spin), 0);
		}
	}

	int i = 0;
	for (auto ptr : blockF_key2index)
		blockF_key2index[ptr.first] = i++;

	/* // check what file belong to which rank
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

			if (MPI_stamp_spin_to_rank(stamp, spin) == El::mpi::Rank())
				std::cout << "I found " << file.path().filename().string() << " belongs to current rank=" << El::mpi::Rank() << "\n";
		}
	}
	*/
}



//////////////////// cache for binomial coefficient, interval transformation ////////////////////////////////


std::vector<std::vector<mpz_class>> binomial_cache;
int binomial_cache_N = -1;

// generate table of Binomial[m,n] with m,n<=N
int init_binomial_coeff(int N)
{
	if (N < binomial_cache_N)return 0;
	binomial_cache_N = N;
	binomial_cache.resize(N + 1);
	for (int m = 0; m <= N; m++)
	{
		binomial_cache[m].resize(m + 1);

		mpz_class binomial_coeff = 1;
		for (int n = 0; n <= m; n++)
		{
			// store Binomial[m,n]
			binomial_cache[m][n] = binomial_coeff;

			// Binomial[m,n+1]=Binomial[m,n] * (m-n)/(1+n)
			binomial_coeff = (binomial_coeff * (m - n)) / (n + 1);
		}
	}
	return 1;
}

// return Binomial[m,n] , assuming m>=n
inline mpz_class binomial_coeff(int m, int n)
{
	if (m > binomial_cache_N)init_binomial_coeff(m);
	return binomial_cache[m][n];
}

// return Binomial[m,n] , assuming m>=n
inline mpz_class binomial_coeff_cached(int m, int n)
{
	return binomial_cache[m][n];
}

// ----------- interval transformation  ---------------

void interval_transformation(std::vector<El::BigFloat> & coeff, const El::BigFloat & a, const El::BigFloat & b, int max_degree)
{

	int N = coeff.size()-1; // degree of the polynomial
	int M = max_degree; // maximum degree of the polynomials in current matrix

	std::vector<El::BigFloat> new_coefficients(M + 1, 0);

	init_binomial_coeff(M);

	for (int n = 0; n <= N; n++)
	{
		El::BigFloat a_pow = as_BigFloat(pow(as_BoostFloat(a), n));
		El::BigFloat b_pow = 1;
		for (int m = 0; m <= n; m++)
		{
			for (int k = 0; k <= M - n; k++)
			{
				new_coefficients[k + m].gmp_float += coeff[n].gmp_float * a_pow.gmp_float * b_pow.gmp_float * binomial_coeff_cached(n, m) * binomial_coeff_cached(M - n, k);
			}
			a_pow = a_pow / a;
			b_pow = b_pow * b;
		}
	}

	coeff.assign(new_coefficients.begin(), new_coefficients.end());
}


////////////////////////////////////////////////////////////////////////////////////



El::BigFloat Fprefactor(int L, const El::BigFloat & x)
{
	using namespace param;
	auto & order = kappa;
	El::BigFloat denominator = 1;
	El::BigFloat Delta = x + (L + dim - 2);

	for (int64_t k = 1; k <= order; ++k)
		denominator *= Delta - (-k - L + 1);
	for (int64_t k = 1; 2 * k <= order; ++k)
		denominator *= Delta -(nu + 1 - k);
	for (int64_t k = 1; k <= std::min(order, L); ++k)
		denominator *= Delta - (1 + 2 * nu + L - k);

	El::BigFloat numerator_BigFloat(as_BigFloat(
		pow(r_crossing_4, as_BoostFloat(Delta))
	));

	//std::cout << "Fprefactor[" << L << "," << x << "]=" << std::setprecision(200) << numerator_BigFloat / denominator << "\n";

	return numerator_BigFloat / denominator;
}

El::BigFloat FSprefactor(int L, const El::BigFloat & x)
{
	using namespace param;
	auto & order = kappa;
	El::BigFloat denominator = 1;
	El::BigFloat Delta = x + (L + dim - 2);
	
	for (int64_t k = 2; k <= order; k = k + 2)
		denominator *= Delta - (-k - L + 1);
	for (int64_t k = 1; 2 * k <= order; ++k)
		denominator *= Delta - (nu + 1 - k);
	for (int64_t k = 2; k <= std::min(order, L); k = k + 2)
		denominator *= Delta - (1 + 2 * nu + L - k);

	El::BigFloat numerator_BigFloat(as_BigFloat(
		pow(r_crossing_4, as_BoostFloat(Delta))
	));

	//std::cout << "FSprefactor[" << L << "," << x << "]=" << std::setprecision(200) << numerator_BigFloat / denominator << "\n";
	//std::cout << "numerator=" << numerator_BigFloat << ", denominator=" << denominator << "\n";

	return numerator_BigFloat / denominator;
}


void load_block_folder(const std::string & block_folder, const std::string & stamp, int spin,
	std::map<std::pair<std::string, int>, std::vector<std::vector<std::vector<El::BigFloat>>> > & blockF);

inline auto & blockF_lookup(const std::string & stamp, int L, int m, int n)
{
	auto pblock = blockF.find(std::make_pair(stamp, L));
	if (pblock == blockF.end())
	{
		load_block_folder(param::block_folder, stamp, L, blockF);
		pblock = blockF.find(std::make_pair(stamp, L));
	}
	
	if (pblock->second.size() < m + 1 || pblock->second.at(m).size() < n + 1)
		MMA_PARSER_ERROR("can't find polynomial for stamp=" << stamp << " L=" << L << " m=" << m << " n=" << n << "\n");

	return pblock->second.at(m).at(n);
}


inline void simpleboot_internal_F(const std::string & stamp, int L, int m, int n, El::BigFloat & x, MMA_ELEMENT & result)
{
	if (param::MPI_F_FS_parallelQ == true && MPI_stamp_spin_to_rank(stamp, L) != El::mpi::Rank())
	{
		result = El::BigFloat(0);
		return;
	}

	//std::cout << "rank=" << El::mpi::Rank() << " : stamp=" << stamp << ", L=" << L << " belong to this rank." << " MPI_F_FS_parallelQ=" << param::MPI_F_FS_parallelQ << "\n";

	Polynomial poly;
	poly.coefficients = blockF_lookup(stamp, L, m, n);

	result = poly(x) * Fprefactor(L, x);
	//if (internal_print_Q) std::cout << "stamp=" << stamp << ", L=" << L << ", m=" << m << ", n=" << n << ", x=" << x << "\nresult=" << result << "\n";

	//std::cout << "F[" << stamp << "," << L << "," << m << "," << n << "," << x << "]=" << result <<
	//	"   poly(x)=" << poly(x) << "   Fprefactor(L, x)=" << Fprefactor(L, x) << "\n";
}

inline void simpleboot_internal_FS(const std::string & stamp, int L, int m, int n, El::BigFloat & x, MMA_ELEMENT & result)
{
	if (param::MPI_F_FS_parallelQ == true && MPI_stamp_spin_to_rank(stamp, L) != El::mpi::Rank())
	{
		result = El::BigFloat(0);
		return;
	}

	//std::cout << "rank=" << El::mpi::Rank() << " : stamp=" << stamp << ", L=" << L << " belong to this rank." << " MPI_F_FS_parallelQ=" << param::MPI_F_FS_parallelQ << "\n";

	Polynomial poly;
	poly.coefficients = blockF_lookup(stamp, L, m, n);

	result = poly(x) * FSprefactor(L, x);

	//std::cout << "FS[" << stamp << "," << L << "," << m << "," << n << "," << x << "]=" << result << 
	//	"   poly(x)=" << poly(x) << "   FSprefactor(" << L << "," << x << ")=" << FSprefactor(L, x) << "\n";
}


int current_matrix_max_polynomial_degree = -1;
void simpleboot_internal_PT(const std::string & stamp, int L, int m, int n, El::BigFloat & a, El::BigFloat & b, MMA_ELEMENT & result)
{
	Polynomial poly;

	poly.coefficients = blockF_lookup(stamp, L, m, n);

	if (current_matrix_max_polynomial_degree == -1)current_matrix_max_polynomial_degree = poly.degree() + (param::maxderivs - m - n);

	if(current_matrix_max_polynomial_degree!= poly.degree() + (param::maxderivs - m - n))
		MMA_PARSER_ERROR("inconsistent polynomial degree : prediction from stamp="<< stamp << ", L=" << L << ", m=" << m << ", n=" << n << " is " << poly.degree() + (param::maxderivs - m - n) << ", while previous prediction is " << param::maxderivs << "\n");

	if (param::maxderivs < m + n) MMA_PARSER_ERROR("incorrect maxderivs in the param file: maxderivs=" << param::maxderivs << "\n");

	interval_transformation(poly.coefficients, a, b, current_matrix_max_polynomial_degree);

	result = std::move(poly);
	//std::cout << "result=" << result << "\n";
	return;
}



void simpleboot_internal_P(const std::string & stamp, int L, int m, int n, El::BigFloat & shift, MMA_ELEMENT & result)
{
	Polynomial poly;
//	result = poly;
//	return;

	poly.coefficients = blockF_lookup(stamp, L, m, n);

	if (shift != El::BigFloat(0))
	{
		poly.shift(shift);
		//std::cout << "shift=" << shift << ", prec=" << shift.Precision() << "\n";
	}
	
	result = std::move(poly);
	//std::cout << "result=" << result << "\n";
	return;
}



// for my purpose now, I only need
// P[stamp, L, m, n, shift] : block derivative polynomial
// F[stamp, L, m, n, x, dim, kappa] : block derivative value with prefactor
// Fs[stamp, L, m, n, x, dim, kappa] : block derivative value with prefactor (shortened pole)
const char * parse_MMA_function(const std::string & name, const char * begin, const char * end, MMA_ELEMENT & result)
{
	const char * pstr = begin;
	if (name.compare("F") == 0 || name.compare("FS") == 0)
	{
		MMA_TOKEN token;
		MMA_ELEMENT element;

		int L, m, n, kappa;
		El::BigFloat dim, x;

		pstr = parse_get_token(pstr, end, token); // F[stamp, L, m, n, x, dim, kappa]
		std::string stamp(std::move(AS_MMA_TOKEN(token, String)));
		pstr = parse_MMA_check_op(pstr, end, token, ',');

		pstr = parse_MMA_token_as_int(pstr, end, token, L);
		pstr = parse_MMA_check_op(pstr, end, token, ',');

		pstr = parse_MMA_token_as_int(pstr, end, token, m);
		pstr = parse_MMA_check_op(pstr, end, token, ',');

		pstr = parse_MMA_token_as_int(pstr, end, token, n);
		pstr = parse_MMA_check_op(pstr, end, token, ',');

		pstr = parse_MMA_expr_as_number(pstr, end, x);
		pstr = parse_MMA_check_op(pstr, end, token, ']');

		if(name.compare("F") == 0)
			simpleboot_internal_F(stamp, L, m, n, x, result);
		else
			simpleboot_internal_FS(stamp, L, m, n, x, result);

		return pstr;
	}

	if (name.compare("P") == 0)
	{
		MMA_TOKEN token;
		MMA_ELEMENT element;

		int L, m, n;
		El::BigFloat x;

		pstr = parse_get_token(pstr, end, token); // P[stamp, L, m, n, shift]
		std::string stamp(std::move(AS_MMA_TOKEN(token, String)));
		pstr = parse_MMA_check_op(pstr, end, token, ',');

		pstr = parse_MMA_token_as_int(pstr, end, token, L);
		pstr = parse_MMA_check_op(pstr, end, token, ',');

		pstr = parse_MMA_token_as_int(pstr, end, token, m);
		pstr = parse_MMA_check_op(pstr, end, token, ',');

		pstr = parse_MMA_token_as_int(pstr, end, token, n);
		pstr = parse_MMA_check_op(pstr, end, token, ',');

		pstr = parse_MMA_expr_as_number(pstr, end, x);
		pstr = parse_MMA_check_op(pstr, end, token, ']');

		simpleboot_internal_P(stamp, L, m, n, x, result);

		return pstr;
	}

	if (name.compare("PT") == 0)
	{
		MMA_TOKEN token;
		MMA_ELEMENT element;

		int L, m, n;
		El::BigFloat a, b;

		pstr = parse_get_token(pstr, end, token); // PT[stamp, L, m, n, a, b]
		std::string stamp(std::move(AS_MMA_TOKEN(token, String)));
		pstr = parse_MMA_check_op(pstr, end, token, ',');

		pstr = parse_MMA_token_as_int(pstr, end, token, L);
		pstr = parse_MMA_check_op(pstr, end, token, ',');

		pstr = parse_MMA_token_as_int(pstr, end, token, m);
		pstr = parse_MMA_check_op(pstr, end, token, ',');

		pstr = parse_MMA_token_as_int(pstr, end, token, n);
		pstr = parse_MMA_check_op(pstr, end, token, ',');

		pstr = parse_MMA_expr_as_number(pstr, end, a);
		pstr = parse_MMA_check_op(pstr, end, token, ',');

		pstr = parse_MMA_expr_as_number(pstr, end, b);
		pstr = parse_MMA_check_op(pstr, end, token, ']');

		simpleboot_internal_PT(stamp, L, m, n, a, b, result);

		return pstr;
	}


	MMA_PARSER_ERROR("parse_MMA_function error : unsupport " << name << " \n");
}


void parse_MMA_symbol(const std::string & name, MMA_ELEMENT & result)
{
	auto pvar = param::var_map.find(name);
	if (pvar == param::var_map.end())
		MMA_PARSER_ERROR("can't find symbol " << name << "\n");

	SET_MMA_ELEMENT(result, pvar->second, Expression);

	//std::cout << "parse_MMA_symbol find symbol " << name << " = " << result << "\n";
	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


MPI_Win the_window;
int *window_data;
int counter_process = 0;

void mpi_counter_init(int init_value = 0)
{
	MPI_Win_allocate(sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &window_data, &the_window);

	int counter_init = init_value;

	MPI_Win_fence(0, the_window);
	if (El::mpi::Rank() == counter_process)
		MPI_Put(&counter_init, 1, MPI_INT, counter_process, 0, 1, MPI_INT, the_window);
	MPI_Win_fence(0, the_window);
}

int mpi_counter_get()
{
	//std::cout << "[Rank=" << El::mpi::Rank() << " mpi_counter_get BEGIN]\n";
	//MPI_Win_fence(0, the_window);
	int counter_value;
	int decrement = 1;
	MPI_Fetch_and_op(&decrement, &counter_value, MPI_INT, counter_process, 0, MPI_SUM, the_window);
	//MPI_Win_fence(0, the_window);
	//std::cout << "[Rank=" << El::mpi::Rank() << " mpi_counter_get   END] " << counter_value << "\n";

	return counter_value;
}

void test_mpi()
{
	mpi_counter_init();

	std::vector<int> matrix_indices;

	int counter;
	do
	{
		counter = mpi_counter_get();
		matrix_indices.push_back(counter);
		std::cout << "current rank : " << El::mpi::Rank() << " get jobid=" << counter << "\n";
	} while (counter < 200);

	std::cout << "[Rank=" << El::mpi::Rank() << "] " << "matrices_valid_indices=";
	for (auto i : matrix_indices) std::cout << i << " ";
	std::cout << "\n";
	
	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void test_parse_MMA_expr();

void test_gmp_mpfr();

void test_PT();

template <typename T>
const char *
sb_parse_vector(const char *begin, const char *end, std::vector<T> &result_vector)
{
	//test_PT();
	//exit(0);

	//std::cout << std::setprecision(50) << std::fixed;
	//test_gmp_mpfr();
	//exit(0);

	const auto open_brace(std::find(begin, end, '{'));
	if (open_brace == end)
	{
		throw std::runtime_error("Missing '{' at beginning of array of numbers");
	}
	auto start_element(std::next(open_brace));

	const auto close_brace(std::find(start_element, end, '}'));
	if (close_brace == end)
	{
		throw std::runtime_error("Missing '}' at end of array of numbers");
	}
	const auto brace_end = std::next(close_brace);

	MMA_ELEMENT element;
	MMA_TOKEN token;
	while (start_element <= close_brace)
	{
		start_element = parse_MMA_expr(start_element, brace_end, element);
		start_element = parse_get_token(start_element, brace_end, token);

		if (element.index() != MMA_ELEMENT_Expression || AS_MMA_ELEMENT(element, Expression).index() != MMA_EXPR_Number)
			MMA_PARSER_ERROR("sb_parse_vector error : expect number , but I got " << element << " before " << std::string(start_element, 20) << "\n");

		if constexpr (std::is_same_v<T, El::BigFloat>)
		{
			result_vector.emplace_back(AS_MMA_EXPR(AS_MMA_ELEMENT(element, Expression), Number));
		}
		else if constexpr (std::is_same_v<T, Boost_Float>)
		{
			result_vector.emplace_back(boost::multiprecision::mpf_float(
				AS_MMA_EXPR(AS_MMA_ELEMENT(element, Expression), Number).gmp_float.get_mpf_t()
			));
		}

		//std::cout << "[" << result_vector.size() << "]=" << result_vector.back() << "\n";

		if (token.index() != MMA_TOKEN_Operator)
			MMA_PARSER_ERROR("sb_parse_vector error : expect ',' or '}' , but I got " << token << " before " << std::string(start_element, 20) << "\n");

		if (AS_MMA_TOKEN(token, Operator) == '}') break;

		if (AS_MMA_TOKEN(token, Operator) != ',')
			MMA_PARSER_ERROR("sb_parse_vector error : expect ',' , but I got " << token << " before " << std::string(start_element, 20) << "\n");
	}

	
	/*
	int prev_prec = std::cout.precision();
	//std::cout << "from text : \n";
	//std::cout << std::string(begin, std::next(close_brace));
	std::cout << "\nRank=" << El::mpi::Rank() << "\nparse vector : {" << std::setprecision(20);
	for (auto & a : result_vector) std::cout << a << ", ";
	std::cout << "}\n" << std::setprecision(prev_prec);	
	exit(0);
	*/

	return std::next(close_brace);
}

template const char *sb_parse_vector(const char *begin, const char *end, std::vector<El::BigFloat> &result_vector);
template const char *sb_parse_vector(const char *begin, const char *end, std::vector<Boost_Float> &result_vector);



const char * sb_parse_polynomial(const char *begin, const char *end, Polynomial &polynomial)
{
	internal_print_Q = true;

	MMA_ELEMENT element;
	const char * pstr = parse_MMA_expr(begin, end, element);

	if (element.index() != MMA_ELEMENT_Expression)
		MMA_PARSER_ERROR("sb_parse_vector error : expect polynomial , but I got " << element << " from " << std::string(begin, 20) << "\n");

	MMA_EXPR & expr = AS_MMA_ELEMENT(element, Expression);  // is this copied

	if (expr.index() == MMA_EXPR_Polynomial)
	{
		polynomial = std::move(AS_MMA_EXPR(AS_MMA_ELEMENT(element, Expression), Polynomial));
	}
	else
	{
		polynomial.coefficients.clear();
		polynomial.coefficients.emplace_back(AS_MMA_EXPR(expr, Number));
	}

	//int prev_prec = std::cout.precision();
	//std::cout << "parse polynomial : from text = " << std::string(begin, pstr) << "\n" << std::setprecision(200) << polynomial << std::setprecision(prev_prec) << "\n";
	//std::cout << "parse polynomial : " << std::setprecision(200) << polynomial << std::setprecision(prev_prec) << "\n";

	return pstr;
}


/////////////////////////////// parse damped rational constant //////////////////////////////////////////////////////////


const char * sb_parse_damped_rational_constant(const char * constant_start, const char * end, Boost_Float & damped_rational_constant)
{
	El::BigFloat damped_rational_constant_BigFloat;
	const char * pstr = parse_MMA_expr_as_number(constant_start, end, damped_rational_constant_BigFloat);
	if (pstr == end) throw std::runtime_error("Missing comma after DampedRational.constant");

	damped_rational_constant = as_BoostFloat(damped_rational_constant_BigFloat);

	return pstr;
}


/////////////////////////////// parse parameter file //////////////////////////////////////////////////////////



int parse_parameter_find_item(const char * begin, const char * end, const std::string & itemname, const char * & begin_item, const char * & end_item, bool required = true)
{
	const std::string item_head("<" + itemname + ">"), item_tail("</" + itemname + ">");

	begin_item = std::search(begin, end, item_head.begin(), item_head.end());
	end_item = std::search(begin, end, item_tail.begin(), item_tail.end());

	if (begin_item == end && required == false) return 0;
	begin_item = begin_item + item_head.size();

	if (begin <= begin_item && begin_item < end_item && end_item <= end) return 1;
	
	MMA_PARSER_ERROR("parse_parameter_file error : can't process item " << itemname << " correctly\n");
	return 0;
}

void load_block_folder(const std::string & block_folder, std::map<std::pair<std::string, int>, std::vector<std::vector<std::vector<El::BigFloat>>> > & blockF);

const char * parse_parameter_file(const char * begin, const char * end)
{
	using namespace param;

	const char * begin_item;
	const char * end_item;
	MMA_TOKEN token;

	parse_parameter_find_item(begin, end, "kappa", begin_item, end_item);
	parse_MMA_token_as_int(begin_item, end_item, token, kappa);

	// this is only used for interval positivity
	if (parse_parameter_find_item(begin, end, "maxderivs", begin_item, end_item, false))
		parse_MMA_token_as_int(begin_item, end_item, token, maxderivs);
	else
		maxderivs = 0;

	El::BigFloat dim_temp; // If I directly using dim, the dim.Precision() is not correct
	parse_parameter_find_item(begin, end, "dim", begin_item, end_item);
	parse_MMA_token_as_float(begin_item, end_item, token, dim_temp);

	dim = std::move(dim_temp); // somehow without std::move, the precision is not correct
	nu = (dim - 2) / 2;
	Boost_Float r_crossing_4_temp = (3 - 2 * sqrt(Boost_Float(2))) * 4;
	r_crossing_4 = r_crossing_4_temp; // the precision is correct without std::move

	//std::cout << "parse_parameter_file : \n";
	//std::cout << "dim_temp=" << dim_temp << ", prec=" << dim_temp.Precision() << "\n";
	//std::cout << "dim=" << dim << ", prec=" << dim.Precision() << "\n";
	//std::cout << "r_crossing_4_temp=" << r_crossing_4_temp << ", prec=" << mpfr_get_prec(r_crossing_4_temp.backend().data()) << "\n";
	//std::cout << "r_crossing_4=" << r_crossing_4 << ", prec=" << mpfr_get_prec(r_crossing_4.backend().data()) << "\n";


	parse_parameter_find_item(begin, end, "block", begin_item, end_item);
	parse_MMA_token_as_string(begin_item, end_item, token, block_folder);

	parse_parameter_find_item(begin, end, "input", begin_item, end_item);
	const char * pstr = parse_MMA_check_op(begin_item, end_item, token, '{');
	char op;
	std::string filename;
	while (1)
	{
		pstr = parse_MMA_token_as_string(pstr, end_item, token, filename);

		std::cout << "find input files : " << filename << "\n";

		input_files.push_back(std::move(filename));
		pstr = parse_MMA_get_op(pstr, end_item, token, op);
		if (op == '}') break;
		if (op != ',') MMA_PARSER_ERROR("parse_parameter_file error : expecting ',' , but I get " << op << " before " << std::string(pstr, 20) << "\n");
	}

	if (parse_parameter_find_item(begin, end, "variables", begin_item, end_item, false) != 0)
	{
		pstr = parse_MMA_check_op(begin_item, end_item, token, '{');
		std::string var_name;
		El::BigFloat var_value;

		while (1)
		{
			pstr = parse_MMA_token_as_string(pstr, end_item, token, var_name);
			pstr = parse_MMA_check_op(pstr, end_item, token, ',');
			pstr = parse_MMA_token_as_float(pstr, end_item, token, var_value);

			var_map.emplace(var_name, var_value);

			pstr = parse_MMA_get_op(pstr, end_item, token, op);
			if (op == '}') break;
			if (op != ',') MMA_PARSER_ERROR("parse_parameter_file error : expecting ',' , but I get " << op << " before " << std::string(pstr, 20) << "\n");
		}
	}

	//std::cout << std::setprecision(50) << std::fixed;

	std::cout << "parameter file processed\n";
	std::cout << "dim=" << dim << "\n";
	std::cout << "kappa=" << kappa << "\n";
	std::cout << "block=" << block_folder << "\n";

	std::cout << "input={";
	for (auto & file : input_files) std::cout << file << " ";
	std::cout << "}\n";

	std::cout << "variables={\n";
	for (const auto&[key, value] : var_map) std::cout << key << " = " << value << "\n";
	std::cout << "}\n";

	//load_block_folder(param::block_folder, blockF);
	generate_blockF_key2index(param::block_folder, blockF_key2index);

	return end;
}

int parse_parameter_file(boost::filesystem::path & param_file)
{
	boost::filesystem::ifstream input_stream(param_file);
	if (!input_stream.good())
	{
		throw std::runtime_error("Unable to open parameter file: " + param_file.string());
	}

	boost::interprocess::file_mapping mapped_file(
		param_file.c_str(), boost::interprocess::read_only);
	boost::interprocess::mapped_region mapped_region(
		mapped_file, boost::interprocess::read_only);

	try
	{
		const char *begin(static_cast<const char *>(mapped_region.get_address())),
			*end(begin + mapped_region.get_size());
		parse_parameter_file(begin, end);
	}
	catch (std::exception &e)
	{
		throw std::runtime_error("Error when parsing parameter file " + param_file.string()
			+ ": " + e.what());
	}
	return 1;
}



/////////////////////////////// start parse input //////////////////////////////////////////////////////////




void read_input(const boost::filesystem::path &input_file,
	std::vector<El::BigFloat> &objectives,
	std::vector<El::BigFloat> &normalization,
	std::vector<Positive_Matrix_With_Prefactor> &matrices,
	std::vector<int> &matrices_valid_indices);

void sb_parse_input(std::vector<El::BigFloat> &objectives,
	std::vector<El::BigFloat> &normalization,
	std::vector<Positive_Matrix_With_Prefactor> &matrices,
	std::vector<int> &matrices_valid_indices
	)
{
	for (auto &filename : param::input_files)
		read_input(filename, objectives, normalization, matrices, matrices_valid_indices);
	return;
}


/////////////////////////////// test code //////////////////////////////////////////////////////////


void test_PT()
{
	Polynomial poly;
	poly.coefficients = { 2,3.1,3.2,-10.32 };

	std::cout << "test interval transformation : p(x) = " << poly << "\n";

	interval_transformation(poly.coefficients, 2.3, 5.3, 20);

	std::cout << "transformed p = " << poly << "\n";
}


void test_gmp_mpfr()
{
	using namespace param;

	std::cout << "------- test gmp begin ----------------\n";

	// set constant
	//r_crossing_4 = (3 - 2 * sqrt(Boost_Float(2))) * 4;

	Boost_Float rstar4 = (3 - 2 * sqrt(Boost_Float(2))) * 4;

	std::cout << "mpfr_get_default_prec() =" << mpfr_get_default_prec() << "\n";
	std::cout << "Boost_Float default precision=" << Boost_Float::default_precision() << "\n";

	std::cout << std::setprecision(200) << std::fixed;
	std::cout << "rstar4=" << rstar4 << ", prec=" << mpfr_get_prec(rstar4.backend().data()) << "\n";

	std::cout << "r_crossing_4=" << r_crossing_4 << ", prec=" << mpfr_get_prec(r_crossing_4.backend().data()) << "\n";

	r_crossing_4 = rstar4;

	std::cout << "r_crossing_4=" << r_crossing_4 << ", prec=" << mpfr_get_prec(r_crossing_4.backend().data()) << "\n";

	std::cout << "dim=" << dim << ", prec=" << dim.Precision() << "\n";
	std::cout << "nu=" << nu << ", prec=" << nu.Precision() << "\n";

	FSprefactor(2, "3.14");
	FSprefactor(1, "3.145");

	Fprefactor(2, "3.14");
	Fprefactor(1, "3.145");

	std::cout << "FSprefactor(0, 0.412)=" << FSprefactor(0, "0.412") << "\n";
	std::cout << "FSprefactor(2, 0)=" << FSprefactor(2, "0") << "\n";

	std::cout << "--------- test gmp end ----------------\n";

	exit(0);
	return;
}




//std::string test_str = "2.3 - 0.3 P[stamp, 1, 1, 2, -0.134] + P[stamp, 1, 1, 2, 3.1*sym+13.2`200]/(1+3) - P[stamp, 0, 1, 2, 0.134] +1.3 , 2+3  ";

std::string test_str = "F[\"Fespsigespsig\", 0, 1, 0, -0.4819999999999999840127884453977458178997\\\n\
03979492187499999999999999999999999999999999999999999999999999999999999999999\\\n\
99999999999999999999999999999999999999999999999999999999999999999999999999999\\\n\
999999999999999999`200.]";


void test_parse_MMA_expr()
{
	std::cout << "test=" << test_str << "\n";

	const char * begin = test_str.c_str();
	const char * end = test_str.c_str() + test_str.size();

	MMA_ELEMENT result;
	parse_MMA_expr(begin, end, result);

	std::cout << "result=" << result << "\n";

	return;
}

void test_MMA_element()
{
	const char * begin = test_str.c_str();
	const char * end = test_str.c_str() + test_str.size();

	std::cout.precision(15);

	const char * pstr = begin;
	int i = 1;
	while (pstr != end && i <= 5000)
	{
		MMA_ELEMENT element;

		std::cout << "current str = " << pstr << "\n";

		pstr = parse_MMA_element(pstr, end, element);
		std::cout << "Find element #" << i << " : " << element << "\n";

		if (element.index() == 0)break;

		i++;
	}
	exit(0);
	return;
}



/////////////////////// trash //////////////////////////////////////

/*
inline void simpleboot_internal_P2(const std::string & stamp, int L, int m, int n, El::BigFloat & shift, MMA_ELEMENT & result)
{
Polynomial poly;
switch (L)
{
case 0:
poly.coefficients.assign({ m,2,4,7.5 });
poly.shift(shift);
break;
case 1:
poly.coefficients.assign({ m + n, 1.7 });
poly.shift(shift);
break;
default:
poly.coefficients.assign({ 1, 0 });
break;
}
result = poly;
return;
}
*/

