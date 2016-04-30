#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <utility>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include "solution.hpp"

class Problem {

private:
	using IndexedChromatid = std::pair<int, std::vector<std::string>>;

	double m_norm_a;
	double m_norm_s;
	unsigned int m_preprocess_time_limit;
	unsigned int m_time_cutoff;

	boost::filesystem::path m_minisam_path;
	std::array<boost::filesystem::path, 2> m_fa_paths;
	std::vector<std::pair<int, boost::filesystem::path>> m_chromatid_paths;

	std::vector<std::string> m_read_names;
	std::vector<std::string> m_read_sequences;
	std::vector<IndexedChromatid> m_chromatids;

	std::vector<Solution> m_expected_solutions;

	boost::filesystem::path read_path(
		std::istream& is, const boost::filesystem::path& base) const
	{
		std::string line;
		while(line.empty()){
			std::getline(is, line);
			while(isspace(line.back())){ line.pop_back(); }
		}
		return boost::filesystem::absolute(
			boost::filesystem::path(line), base);
	}

	void read_queries(){
#ifdef VERBOSE
		std::cerr << "Load: " << m_fa_paths[0].c_str() << std::endl;
		std::cerr << "Load: " << m_fa_paths[1].c_str() << std::endl;
#endif
		std::ifstream ifs1(m_fa_paths[0].c_str());
		std::ifstream ifs2(m_fa_paths[1].c_str());
		std::string name1, name2, seq1, seq2;
		while(std::getline(ifs1, name1) && std::getline(ifs2, name2)){
			while(isspace(name1.back())){ name1.pop_back(); }
			while(isspace(name2.back())){ name2.pop_back(); }
			std::getline(ifs1, seq1);
			std::getline(ifs2, seq2);
			while(isspace(seq1.back())){ seq1.pop_back(); }
			while(isspace(seq2.back())){ seq2.pop_back(); }
			m_read_names.push_back(name1.substr(1));
			m_read_names.push_back(name2.substr(1));
			m_read_sequences.push_back(seq1);
			m_read_sequences.push_back(seq2);
		}
	}

	std::vector<std::string>
	read_chromatid(const boost::filesystem::path& path) const {
#ifdef VERBOSE
		std::cerr << "Load: " << path.c_str() << std::endl;
#endif
		std::ifstream ifs(path.c_str());
		std::string line;
		std::getline(ifs, line);
#ifdef VERBOSE
		std::cerr << line << std::endl;
#endif
		std::vector<std::string> rows;
		while(std::getline(ifs, line)){
			if(line.empty()){ continue; }
			while(isspace(line.back())){ line.pop_back(); }
			rows.emplace_back(std::move(line));
		}
		return rows;
	}

	size_t chromatid_length(const std::vector<std::string> &vs) const {
		size_t sum = 0;
		for(const auto &s : vs){ sum += s.size(); }
		return sum;
	}

	std::vector<Solution> read_minisam(
		const boost::filesystem::path& path) const
	{
#ifdef VERBOSE
		std::cerr << "Load: " << path.c_str() << std::endl;
#endif
		std::ifstream ifs(path.c_str());
		std::string line1, line2;
		std::vector<Solution> result;
		while(std::getline(ifs, line1)){
			for(char& c : line1){
				if(c == ','){ c = ' '; }
			}
			std::getline(ifs, line2);
			for(char& c : line2){
				if(c == ','){ c = ' '; }
			}
			std::istringstream iss1(line1), iss2(line2);
			std::array<std::istringstream *, 2> iss = { &iss1, &iss2 };
			std::string name[2], strand[2];
			int chromatid_id[2], start[2], end[2];
			for(int i = 0; i < 2; ++i){
				*(iss[i]) >> name[i] >> chromatid_id[i];
				*(iss[i]) >> start[i] >> end[i] >> strand[i];
			}
			result.emplace_back(
				chromatid_id[0], 0.0,
				start[0] - 1, end[0], strand[0][0] == '-',
				start[1] - 1, end[1], strand[1][0] == '-');
		}
		return result;
	}

public:
	Problem()
		: m_norm_a(0.0)
		, m_norm_s(0.0)
		, m_preprocess_time_limit(0)
		, m_time_cutoff(0)
		, m_minisam_path()
		, m_fa_paths()
		, m_chromatid_paths()
		, m_read_names()
		, m_read_sequences()
		, m_chromatids()
		, m_expected_solutions()
	{ }

	explicit Problem(const char *filename)
		: m_norm_a(0.0)
		, m_norm_s(0.0)
		, m_preprocess_time_limit(0)
		, m_time_cutoff(0)
		, m_minisam_path()
		, m_fa_paths()
		, m_chromatid_paths()
		, m_read_names()
		, m_read_sequences()
		, m_chromatids()
		, m_expected_solutions()
	{
		// Read meta file
#ifdef VERBOSE
		std::cerr << "Load: " << filename << std::endl;
#endif
		const boost::filesystem::path metafile_path(filename);
		const auto base_path = metafile_path.parent_path();
		std::ifstream ifs(filename);
		ifs >> m_norm_a >> m_norm_s;
		ifs >> m_preprocess_time_limit >> m_time_cutoff;
		m_minisam_path = read_path(ifs, base_path);
		m_fa_paths[0]  = read_path(ifs, base_path);
		m_fa_paths[1]  = read_path(ifs, base_path);
		std::string line;
		while(std::getline(ifs, line)){
			while(isspace(line.back())){ line.pop_back(); }
			std::istringstream iss(line);
			int index;
			iss >> index;
			while(isspace(iss.peek())){ iss.ignore(); }
			const auto path = read_path(iss, base_path);
			m_chromatid_paths.emplace_back(index, path);
		}
		ifs.close();
		// Read resources
		read_queries();
		for(const auto& p : m_chromatid_paths){
			m_chromatids.emplace_back(p.first, read_chromatid(p.second));
		}
		m_expected_solutions = read_minisam(m_minisam_path);
	}

	double norm_a() const { return m_norm_a; }
	double norm_s() const { return m_norm_s; }

	const std::vector<std::pair<int, std::vector<std::string>>>& chromatids() const {
		return m_chromatids;
	}

	const std::vector<std::string>& read_names() const {
		return m_read_names;
	}
	const std::vector<std::string>& read_sequences() const {
		return m_read_sequences;
	}

	const std::vector<Solution>& expected_solutions() const {
		return m_expected_solutions;
	}

	const unsigned int time_cutoff() const { return m_time_cutoff; }

	void print_information(std::ostream& os) const {
		os << "## Problem Description" << std::endl;
		os << "- NormA: " << m_norm_a << std::endl;
		os << "- NormS: " << m_norm_s << std::endl;
		os << "- Preprocessing time limit: "
           << m_preprocess_time_limit << " [s]" << std::endl;
		os << "- TimeCutoff: " << m_time_cutoff << " [s]" << std::endl;
		os << "- # of queries: " << m_read_names.size() << std::endl;
		os << std::endl;
		os << "### Chromatids" << std::endl;
		const size_t n = m_chromatid_paths.size();
		for(size_t i = 0; i < n; ++i){
			const auto& p = m_chromatid_paths[i];
			const auto len = chromatid_length(m_chromatids[i].second);
			os << "- " << p.first << ": " << p.second.filename()
			   << " (" << len << " bytes)" << std::endl;
		}
		os << std::endl;
	}

};

