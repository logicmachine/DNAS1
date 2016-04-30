#pragma once

#include <vector>
#include <string>
#include <algorithm>
#include <cstdint>

class Query {

private:
	std::array<std::vector<uint8_t>, 2> m_fragments;
	std::array<std::vector<uint8_t>, 2> m_rev_fragments;

	std::vector<uint8_t> translate(const std::string& s) const {
		const int n = s.size();
		std::vector<uint8_t> t(n);
		for(int i = 0; i < n; ++i){
			switch(s[i]){
				case 'A': t[i] = 0; break;
				case 'T': t[i] = 1; break;
				case 'C': t[i] = 2; break;
				case 'G': t[i] = 3; break;
			}
		}
		return t;
	}

	std::vector<uint8_t> reverse_complement(
		const std::vector<uint8_t>& s) const
	{
		std::vector<uint8_t> t(s);
		std::reverse(t.begin(), t.end());
		for(uint8_t& c : t){ c ^= 1; }
		return t;
	}

public:
	Query(const std::string& s1, const std::string& s2)
		: m_fragments{ translate(s1), translate(s2) }
		, m_rev_fragments{
			reverse_complement(m_fragments[0]),
			reverse_complement(m_fragments[1])}
	{ }

	const std::vector<uint8_t>& fragments(int i) const {
		return m_fragments[i];
	}
	const std::vector<uint8_t>& rev_fragments(int i) const {
		return m_rev_fragments[i];
	}

};

