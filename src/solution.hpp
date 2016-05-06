#pragma once

#include <iomanip>
#include <array>
#include <string>
#include <sstream>

class Solution {

public:
	struct Position {
		int start;
		int end;
		bool is_reversed;

		Position()
			: start(0)
			, end(1)
			, is_reversed(false)
		{ }

		Position(int s, int e, bool is_rev)
			: start(s)
			, end(e)
			, is_reversed(is_rev)
		{ }
	};

private:
	int m_chromatid_id;
	double m_confidence;
	std::array<Position, 2> m_positions;

public:
	Solution()
		: m_chromatid_id(20)
		, m_confidence(0.0)
		, m_positions()
	{ }

	Solution(
		int chromatid_id, double confidence,
		int start1, int end1, bool is_rev1,
		int start2, int end2, bool is_rev2)
		: m_chromatid_id(chromatid_id)
		, m_confidence(confidence)
		, m_positions({
			Position(start1, end1, is_rev1),
			Position(start2, end2, is_rev2) })
	{ }

	int chromatid_id() const noexcept { return m_chromatid_id; }
	Solution& chromatid_id(int id) noexcept {
		m_chromatid_id = id;
		return *this;
	}

	double confidence() const noexcept { return m_confidence; }
	Solution& confidence(double c) noexcept {
		m_confidence = c;
		return *this;
	}

	const Position& positions(int i) const noexcept {
		return m_positions[i];
	}
	Position& positions(int i) noexcept {
		return m_positions[i];
	}

	std::string to_string(int problem_id, int side) const {
		std::ostringstream oss;
		oss << std::setiosflags(std::ios::fixed);
		oss << std::setprecision(10);
		oss << "sim" << (problem_id + 1) << "/" << (side + 1) << ",";
		oss << m_chromatid_id << ",";
		oss << (m_positions[side].start + 1) << ",";
		oss << m_positions[side].end << ",";
		oss << (m_positions[side].is_reversed ? "-" : "+") << ",";
		oss << m_confidence;
		return oss.str();
	}

};

