#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <memory>
#include <unordered_map>
#include <cassert>
#include "query.hpp"
#include "solution.hpp"

class DNASequencing {

private:
	using Position = std::pair<int, int>;
	using ChromatidRange = std::tuple<int, int, int>;

	static const int INVALID_BASE = 255u;
	static const int BLOCK_SIZE = 32;
	static const int RANGE_PREFIX_LENGTH = 2;
	static const int RANGE_MAXSIZE = 180;
	static const int MIN_INTERVAL = 200;
	static const int MAX_INTERVAL = 500;

	std::array<std::vector<uint8_t>, 25> m_chromatids;
	std::unordered_multimap<uint64_t, Position> m_block_table;

	void build_block_table(int cid, const std::vector<uint8_t>& v){
		const std::size_t n = v.size(), m = BLOCK_SIZE;
		for(std::size_t i = 0; i + m <= n; i += m){
			bool accept = true;
			uint64_t h = 0;
			for(std::size_t j = 0; accept && j < m; ++j){
				if(v[i + j] == INVALID_BASE){ accept = false; }
				h = ((h << 2) | v[i + j]);
			}
			if(accept){
				m_block_table.emplace(h, Position(cid, i));
			}
		}
	}

	std::vector<ChromatidRange> enumerate_block_matches(
		const std::vector<uint8_t>& pattern) const
	{
		const int n = pattern.size(), m = BLOCK_SIZE;
		uint64_t h = 0;
		for(int i = 0; i + 1 < m; ++i){ h = (h << 2) | pattern[i]; }
		std::vector<Position> heads;
		for(int i = 0; i + m <= n; ++i){
			const auto j = i + m - 1;
			h = ((h << 2) | pattern[j]);
			const auto range = m_block_table.equal_range(h);
			for(auto it = range.first; it != range.second; ++it){
				const auto p = it->second;
				const int cid = p.first, offset = std::max(0, p.second - i);
				heads.emplace_back(cid, offset);
			}
		}
		std::sort(heads.begin(), heads.end());
		std::vector<ChromatidRange> result;
		for(const auto& pos : heads){
			const int chromatid_len = m_chromatids[pos.first].size();
			if(result.empty() ||
			   std::get<0>(result.back()) != pos.first ||
			   std::get<2>(result.back()) + RANGE_PREFIX_LENGTH < pos.second)
			{
				result.emplace_back(
					pos.first,
					std::max(0, pos.second - RANGE_PREFIX_LENGTH),
					std::min(chromatid_len, pos.second + RANGE_MAXSIZE));
			}else{
				std::get<2>(result.back()) =
					std::min(chromatid_len, pos.second + RANGE_MAXSIZE);
			}
		}
		return result;
	}

	std::vector<std::pair<int, int>> smith_waterman(
		const std::vector<uint8_t>& a,
		const std::vector<uint8_t>& b) const
	{
		const int INF = 1000000000;
		const int n = a.size(), m = b.size();
		std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, INF));
		for(int i = 0; i <= n; ++i){ dp[i][m] = 0; }
		for(int j = m; j > 0; --j){ dp[n][j - 1] = dp[n][j] + 1; }
		for(int i = n - 1; i >= 0; --i){
			for(int j = m -1; j >= 0; --j){
				const int x = dp[i][j + 1] + 1;
				const int y = dp[i + 1][j] + 1;
				const int xy = dp[i + 1][j + 1] + (a[i] == b[j] ? 0 : 1);
				dp[i][j] = std::min(std::min(x, y), xy);
			}
		}
		std::vector<std::pair<int, int>> result(n);
		for(int k = 0; k < n; ++k){
			int i = k, j = 0;
			while(i < n && j < m){
				const int x = dp[i][j + 1];
				const int y = dp[i + 1][j];
				const int xy = dp[i + 1][j + 1];
				const auto minval = std::min(std::min(x, y), xy);
				if(minval == xy){
					++i; ++j;
				}else if(minval == x){
					++j;
				}else if(minval == y){
					++i;
				}
			}
			result[k] = std::make_pair(dp[k][0], i);
		}
		return result;
	}

	double compute_confidence(double align, double delta) const {
		const double inv_vcm[2][2] = {
			{  1.0,                  -0.003784544765752076 },
			{ -0.003784544765752076,  1.0                  }
		};
		align = align / 7.799925841261347;
		delta = (delta - 300.1853710156314) / 133.77059238237817;
		const double t0 = align * inv_vcm[0][0] + delta * inv_vcm[1][0];
		const double t1 = align * inv_vcm[0][1] + delta * inv_vcm[1][1];
		return sqrt(t0 * align + t1 * delta);
	}

	void solve_pair(
		Solution& solution,
		std::priority_queue<double, std::vector<double>, std::greater<double>>& pq,
		const std::vector<uint8_t>& head,
		const std::vector<uint8_t>& tail,
		bool is_strand) const
	{
		const auto hranges = enumerate_block_matches(head);
		const auto tranges = enumerate_block_matches(tail);
		const int n = hranges.size(), m = tranges.size();
		for(int i = 0, j = 0; i < n; ++i){
			while(j < m){
				if(std::get<0>(hranges[i]) < std::get<0>(tranges[j])){
					break;
				}else if(std::get<0>(hranges[i]) == std::get<0>(tranges[j])){
					const auto t = std::get<1>(hranges[i]) + MIN_INTERVAL;
					if(t <= std::get<2>(tranges[j])){ break; }
				}
				++j;
			}

			const auto& chromatid = m_chromatids[std::get<0>(hranges[i])];
			const std::vector<uint8_t> head_substr(
				chromatid.begin() + std::get<1>(hranges[i]),
				chromatid.begin() + std::get<2>(hranges[i]));
			std::vector<std::pair<int, int>> head_scores;

			int k = j;
			while(k < m){
				const auto t = std::get<2>(hranges[i]) + MAX_INTERVAL;
				if(t <= std::get<1>(tranges[k])){ break; }
				if(head_scores.empty()){
					head_scores = smith_waterman(head_substr, head);
				}
				const std::vector<uint8_t> tail_substr(
					chromatid.begin() + std::get<1>(tranges[k]),
					chromatid.begin() + std::get<2>(tranges[k]));
				const auto tail_scores = smith_waterman(tail_substr, tail);
				const int hlen = head_scores.size();
				const int tlen = tail_scores.size();
				double best_score = 0;
				Solution::Position best_hpos, best_tpos;
				for(int ii = 0; ii < hlen; ++ii){
					for(int jj = 0; jj < tlen; ++jj){
						const int hpos = std::get<1>(hranges[i]) + ii;
						const int tpos = std::get<1>(tranges[k]) + jj;
						const int sw_score =
							head_scores[ii].first + tail_scores[jj].first;
						const int interval =
							tpos - (hpos - ii + head_scores[ii].second);
						const double score = 1.0 / compute_confidence(sw_score, interval);
						if(score >= best_score){
							best_score = score;
							best_hpos = Solution::Position(
								hpos, hpos - ii + head_scores[ii].second, false);
							best_tpos = Solution::Position(
								tpos, tpos - jj + tail_scores[jj].second, true);
						}
					}
				}
				if(pq.size() < 32 || pq.top() < best_score){
					pq.push(best_score);
					if(pq.size() > 32){ pq.pop(); }
				}
				if(solution.confidence() < best_score){
					const auto hpos = best_hpos;
					const auto tpos = best_tpos;
					solution.chromatid_id(std::get<0>(hranges[i]));
					if(!is_strand){
						solution.positions(0) = hpos;
						solution.positions(1) = tpos;
					}else{
						solution.positions(0) = tpos;
						solution.positions(1) = hpos;
					}
					solution.confidence(best_score);
				}
				++k;
			}
		}
	}

	Solution solve(const Query& query) const {
		Solution solution;
		std::priority_queue<double, std::vector<double>, std::greater<double>> pq;
		solve_pair(solution, pq, query.fragments(0), query.rev_fragments(1), false);
		solve_pair(solution, pq, query.fragments(1), query.rev_fragments(0), true);
		double sum = 0;
		while(!pq.empty()){
			sum += pq.top();
			pq.pop();
		}
		solution.confidence(solution.confidence() / sum);
		//std::cout << std::endl;
		return solution;
	}

public:
	int initTest(int){
		return 0;
	}

	int passReferenceGenome(
		int chromatid_sequence_id,
		std::vector<std::string> chromatid_sequence)
	{
		std::size_t total_length = 0;
		for(const auto& s: chromatid_sequence){ total_length += s.size(); }
		std::vector<uint8_t> combined;
		combined.reserve(total_length);
		for(const auto& s : chromatid_sequence){
			for(const auto c : s){
				switch(c){
					case 'A': combined.push_back(0); break;
					case 'T': combined.push_back(1); break;
					case 'C': combined.push_back(2); break;
					case 'G': combined.push_back(3); break;
					default: combined.push_back(INVALID_BASE);
				}
			}
		}
		m_chromatids[chromatid_sequence_id] = std::move(combined);
		build_block_table(
			chromatid_sequence_id, m_chromatids[chromatid_sequence_id]);
		return 0;
	}

	int preProcessing(){
		return 0;
	}

	std::vector<std::string> getAlignment(
		int n,
		double norm_a,
		double norm_s,
		std::vector<std::string> read_name,
		std::vector<std::string> read_sequence)
	{
		std::vector<Solution> result(n / 2);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
		for(int i = 0; i < n; i += 2){
#ifndef RUN_ONE_BY_ONE
			if(i % 100 == 0){ std::cout << i << " / " << n << std::endl; }
#endif
			const Query q(read_sequence[i], read_sequence[i + 1]);
			result[i / 2] = solve(q);
		}

		std::vector<std::string> stringified(n);
		for(int i = 0; i < n / 2; ++i){
			stringified[i * 2 + 0] = result[i].to_string(i, 0);
			stringified[i * 2 + 1] = result[i].to_string(i, 1);
		}
		return stringified;
	}

};
