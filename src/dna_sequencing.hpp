#pragma once

#include <iostream>
#include "query.hpp"
#include "solution.hpp"

class CombinedChromatids {

public:
	static const uint32_t A = 0u;
	static const uint32_t T = 1u;
	static const uint32_t C = 2u;
	static const uint32_t G = 3u;
	static const uint32_t N = 15u;

private:
	static const uint32_t PADDING_LENGTH = 1024u;
	std::vector<uint32_t> m_sequence;
	std::vector<std::pair<uint32_t, int>> m_head_pointers;

public:
	CombinedChromatids()
		: m_sequence()
		, m_head_pointers()
	{ }

	template <typename Iterator>
	CombinedChromatids(Iterator begin, Iterator end)
		: m_sequence()
		, m_head_pointers()
	{
		uint32_t total_length = 0;
		for(auto it = begin; it != end; ++it){
			m_head_pointers.emplace_back(total_length, it->first);
			for(const auto& s : it->second){ total_length += s.size(); }
			total_length += PADDING_LENGTH;
			total_length += 32u - (total_length & 31u);
		}
		uint32_t written = 0u;
		m_sequence = std::vector<uint32_t>(total_length >> 3, 0xffffffffu);
		for(auto it = begin; it != end; ++it){
			for(const auto& s : it->second){
				for(const char c : s){
					const auto index = written >> 3;
					const auto shift = (written & 7) * 4;
					const auto mask = ~(0x0fu << shift);
					switch(c){
						case 'A':
							m_sequence[index] &= mask;
							m_sequence[index] |= (A << shift);
							break;
						case 'T':
							m_sequence[index] &= mask;
							m_sequence[index] |= (T << shift);
							break;
						case 'C':
							m_sequence[index] &= mask;
							m_sequence[index] |= (C << shift);
							break;
						case 'G':
							m_sequence[index] &= mask;
							m_sequence[index] |= (G << shift);
							break;
					}
					++written;
				}
			}
			written += PADDING_LENGTH;
			written += 32u - (written & 31u);
		}
	}

	size_t size() const {
		return m_sequence.size() * 8;
	}

	uint32_t operator[](uint32_t i) const {
		return (m_sequence[i >> 3] >> ((i & 7) * 4)) & 0x0fu;
	}

	std::pair<int, uint32_t> translate(uint32_t i) const {
		auto it = std::lower_bound(
			m_head_pointers.begin(), m_head_pointers.end(),
			std::make_pair(i + 1, 0));
		--it;
		return std::make_pair(it->second, i - it->first);
	}

	std::vector<uint8_t> subchromatid(uint32_t p, uint32_t len) const {
		std::vector<uint8_t> v(len);
		for(uint32_t i = 0; i < len; ++i){ v[i] = (*this)[p + i]; }
		return v;
	}

};

class ChromatidBlockMap {

public:
	using ValueIterator = std::vector<uint32_t>::const_iterator;
	using ValueRange = std::pair<ValueIterator, ValueIterator>;

	static const int BLOCK_SIZE = 25;
	static const int BLOCK_STEP = 13;

private:
	static const int JUMP_TABLE_DEPTH = 30;
	std::vector<uint32_t> m_key_lower_table;
	std::vector<uint32_t> m_range_offsets;
	std::vector<uint32_t> m_positions;
	std::vector<uint32_t> m_jump_table;

	size_t select_pivot(const uint64_t *keys, size_t n) const {
		const auto a = keys[0], b = keys[n / 2], c = keys[n - 1];
		const auto median = a + b + c - std::min({ a, b, c }) - std::max({ a, b, c });
		if(a == median){ return 0; }
		if(b == median){ return n / 2; }
		return n - 1;
	}

	void key_value_sort(
		uint64_t *keys, uint32_t *values, size_t n) const
	{
		if(n <= 1){ return; }
		if(n == 2){
			if(keys[0] > keys[1]){
				std::swap(  keys[0],   keys[1]);
				std::swap(values[0], values[1]);
			}
			return;
		}
		const auto pivot_index = select_pivot(keys, n);
		std::swap(  keys[pivot_index],   keys[n - 1]);
		std::swap(values[pivot_index], values[n - 1]);
		const auto x = keys[n - 1];
		size_t i = static_cast<size_t>(-1);
		for(size_t j = 0; j < n - 1; ++j){
			if(keys[j] <= x){
				++i;
				std::swap(  keys[i],   keys[j]);
				std::swap(values[i], values[j]);
			}
		}
		++i;
		std::swap(  keys[i],   keys[n - 1]);
		std::swap(values[i], values[n - 1]);
		key_value_sort(keys, values, i);
		key_value_sort(keys + i + 1, values + i + 1, n - (i + 1));
	}

public:
	ChromatidBlockMap()
		: m_key_lower_table()
		, m_range_offsets()
		, m_positions()
		, m_jump_table()
	{ }

	ChromatidBlockMap(const CombinedChromatids& cc)
		: m_key_lower_table()
		, m_range_offsets()
		, m_positions()
		, m_jump_table()
	{
		const uint64_t hmask = (4ull << (BLOCK_SIZE * 2)) - 1;
		const uint32_t n = cc.size();
		std::vector<uint64_t> key_work;
		key_work.reserve(n * 2 / BLOCK_STEP);
		m_positions.reserve(n * 2 / BLOCK_STEP);
		for(uint32_t i = 0; i + BLOCK_SIZE + 1 <= n; i += BLOCK_STEP){
			bool is_valid = true;
			uint64_t h = 0;
			for(uint32_t j = 0; is_valid && j <= BLOCK_SIZE; ++j){
				const auto x = cc[i + j];
				if(x == CombinedChromatids::N){
					is_valid = false;
				}else{
					h = (h << 2) | x;
				}
			}
			if(is_valid){
				key_work.push_back(h >> 2);
				m_positions.push_back(i);
				key_work.push_back(h & (hmask >> 2));
				m_positions.push_back(i);
			}
		}
		key_work.shrink_to_fit();
		m_positions.shrink_to_fit();
		key_value_sort(key_work.data(), m_positions.data(), key_work.size());

		size_t tail = 0;
		for(size_t i = 0; i < key_work.size(); ){
			size_t j = i;
			while(j < key_work.size() && key_work[i] == key_work[j]){ ++j; }
			if(j - i < 2500){
				for(uint32_t k = i; k < j; ++k){
					key_work[tail]    = key_work[k];
					m_positions[tail] = m_positions[k];
					++tail;
				}
			}
			i = j;
		}
		key_work.resize(tail);
		key_work.shrink_to_fit();
		m_positions.resize(tail);
		m_positions.shrink_to_fit();

		const int key_shift = (BLOCK_SIZE * 2 - JUMP_TABLE_DEPTH);
		const uint32_t lower_mask = (1u << key_shift) - 1u;
		const uint32_t m = key_work.size();
		uint32_t num_keys = 0;
		for(uint32_t i = 0; i < m; ++i){
			if(i == 0 || key_work[i] != key_work[i - 1]){ ++num_keys; }
		}
		m_key_lower_table = std::vector<uint32_t>(num_keys);
		m_range_offsets = std::vector<uint32_t>(num_keys + 1);
		m_jump_table = std::vector<uint32_t>((1 << JUMP_TABLE_DEPTH) + 1, 0xffffffffu);
		for(uint32_t i = 0, j = 0; i < m; ++i){
			const auto key = key_work[i];
			if(i == 0 || key != key_work[i - 1]){
				if(i == 0 || (key >> key_shift) != (key_work[i - 1] >> key_shift)){
					m_jump_table[key >> key_shift] = j;
				}
				m_key_lower_table[j] = static_cast<uint32_t>(key_work[i] & lower_mask);
				m_range_offsets[j] = i;
				++j;
			}
		}
		m_range_offsets.back() = m;
		m_jump_table.back() = num_keys;
		for(uint32_t i = (1 << JUMP_TABLE_DEPTH) - 1; i > 0; --i){
			m_jump_table[i] = std::min(m_jump_table[i + 1], m_jump_table[i]);
		}
		m_jump_table[0] = 0;
	}

	ValueRange equal_range(uint64_t key) const {
		const int shift = (BLOCK_SIZE * 2) - JUMP_TABLE_DEPTH;
		const auto lo = m_jump_table[key >> shift];
		const auto hi = m_jump_table[(key >> shift) + 1];
		const uint32_t key_lower = static_cast<uint32_t>(key & ((1u << shift) - 1));
		const auto begin = m_key_lower_table.begin() + lo;
		const auto end = m_key_lower_table.begin() + hi;
		const auto it = std::lower_bound(begin, end, key_lower);
		if(it == end || *it != key_lower){
			return std::make_pair(m_positions.begin(), m_positions.begin());
		}
		const auto idx = it - m_key_lower_table.begin();
		return std::make_pair(
			m_positions.begin() + m_range_offsets[idx],
			m_positions.begin() + m_range_offsets[idx + 1]);
	}

};

class DNASequencing {

public:
	struct CandidatePair {
		double weight;
		uint32_t head_position;
		uint32_t tail_position;
		bool strand;

		CandidatePair()
			: weight(0.0)
			, head_position(0)
			, tail_position(0)
			, strand(false)
		{ }

		CandidatePair(double w, uint32_t h, uint32_t t, bool st)
			: weight(w)
			, head_position(h)
			, tail_position(t)
			, strand(st)
		{ }
	};

private:
	std::vector<std::pair<int, std::vector<std::string>>> m_raw_chromatids;

	CombinedChromatids m_chromatids;
	ChromatidBlockMap m_block_map;

	template <typename T, int DEPTH = sizeof(T)>
	void radix_sort(T *data, size_t n) const {
		if(n < 256){
			std::sort(data, data + n);
			return;
		}
		const size_t CACHE_BLOCK_SIZE = 64;
		std::vector<T> work_buffer(n);
		T cache[256][CACHE_BLOCK_SIZE];
		size_t cached_count[256];
		size_t offset[257];
		T *cur = data, *next = work_buffer.data();
		for(int d = 0; d < DEPTH; ++d){
			const int shift = d * 8;
			for(size_t i = 0; i <= 256; ++i){ offset[i] = 0; }
			for(size_t i = 0; i < n; ++i){
				++offset[static_cast<int>((cur[i] >> shift) & 0xff) + 1];
			}
			for(size_t i = 0; i < 256; ++i){
				offset[i + 1] += offset[i];
				cached_count[i] = 0;
			}
			for(size_t i = 0; i < n; ++i){
				const auto g = static_cast<int>((cur[i] >> shift) & 0xff);
				cache[g][cached_count[g]++] = cur[i];
				if(cached_count[g] == CACHE_BLOCK_SIZE){
					for(size_t j = 0; j < CACHE_BLOCK_SIZE; ++j){
						next[offset[g]++] = cache[g][j];
					}
					cached_count[g] = 0;
				}
			}
			for(size_t g = 0; g < 256; ++g){
				for(size_t j = 0; j < cached_count[g]; ++j){
					next[offset[g]++] = cache[g][j];
				}
			}
			std::swap(cur, next);
		}
		if(DEPTH & 1){
			for(size_t i = 0; i < n; ++i){ next[i] = cur[i]; }
		}
	}

	std::vector<uint32_t> enumerate_block_matches(const std::vector<uint8_t>& q) const {
		const int BLOCK_SIZE = ChromatidBlockMap::BLOCK_SIZE;
		const uint64_t hmask = (1ull << (BLOCK_SIZE * 2)) - 1;
		const int n = q.size();
		uint64_t h = 0;
		for(int i = 0; i + 1 < BLOCK_SIZE; ++i){ h = (h << 2) | q[i]; }
		std::vector<uint32_t> result;
		for(int i = 0; i + BLOCK_SIZE <= n; ++i){
			const auto j = i + BLOCK_SIZE - 1;
			h = ((h << 2) | q[j]) & hmask;
			if(i % 2 == 0){
				const auto range = m_block_map.equal_range(h);
				for(auto it = range.first; it != range.second; ++it){
					result.push_back(*it - i);
				}
			}
		}
		return result;
	}

	std::vector<std::pair<uint32_t, int>> merge_block_matches(
		const std::vector<uint32_t>& matches) const
	{
		static const uint32_t MERGE_THRESHOLD = 30;
		const int n = matches.size();
		std::vector<std::pair<uint32_t, int>> result;
		for(int i = 0; i < n; ){
			const uint32_t merge_limit = matches[i] + MERGE_THRESHOLD;
			int j = i;
			while(j < n && matches[j] < merge_limit){ ++j; }
			result.emplace_back(matches[i], j - i);
			i = j;
		}
		return result;
	}

	void prune_and_sort_block_matches(
		std::vector<uint32_t>& hmatches,
		std::vector<uint32_t>& tmatches) const
	{
		const int THRESHOLD = 20;
		const int MIN_INTERVAL =  150;
		const int MAX_INTERVAL = 1050;
		const int n = hmatches.size(), m = tmatches.size();
		if(n <= THRESHOLD && m <= THRESHOLD){
			radix_sort(hmatches.data(), n);
			radix_sort(tmatches.data(), m);
		}else if(n <= THRESHOLD){
			radix_sort(hmatches.data(), n);
			int tail = m - 1;
			for(int i = 0; i <= tail; ++i){
				const uint32_t lo = tmatches[i] - MAX_INTERVAL;
				const uint32_t hi = tmatches[i] - MIN_INTERVAL;
				bool found = false;
				for(const auto q : hmatches){
					if(lo <= q && q < hi){
						found = true;
						break;
					}
				}
				if(!found){
					std::swap(tmatches[i], tmatches[tail]);
					--i;
					--tail;
				}
			}
			tmatches.resize(tail + 1);
			radix_sort(tmatches.data(), tmatches.size());
		}else if(m <= THRESHOLD){
			radix_sort(tmatches.data(), m);
			int tail = n - 1;
			for(int i = 0; i <= tail; ++i){
				const uint32_t lo = hmatches[i] + MIN_INTERVAL;
				const uint32_t hi = hmatches[i] + MAX_INTERVAL;
				bool found = false;
				for(const auto q : tmatches){
					if(lo <= q && q < hi){
						found = true;
						break;
					}
				}
				if(!found){
					std::swap(hmatches[i], hmatches[tail]);
					--i;
					--tail;
				}
			}
			hmatches.resize(tail + 1);
			radix_sort(hmatches.data(), hmatches.size());
		}else{
			radix_sort(hmatches.data(), n);
			radix_sort(tmatches.data(), m);
		}
	}

	std::vector<CandidatePair> enumerate_candidates(
		const std::vector<uint8_t>& head,
		const std::vector<uint8_t>& tail,
		bool strand) const
	{
		const uint32_t MIN_BLOCK_INTERVAL = 200u;
		const uint32_t MAX_BLOCK_INTERVAL = 1000u;

		auto hmatches = enumerate_block_matches(head);
		if(hmatches.empty()){ return std::vector<CandidatePair>(); }
		auto tmatches = enumerate_block_matches(tail);
		if(tmatches.empty()){ return std::vector<CandidatePair>(); }

		prune_and_sort_block_matches(hmatches, tmatches);
		if(hmatches.empty() || tmatches.empty()){
			return std::vector<CandidatePair>();
		}
		const auto hblocks = merge_block_matches(hmatches);
		const auto tblocks = merge_block_matches(tmatches);
		const int n = hblocks.size(), m = tblocks.size();
		std::vector<CandidatePair> result;
		for(int i = 0, j = 0; i < n; ++i){
			while(j < m && tblocks[j].first < hblocks[i].first + MIN_BLOCK_INTERVAL){ ++j; }
			for(int k = j; k < m; ++k){
				const auto& head = hblocks[i];
				const auto& tail = tblocks[k];
				if(head.first + MAX_BLOCK_INTERVAL <= tail.first){ break; }
				result.emplace_back(
					head.second * tail.second, head.first, tail.first, strand);
			}
		}
		return result;
	}

	int levenstein(const std::vector<uint8_t>& a, const std::vector<uint8_t>& b) const {
		const int INF = 1000000000;
		const int n = a.size(), m = b.size();
		std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, INF));
		for(int i = 0; i <= n; ++i){ dp[i][0] = 0; }
		for(int j = 1; j <= m; ++j){ dp[0][j] = j; }
		for(int i = 1; i <= n; ++i){
			for(int j = 1; j <= m; ++j){
				dp[i][j] = std::min(
					std::min(dp[i - 1][j], dp[i][j - 1]) + 1,
					dp[i - 1][j - 1] + (a[i - 1] == b[j - 1] ? 0 : 1));
			}
		}
		int score = INF;
		for(int i = 0; i <= n; ++i){ score = std::min(score, dp[i][m]); }
		return score;
	}

	double compute_fine_confidence(
		const CandidatePair& candidate,
		const std::vector<uint8_t>& head,
		const std::vector<uint8_t>& tail) const
	{
		const uint32_t HEAD_PADDING = 20u, SUBCHROMATID_LENGTH = 200u;
		const auto hpos = candidate.head_position;
		const auto tpos = candidate.tail_position;
		const uint32_t hh = (hpos > HEAD_PADDING ? hpos - HEAD_PADDING : 0);
		const uint32_t ht = std::min<uint32_t>(hh + SUBCHROMATID_LENGTH, m_chromatids.size());
		const uint32_t th = (tpos > HEAD_PADDING ? tpos - HEAD_PADDING : 0);
		const uint32_t tt = std::min<uint32_t>(th + SUBCHROMATID_LENGTH, m_chromatids.size());
		const auto sh = levenstein(m_chromatids.subchromatid(hh, ht - hh), head);
		const auto st = levenstein(m_chromatids.subchromatid(th, tt - th), tail);
		const int align = sh + st, interval = tpos - hpos;
		const double n_align = align, n_interval = (interval - 450) * 0.005;
		const double confidence = 1.0 / (1.0 + sqrt(n_align * n_align + n_interval * n_interval));
		return confidence;
	}

	Solution solve(const Query& q){
		const int TENTATIVE_LENGTH = 150;
		const double WEIGHT_THRESHOLD_MULTIPLY = 0.3;
		const double COARSE_LO_THRESHOLD = 0.1;
		const double COARSE_HI_THRESHOLD = 0.7;

		const auto c0 = enumerate_candidates(q.fragments(0), q.rev_fragments(1), false);
		const auto c1 = enumerate_candidates(q.fragments(1), q.rev_fragments(0), true);
		if(c0.empty() && c1.empty()){ return Solution(); }

		double max_weight = 0.0;
		for(const auto& x : c0){ max_weight = std::max(max_weight, x.weight); }
		for(const auto& x : c1){ max_weight = std::max(max_weight, x.weight); }
		const double weight_threshold = max_weight * WEIGHT_THRESHOLD_MULTIPLY;

		std::vector<CandidatePair> merged;
		for(const auto& x : c0){
			if(x.weight >= weight_threshold){ merged.push_back(x); }
		}
		for(const auto& x : c1){
			if(x.weight >= weight_threshold){ merged.push_back(x); }
		}
		std::sort(
			merged.begin(), merged.end(),
			[](const CandidatePair& a, const CandidatePair& b) -> bool {
				return a.weight > b.weight;
			});
		double weight_sum = 0.0;
		for(const auto& x : merged){
			weight_sum += x.weight;
		}

		Solution solution;
		const auto x = merged.front();
		const auto c_confidence = x.weight / weight_sum;
		if(COARSE_LO_THRESHOLD < c_confidence && c_confidence < COARSE_HI_THRESHOLD){
			// Use fine (slow) algorithm
			const int n = merged.size();
			std::vector<double> fine_confidences(n);
			double fc_sum = 0.0;
			for(int i = 0; i < n; ++i){
				if(!merged[i].strand){
					fine_confidences[i] = compute_fine_confidence(
						merged[i], q.fragments(0), q.rev_fragments(1));
				}else{
					fine_confidences[i] = compute_fine_confidence(
						merged[i], q.fragments(1), q.rev_fragments(0));
				}
				fc_sum += fine_confidences[i];
			}
			const auto argmax =
				std::max_element(fine_confidences.begin(), fine_confidences.end());
			const double f_confidence = *argmax / fc_sum;
			const auto& y = merged[argmax - fine_confidences.begin()];
			const auto hpos = m_chromatids.translate(y.head_position);
			const auto tpos = m_chromatids.translate(y.tail_position);
			solution.chromatid_id(hpos.first);
			if(!y.strand){
				solution.positions(0) = Solution::Position(
					hpos.second, hpos.second + TENTATIVE_LENGTH, false);
				solution.positions(1) = Solution::Position(
					tpos.second, tpos.second + TENTATIVE_LENGTH, true);
			}else{
				solution.positions(0) = Solution::Position(
					tpos.second, tpos.second + TENTATIVE_LENGTH, true);
				solution.positions(1) = Solution::Position(
					hpos.second, hpos.second + TENTATIVE_LENGTH, false);
			}
			solution.confidence(
				f_confidence * (COARSE_HI_THRESHOLD - COARSE_LO_THRESHOLD) + COARSE_LO_THRESHOLD);
		}else{
			// Use result of coarse algorithm
			const auto hpos = m_chromatids.translate(x.head_position);
			const auto tpos = m_chromatids.translate(x.tail_position);
			solution.chromatid_id(hpos.first);
			if(!x.strand){
				solution.positions(0) = Solution::Position(
					hpos.second, hpos.second + TENTATIVE_LENGTH, false);
				solution.positions(1) = Solution::Position(
					tpos.second, tpos.second + TENTATIVE_LENGTH, true);
			}else{
				solution.positions(0) = Solution::Position(
					tpos.second, tpos.second + TENTATIVE_LENGTH, true);
				solution.positions(1) = Solution::Position(
					hpos.second, hpos.second + TENTATIVE_LENGTH, false);
			}
			solution.confidence(c_confidence);
		}
		return solution;
	}

public:
	DNASequencing()
		: m_raw_chromatids()
		, m_chromatids()
		, m_block_map()
	{ }

	int initTest(int){
		return 0;
	}

	int passReferenceGenome(
		int chromatid_sequence_id,
		std::vector<std::string> chromatid_sequence)
	{
		m_raw_chromatids.emplace_back(
			chromatid_sequence_id, std::move(chromatid_sequence));
		return 0;
	}

	int preProcessing(){
		m_chromatids = CombinedChromatids(
			m_raw_chromatids.begin(), m_raw_chromatids.end());
		m_raw_chromatids.clear();
		m_block_map = ChromatidBlockMap(m_chromatids);
		return 0;
	}

	std::vector<std::string> getAlignment(
		int n, double norm_a, double norm_s,
		std::vector<std::string> read_name,
		std::vector<std::string> read_sequence)
	{
		std::vector<Solution> result(n / 2);
#pragma omp parallel for schedule(dynamic)
		for(int i = 0; i < n; i += 2){
#ifndef RUN_ONE_BY_ONE
			if(i % 10000 == 0){ std::cout << i << " / " << n << std::endl; }
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

