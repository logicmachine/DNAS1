#include <iostream>
#include <array>
#include <string>
#include <sstream>
#include <chrono>
#include "submission_wrapper.hpp"
#include "problem.hpp"
#include "solution.hpp"

static const int MAX_POSITION_DISTANCE = 300;
static const double MAX_AUC = 0.999999;

static Solution read_solution(std::string line1, std::string line2){
	for(char& c : line1){
		if(c == ','){ c = ' '; }
	}
	for(char& c : line2){
		if(c == ','){ c = ' '; }
	}
	std::istringstream iss1(line1), iss2(line2);
	std::array<std::istringstream *, 2> iss = { &iss1, &iss2 };
	std::string name[2], strand[2];
	int chromatid_id[2], start[2], end[2];
	double confidence[2];
	for(int i = 0; i < 2; ++i){
		*(iss[i]) >> name[i] >> chromatid_id[i];
		*(iss[i]) >> start[i] >> end[i] >> strand[i] >> confidence[i];
	}
	return Solution(
		chromatid_id[0], confidence[0],
		start[0] - 1, end[0], strand[0][0] == '-',
		start[1] - 1, end[1], strand[1][0] == '-');
}

static std::pair<bool, Solution> verify_solution(
	const Problem& problem,
	std::size_t read_id,
	const std::vector<std::string>& actual)
{
	const auto& ex = problem.expected_solutions()[read_id];
	const auto& ac = read_solution(actual[0], actual[1]);
	bool fail = false;
	for(int j = 0; j < 2; ++j){
		const auto e = ex.positions(j);
		const auto a = ac.positions(j);
		const int d = abs(static_cast<int>(e.start - a.start));
		if(e.is_reversed != a.is_reversed || d >= MAX_POSITION_DISTANCE){
			fail = true;
		}
	}
	return std::make_pair(fail, ac);
}

static double compute_auc(std::vector<std::pair<double, int>> results){
	const int n = results.size();
	std::sort(
		results.begin(), results.end(),
		std::greater<std::pair<double, int>>());
	std::vector<int> cumul_si = { results[0].second };
	std::vector<int> pos = { 0 };
	for(int i = 1; i < n; ++i){
		if(results[i].first == results[i - 1].first){
			cumul_si.back() += results[i].second;
			pos.back() = i;
		}else{
			const auto cumul = cumul_si.back() + results[i].second;
			cumul_si.push_back(cumul);
			pos.push_back(i);
		}
	}
	const int m = cumul_si.size();
	const double inv_n = 1.0 / n, inv_np1 = 1.0 / (n + 1);
	const double lf_multiplier = 1.0 / log(n + 1);
	double auc = 0.0;
	for(int i = 0; i < m; ++i){
		const double fi = 1.0 * (2 + pos[i] - cumul_si[i]) * inv_np1;
		const double fi1 = (i == m - 1) ?
			1.0 : 1.0 * (2 + pos[i + 1] - cumul_si[i + 1]) * inv_np1;
		const double lfi = lf_multiplier * log(fi);
		const double lfi1 = lf_multiplier * log(fi1);
		auc += cumul_si[i] * (lfi1 - lfi) * inv_n;
	}
	return auc;
}

static int print_scores(
	const Problem& problem,
	const std::vector<std::string>& actual,
	double preprocess_time,
	double align_time)
{
	std::cout << "## Incorrect Answers" << std::endl;
	const auto& expected = problem.expected_solutions();
	std::vector<std::pair<double, int>> verify_results;
	std::size_t correct_count = 0, confident_count = 0;
	for(std::size_t i = 0; i < expected.size(); ++i){
		const auto p = verify_solution(
			problem, i, { actual[i * 2], actual[i * 2 + 1] });
		const auto fail = p.first;
		const auto& ac = p.second;
		if(!fail){
			correct_count += 2;
			verify_results.emplace_back(ac.confidence(), 1);
			verify_results.emplace_back(ac.confidence(), 1);
		}else{
			verify_results.emplace_back(ac.confidence(), 0);
			verify_results.emplace_back(ac.confidence(), 0);
		}
	}
	std::cout << std::endl;
	std::sort(verify_results.begin(), verify_results.end());

	const double auc = compute_auc(verify_results);
	const double accuracy = log(1 - std::min(auc, MAX_AUC)) / problem.norm_a();

	const double speed =
		(1.0 - align_time / problem.time_cutoff()) / problem.norm_s();

	std::cout << "## Score Information" << std::endl;
	std::cout << "- Preprocessing time: "
	          << preprocess_time << " [s]" << std::endl;
	std::cout << "- Alignment time: "
              << align_time << " [s]" << std::endl;
	std::cout << "- Speed: " << speed << std::endl;
	std::cout << "- Correct answers: "
	          << correct_count << "/"
	          << expected.size() * 2
	          << " (" << confident_count << ")" << std::endl;
	std::cout << "- AuC: " << auc << std::endl;
	std::cout << "- Accuracy: " << accuracy << std::endl;
	std::cout << "- Score: " << accuracy * speed << std::endl;
	std::cout << std::endl;
}

int main(int argc, char *argv[]){
	if(argc < 2){
		std::cout << "Usage: " << argv[0] << " input" << std::endl;
		return 0;
	}
	std::cout << std::setiosflags(std::ios::fixed);
	const Problem problem(argv[1]);

	const auto preprocess_begin_time = std::chrono::steady_clock::now();
	SubmissionWrapper sequencer;
	sequencer.initTest(0);
	for(const auto& p : problem.chromatids()){
		sequencer.passReferenceGenome(p.first, p.second);
	}
	sequencer.preProcessing();
	const auto preprocess_end_time = std::chrono::steady_clock::now();
	const auto preprocess_time =
		std::chrono::duration_cast<std::chrono::milliseconds>(
			preprocess_end_time - preprocess_begin_time).count() * 1e-3;

#ifdef RUN_ONE_BY_ONE
	const auto n = problem.read_names().size();
	for(std::size_t i = 0; i < n; i += 2){
		const std::vector<std::string> read_names = {
			problem.read_names()[i + 0],
			problem.read_names()[i + 1]
		};
		const std::vector<std::string> read_sequences = {
			problem.read_sequences()[i + 0],
			problem.read_sequences()[i + 1]
		};
		const auto actual = sequencer.getAlignment(
			read_names.size(),
			problem.norm_a(),
			problem.norm_s(),
			read_names,
			read_sequences);
		const auto p = verify_solution(problem, i / 2, actual);
		if(p.first){
			const auto& ex = problem.expected_solutions()[i / 2];
			std::cout << "(" << read_names[0] << ", " << read_names[1] << ")" << std::endl;
			for(int j = 0; j < 2; ++j){
				const auto& pos = ex.positions(j);
				std::cout << "  + " << pos.start << ", " << pos.end << ", " << pos.is_reversed << std::endl;
				std::cout << "  - " << actual[j] << std::endl;
			}
		}
	}
#else
	const auto align_begin_time = std::chrono::steady_clock::now();
	const auto actual = sequencer.getAlignment(
		problem.read_names().size(),
		problem.norm_a(),
		problem.norm_s(),
		problem.read_names(),
		problem.read_sequences());
	const auto align_end_time = std::chrono::steady_clock::now();
	const auto align_time =
		std::chrono::duration_cast<std::chrono::milliseconds>(
			align_end_time - align_begin_time).count() * 1e-3;

	problem.print_information(std::cout);
	print_scores(problem, actual, preprocess_time, align_time);
#endif
	return 0;
}

