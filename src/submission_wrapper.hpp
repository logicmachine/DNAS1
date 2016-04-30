#pragma once

#include <vector>
#include <string>
#include <memory>

class DNASequencing;

class SubmissionWrapper {

private:
	std::unique_ptr<DNASequencing> m_impl;

public:
	SubmissionWrapper();
	SubmissionWrapper(SubmissionWrapper&&) noexcept;
	~SubmissionWrapper();

	SubmissionWrapper &operator=(SubmissionWrapper&&) noexcept;

	int initTest(int test_difficuluty);

	int passReferenceGenome(
		int chromatid_sequence_id,
		std::vector<std::string> chromatid_sequence);

	int preProcessing();

	std::vector<std::string> getAlignment(
		int n,
		double norm_a,
		double norm_s,
		std::vector<std::string> read_name,
		std::vector<std::string> read_seqence);

};
