#include "submission_wrapper.hpp"
#include "dna_sequencing.hpp"

SubmissionWrapper::SubmissionWrapper()
	: m_impl(new DNASequencing())
{ }

SubmissionWrapper::SubmissionWrapper(SubmissionWrapper&&) noexcept = default;

SubmissionWrapper::~SubmissionWrapper() = default;

SubmissionWrapper&
SubmissionWrapper::operator=(SubmissionWrapper&&) noexcept = default;


int SubmissionWrapper::initTest(int test_difficulty){
	return m_impl->initTest(test_difficulty);
}

int SubmissionWrapper::passReferenceGenome(
	int chromatid_sequence_id,
	std::vector<std::string> chromatid_sequence)
{
	return m_impl->passReferenceGenome(
		chromatid_sequence_id, std::move(chromatid_sequence));
}

int SubmissionWrapper::preProcessing(){
	return m_impl->preProcessing();
}

std::vector<std::string> SubmissionWrapper::getAlignment(
	int n,
	double norm_a,
	double norm_s,
	std::vector<std::string> read_name,
	std::vector<std::string> read_sequence)
{
	return m_impl->getAlignment(
		n, norm_a, norm_s, std::move(read_name), std::move(read_sequence));
}
