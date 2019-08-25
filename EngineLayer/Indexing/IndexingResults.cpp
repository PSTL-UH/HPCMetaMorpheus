#include "IndexingResults.h"
#include "IndexingEngine.h"

using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;
namespace EngineLayer
{
	namespace Indexing
	{

		IndexingResults::IndexingResults(std::vector<PeptideWithSetModifications*> &peptideIndex, std::vector<std::vector<int>&> &fragmentIndex, std::vector<std::vector<int>&> &precursorIndex, IndexingEngine *indexParams) : MetaMorpheusEngineResults(indexParams)
		{
			setPeptideIndex(peptideIndex);
			setFragmentIndex(fragmentIndex);
			setPrecursorIndex(precursorIndex);
		}

		std::vector<std::vector<int>> IndexingResults::getFragmentIndex() const
		{
			return privateFragmentIndex;
		}

		void IndexingResults::setFragmentIndex(const std::vector<std::vector<int>> &value)
		{
			privateFragmentIndex = value;
		}

		std::vector<std::vector<int>> IndexingResults::getPrecursorIndex() const
		{
			return privatePrecursorIndex;
		}

		void IndexingResults::setPrecursorIndex(const std::vector<std::vector<int>> &value)
		{
			privatePrecursorIndex = value;
		}

		std::vector<PeptideWithSetModifications*> IndexingResults::getPeptideIndex() const
		{
			return privatePeptideIndex;
		}

		void IndexingResults::setPeptideIndex(const std::vector<PeptideWithSetModifications*> &value)
		{
			privatePeptideIndex = value;
		}

		std::string IndexingResults::ToString()
		{
			auto sb = new StringBuilder();

			sb->appendLine(MetaMorpheusEngineResults::ToString());
			sb->appendLine("\t\tfragmentIndexDict.Count: " + std::to_string(getFragmentIndex().size()));
			if (getPrecursorIndex().size() > 0)
			{
				sb->appendLine("\t\tprecursorIndexDict.Count: " + std::to_string(getPrecursorIndex().size()));
			}
			sb->appendLine("\t\tpeptideIndex.Count: " + std::to_string(getPeptideIndex().size()));

			delete sb;
			return sb->toString();
		}
	}
}
