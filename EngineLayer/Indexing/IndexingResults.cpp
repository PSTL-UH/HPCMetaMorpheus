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

		std::wstring IndexingResults::ToString()
		{
			auto sb = new StringBuilder();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			sb->appendLine(MetaMorpheusEngineResults::ToString());
			sb->appendLine(L"\t\tfragmentIndexDict.Count: " + std::to_wstring(getFragmentIndex().size()));
			if (getPrecursorIndex().size() > 0)
			{
				sb->appendLine(L"\t\tprecursorIndexDict.Count: " + std::to_wstring(getPrecursorIndex().size()));
			}
			sb->appendLine(L"\t\tpeptideIndex.Count: " + std::to_wstring(getPeptideIndex().size()));

			delete sb;
			return sb->toString();
		}
	}
}
