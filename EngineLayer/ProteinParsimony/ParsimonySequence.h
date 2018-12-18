#pragma once

#include <string>
#include <any>

using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
	namespace ProteinParsimony
	{
		class ParsimonySequence
		{
		private:
			std::wstring privateSequence;
			Protease *privateProtease;

		public:
			ParsimonySequence(PeptideWithSetModifications *pwsm, bool TreatModPeptidesAsDifferentPeptides);

				std::wstring getSequence() const;
				Protease *getProtease() const;

			bool Equals(std::any obj) override;

			int GetHashCode() override;
		};
	}
}
